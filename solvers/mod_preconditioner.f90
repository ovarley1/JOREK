module mod_preconditioner

  private
  public initialize_preconditioner, reset_preconditioner, update_pc_rhs, gather_solution

  contains

!> Initialize preconditioner structure (PC)
!! call subroutines for setting mode families, creating communicators, distributing tasks
  subroutine initialize_preconditioner(pc,comm_glob)
    use phys_module, only: autodistribute_modes, n_mode_families, autodistribute_ranks, centralize_harm_mat
    use phys_module, only: freeze_harmonics, linear_pc_evolving
    use mod_parameters, only: n_tor
    use mod_linear_compact, only: n_act
    use data_structure, only: type_PRECOND
    use mpi_mod
    implicit none

    type(type_PRECOND) :: pc
    integer            :: comm_glob, my_id, n_cpu, ierr
    integer            :: i
    character(len=256) :: s

    pc%comm = comm_glob

    call MPI_COMM_RANK(pc%comm, my_id, ierr)
    call MPI_COMM_SIZE(pc%comm, n_cpu, ierr)

    pc%my_id = my_id
    pc%n_cpu = n_cpu
    if (pc%my_id.eq.0) write(*,*) "Initializing preconditioner"

    pc%autodistribute_ranks = autodistribute_ranks
    pc%autodistribute_modes = autodistribute_modes
    !pc%mat%row_distributed  = .not.centralize_harm_mat

    if (pc%autodistribute_modes) then
      pc%n_mode_families = (n_act + 1)/2 ! n_act is n_tor, or the number of active harmonics with linear_compact
    else
      pc%n_mode_families = n_mode_families
    endif

    ! --- Linear run by mode family: the families come from the evolving harmonics only, frozen ones get no PC rows
#ifndef DIRECT_CONSTRUCTION
    if (freeze_harmonics .and. linear_pc_evolving) call set_evolving_families(pc)
#else
    if (freeze_harmonics .and. linear_pc_evolving .and. (my_id.eq.0)) then
      write(*,*) "linear_pc_evolving is not supported with USE_DIRECT_CONSTRUCTION=1: ignored"
    endif
#endif

    call distribute_ranks(n_cpu, pc)

    call create_communicators(pc)
    pc%mat%comm = pc%MPI_COMM_N ! communicator for PC matrix distribution

    call distribute_modes(pc)

    pc%initialized = .true.

    if (my_id.eq.0) then
      do i=1, pc%n_mode_families
        write(s,'(A17,i4,A12)') " mode_family_id: ", i, " MPI ranks: "
        write(*,*) trim(s), pc%mode_families_ranks(i,1:pc%ranks_per_family(i))
      enddo
      do i=1, pc%n_mode_families
        write(s,'(A17,i4,A9,f6.2)') " mode_family_id: ", i, " weight: ", pc%row_factor
        write(*,*) trim(s), " modes:", pc%mode_families_modes(i,1:pc%modes_per_family(i))
      enddo
    endif

    return

  end subroutine initialize_preconditioner

#ifndef DIRECT_CONSTRUCTION
!> Linear run by mode family: one PC family per evolving cos/sin pair (i_tor = 2k, 2k+1, same n) and none for the frozen
!! harmonics. A harmonic in no family has no PC rows, so GMRES keeps it at exactly 0 (see gather_solution and prec in
!! mod_gmres2) and the PC is built and factorised for the evolving harmonics only. Sets the family arrays that
!! distribute_modes reads, and pc%n_mode_families as a local count (the input n_mode_families is left alone).
!! Runs on every rank, from harmonic_frozen and n_tor, which are identical on every rank (set by setup_linear_family).
!! With manual ranks, ranks_per_family(1:nf) must sum to the number of MPI ranks.
  subroutine set_evolving_families(pc)
    use phys_module, only: n_fam_max, harmonic_frozen, autodistribute_modes, autodistribute_ranks, ranks_per_family
    use phys_module, only: modes_per_family, mode_families_modes, weights_per_family
    use mod_parameters, only: n_tor
    use mod_linear_compact, only: harm_slot
    use data_structure, only: type_PRECOND
    use mpi_mod
    implicit none

    type(type_PRECOND) :: pc
    integer            :: k, nf, ierr

    ! --- Count the evolving pairs. n_tor is odd (n=0, then a cos/sin pair per n) and i_tor=1 (n=0) is always frozen.
    nf = 0
    do k = 1, (n_tor - 1)/2
      if (harmonic_frozen(2*k) .neqv. harmonic_frozen(2*k + 1)) then
        if (pc%my_id.eq.0) write(*,*) "FATAL: linear_pc_evolving: the cos/sin pair of harmonics is split by the freeze:", 2*k
        call MPI_Abort(MPI_COMM_WORLD, 6, ierr)
        stop
      endif
      if (.not. harmonic_frozen(2*k)) nf = nf + 1
    enddo

    if ((nf < 1) .or. (nf > n_fam_max)) then
      if (pc%my_id.eq.0) write(*,*) "FATAL: linear_pc_evolving: number of evolving families out of range:", nf, n_fam_max
      call MPI_Abort(MPI_COMM_WORLD, 6, ierr)
      stop
    endif

    ! --- Every family needs at least one rank. distribute_ranks would exit with status 0 on a manual mismatch
    if (pc%n_cpu < nf) then
      if (pc%my_id.eq.0) then
        write(*,*) "FATAL: linear_pc_evolving needs at least one MPI rank per family. Families, ranks:", nf, pc%n_cpu
      endif
      call MPI_Abort(MPI_COMM_WORLD, 6, ierr)
      stop
    endif
    if (.not. autodistribute_ranks) then
      if ((minval(ranks_per_family(1:nf)) < 1) .or. (sum(ranks_per_family(1:nf)) /= pc%n_cpu)) then
        if (pc%my_id.eq.0) then
          write(*,*) "FATAL: linear_pc_evolving with autodistribute_ranks=.f. needs ranks_per_family(1:nf) >= 1, summing to"
          write(*,*) "       the number of MPI ranks. Families, ranks, ranks_per_family:", nf, pc%n_cpu, ranks_per_family(1:nf)
        endif
        call MPI_Abort(MPI_COMM_WORLD, 6, ierr)
        stop
      endif
    endif

    ! --- Overwrite what the deck set: one family per evolving n, modes (2k, 2k+1), weight 1 (no overlap between families).
    ! --- The modes are the positions of the harmonics in the blocks of the global matrix: harm_slot is the identity, except with
    ! --- linear_compact where the frozen harmonics have no position and the evolving pair j sits in slots (2j-1, 2j).
    nf = 0
    do k = 1, (n_tor - 1)/2
      if (harmonic_frozen(2*k)) cycle
      nf = nf + 1
      modes_per_family(nf)        = 2
      mode_families_modes(nf,1:2) = (/ harm_slot(2*k), harm_slot(2*k + 1) /)
      weights_per_family(nf)      = 1.d0
    enddo

    if (pc%my_id.eq.0) then
      if (.not. autodistribute_modes) write(*,*) "linear_pc_evolving: the manual mode families of the deck are overridden"
      write(*,'(A,I0,A,I0,A)') " linear_pc_evolving: ", nf, " PC families, one per evolving n; the ", count(harmonic_frozen), &
                               " frozen harmonics get no PC rows"
    endif

    pc%autodistribute_modes = .false.
    pc%n_mode_families      = nf

  end subroutine set_evolving_families
#endif

!> Set up MPI communicators for mode families and corresponding masters
  subroutine create_communicators(pc)
    use data_structure, only: type_PRECOND
    use mpi_mod
    implicit none

    type(type_PRECOND) :: pc

    integer, allocatable :: i_tor(:), ranks_tmp(:)
    integer :: my_id, n_cpu, ierr

    my_id = pc%my_id
    n_cpu = pc%n_cpu

    call MPI_COMM_SPLIT(pc%comm, pc%family_id, my_id, pc%MPI_COMM_N, ierr)
    if (ierr.ne.0) then
      write(*,*) "Error in creating MPI_COMM_N"
      call MPI_Abort(MPI_COMM_WORLD, 0, ierr)
    endif
    call MPI_COMM_RANK(pc%MPI_COMM_N, pc%my_id_n, ierr)
    call MPI_COMM_SIZE(pc%MPI_COMM_N, pc%n_cpu_n, ierr)

    allocate(i_tor(n_cpu)); i_tor = 0
    i_tor(my_id + 1) = pc%my_id_n
    call MPI_Allreduce(MPI_IN_PLACE,i_tor,n_cpu,MPI_INT,MPI_SUM,pc%comm,ierr)
    call MPI_COMM_SPLIT(pc%comm,i_tor(my_id+1),my_id,pc%MPI_COMM_TRANS,ierr)

    pc%n_masters = pc%n_mode_families
    allocate(ranks_tmp(pc%n_masters)); ranks_tmp=0;

    if (pc%my_id_n.eq.0) ranks_tmp(pc%family_id) = my_id
    call MPI_AllReduce(MPI_IN_PLACE,ranks_tmp,pc%n_masters,MPI_INT,MPI_SUM,pc%comm,ierr)
    call MPI_COMM_GROUP(pc%comm,pc%MPI_GROUP_WORLD,ierr)
    call MPI_GROUP_INCL(pc%MPI_GROUP_WORLD,pc%n_masters,ranks_tmp,pc%MPI_GROUP_MASTER,ierr)
    call MPI_COMM_CREATE(pc%comm,pc%MPI_GROUP_MASTER,pc%MPI_COMM_MASTER,ierr)

    if (pc%my_id_n .eq. 0) then
     call MPI_COMM_RANK(pc%MPI_COMM_MASTER, pc%my_id_master, ierr)
    endif

    if ((my_id.eq.0).and.(pc%my_id_n.ne.0)) then
      write(*,*) "Error in creating communicators: my_id==0 must have my_id_n==0"
      call MPI_Abort(MPI_COMM_WORLD, 0, ierr)
    endif

    deallocate(i_tor, ranks_tmp)

    return

  end subroutine create_communicators

  !> Distribute toroidal modes among mode families
  subroutine distribute_modes(pc)
    use data_structure, only: type_PRECOND
    use phys_module, only: modes_per_family, mode_families_modes, weights_per_family
    implicit none

    type(type_PRECOND) :: pc
    integer            :: i, j, n_fam_max

    allocate(pc%modes_per_family(pc%n_mode_families))

    if (pc%autodistribute_modes) then
      pc%row_factor = 1.0
      pc%modes_per_family(1) = 1
      if (pc%n_mode_families>1) pc%modes_per_family(2:pc%n_mode_families) = 2
    else
      do i = 1, pc%n_mode_families
        pc%row_factor = weights_per_family(i)
        pc%modes_per_family(i) = modes_per_family(i)
      enddo
    endif

    n_fam_max = 1
    do i = 1, pc%n_mode_families
      n_fam_max = max(n_fam_max,pc%modes_per_family(i))
    enddo

    allocate(pc%mode_families_modes(pc%n_mode_families,n_fam_max))
    pc%mode_families_modes(:,:) = -1

    if (pc%autodistribute_modes) then
      pc%mode_families_modes(1,1) = 1
      if (pc%n_mode_families.gt.1) then
        do i = 2, pc%n_mode_families
          pc%mode_families_modes(i,1) =  (i - 1)*2
          pc%mode_families_modes(i,2) =  (i - 1)*2 + 1
        enddo
      endif
    else
      do i = 1, pc%n_mode_families
        do j = 1, modes_per_family(i)
          pc%mode_families_modes(i,j) = mode_families_modes(i,j)
        enddo
      enddo
    endif

    pc%mode_set_n = pc%modes_per_family(pc%family_id)
    allocate(pc%mode_set(pc%mode_set_n))
    pc%mode_set(1:pc%mode_set_n) = pc%mode_families_modes(pc%family_id,1:pc%mode_set_n)

  end subroutine distribute_modes

  !> Distribute MPI ranks among mode families
  subroutine distribute_ranks(n_cpu,pc)
    use data_structure, only: type_PRECOND
    use phys_module, only: ranks_per_family
    implicit none

    type(type_PRECOND) :: pc
    integer, intent(in)  :: n_cpu
    integer :: mcpu, r, i, j
    integer, dimension(:), pointer :: rank_id => Null()

    allocate(pc%rank_range(pc%n_mode_families + 1))
    allocate(pc%ranks_per_family(pc%n_mode_families))
    allocate(pc%mode_families_ranks(pc%n_mode_families,n_cpu))
    allocate(rank_id(n_cpu))

    do i = 1, pc%n_mode_families
      do j = 1, n_cpu
        pc%mode_families_ranks(i,j) = -1
      enddo
    enddo

    mcpu = n_cpu/pc%n_mode_families
    r = mod(n_cpu,pc%n_mode_families)

    if (pc%autodistribute_ranks) then
      do i=1, pc%n_mode_families
        pc%ranks_per_family(i) = mcpu
        if ((r.gt.0).and.(i.le.r))  pc%ranks_per_family(i) = pc%ranks_per_family(i) + 1 ! add extra rank if avaiable
      enddo
    else
      do i = 1, pc%n_mode_families
        pc%ranks_per_family(i) = ranks_per_family(i)
      enddo
    endif

    pc%rank_range(1) = 1
    do i = 2, pc%n_mode_families+1
      pc%rank_range(i) = pc%rank_range(i-1) + pc%ranks_per_family(i-1)
    enddo

    ! check for consistency
    r = 0
    do i= 2, pc%n_mode_families + 1
      r = r + pc%rank_range(i) - pc%rank_range(i-1)
    enddo
    if (r.ne.n_cpu) then
      write(*,*) "Error in distribution of ranks"
      call exit(0)
    endif

    do i = 1, n_cpu
      do j = 2, pc%n_mode_families + 1
        if ((i.ge.pc%rank_range(j-1)).and.(i.lt.pc%rank_range(j))) then
          rank_id(i) = j - 1
          exit
        endif
      enddo
    enddo

    do j=1,pc%n_mode_families
      r = 0
      do i = 1, n_cpu
        if (rank_id(i).eq.j) then
          r = r + 1
          pc%mode_families_ranks(j,r) = i - 1
        endif
      enddo
    enddo

    pc%family_id = rank_id(pc%my_id + 1)
    deallocate(rank_id)

    return

  end subroutine distribute_ranks
  
  !> Distribute RHS vector
  subroutine update_pc_rhs(pc,rhs_vec)
    use data_structure, only: type_PRECOND, type_RHS
    use mod_integer_types

    implicit none
    
    type(type_RHS)     :: rhs_vec
    type(type_PRECOND) :: pc

    integer(kind=int_all) :: i

    do i=1, pc%rhs%n
      pc%rhs%val(i) = rhs_vec%val(pc%row_index(i))
    enddo

  end subroutine update_pc_rhs
  
  !> Collect the solution vector from the individual mode groups
  subroutine gather_solution(pc,sol_vec)
    use data_structure, only: type_PRECOND, type_RHS
    use mod_integer_types
    use mpi_mod

    implicit none
    
    type(type_RHS)        :: sol_vec
    type(type_PRECOND)    :: pc

    integer               :: counts, ierr
    integer(kind=int_all) :: i  
    
    sol_vec%val(1:sol_vec%n) = 0.d0
    if (pc%my_id_n.eq.0) then
      do i = 1, pc%rhs%n
        sol_vec%val(pc%row_index(i)) = pc%rhs%val(i)*pc%row_factor
      enddo
    endif
    counts = sol_vec%n
    call MPI_AllReduce(MPI_IN_PLACE,sol_vec%val,counts,MPI_DOUBLE_PRECISION,MPI_SUM,pc%comm,ierr)
    
    return
  end subroutine gather_solution  

!> Deallocate arrays and reset to the default values
  subroutine reset_preconditioner(pc)
    use data_structure, only: type_PRECOND
    implicit none

    type(type_PRECOND) :: pc !, pc_def

    if (pc%initialized) then

      if (pc%structured) then

        deallocate(pc%rhs%val)
        pc%rhs%val => Null()
        deallocate(pc%row_index)
        pc%row_index => Null()
        deallocate(pc%send_counts, pc%recv_counts)
        pc%send_counts => Null()
        pc%recv_counts => Null()
        deallocate(pc%send_disp, pc%recv_disp)
        pc%send_disp => Null()
        pc%recv_disp => Null()
        deallocate(pc%istart, pc%ifinish)
        pc%istart => Null()
        pc%ifinish => Null()
        deallocate(pc%n_per_rank)
        pc%n_per_rank => Null()

        deallocate(pc%mat%val)
        deallocate(pc%mat%irn)
        deallocate(pc%mat%jcn)
        pc%mat%val => Null()
        pc%mat%irn => Null()
        pc%mat%jcn => Null()

        if (pc%mat%scaled) then
          deallocate(pc%mat%column_scaling)
          pc%mat%column_scaling => Null()
        endif

        pc%mat%scaled = .false.
        pc%mat%row_distributed = .false.
        pc%mat%col_distributed = .false.
        pc%mat%indexing = 1
        pc%mat%block_size = 1
        pc%structured = .false.

      endif

      deallocate(pc%mode_families_ranks)
      pc%mode_families_ranks => Null()
      deallocate(pc%mode_families_modes)
      pc%mode_families_modes => Null()
      deallocate(pc%mode_set)
      pc%mode_set => Null()

      pc%initialized = .false.

    endif

    return

  end subroutine reset_preconditioner


end module mod_preconditioner
