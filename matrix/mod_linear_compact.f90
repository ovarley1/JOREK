!> Linear run by mode family, compaction of the global system (linear_compact): the global matrix, its right-hand side and solution
!! vectors and the preconditioner hold the evolving harmonics only, and the frozen harmonics have no rows and no columns.
!!
!! The harmonic order inside a node block of the system is the order of the active harmonics, so a pair (cos, sin) of an evolving n
!! sits in two neighbouring slots (2j-1, 2j). The element matrices, the nodes (values, deltas) and everything else keep the full
!! order i_tor = 1..n_tor. Index of (node index_node, variable k, harmonic i_tor) in the system:
!!   n_act*n_var*(index_node-1) + n_act*(k-1) + harm_slot(i_tor)    (harm_slot = 0: frozen, not in the system)
!! With compact_on = .false. the defaults below are the identity (n_act = n_tor), which is the full order of the plain code.
module mod_linear_compact

  use mod_parameters, only : n_tor

  implicit none

  integer, private :: n_act_i                              ! implied-do variable of the identity maps below

  logical :: compact_on = .false.                         !< The global system holds the active harmonics only. Set by setup_linear_compact
  integer :: n_act      = n_tor                           !< Number of harmonics per node block of the system (n_tor if not compact_on)
  integer :: harm_active(n_tor) = (/ (n_act_i, n_act_i = 1, n_tor) /) !< Harmonic i_tor of slot s (s = 1..n_act)
  integer :: harm_slot(n_tor)   = (/ (n_act_i, n_act_i = 1, n_tor) /) !< Slot of harmonic i_tor, 0 if the harmonic is frozen

contains

  !> Check that the run supports linear_compact and build the slot map from harmonic_frozen. Runs on every rank (harmonic_frozen and
  !! the inputs are identical on every rank), after setup_linear_family and before global_matrix_structure. Without a linear run
  !! (freeze_harmonics false) every harmonic is active: an identity map that only exercises the compact code paths.
  subroutine setup_linear_compact(my_id)

    use phys_module, only : linear_compact, linear_pc_evolving, freeze_harmonics, harmonic_frozen
    use phys_module, only : freeboundary, keep_n0_const, gmres, use_newton, eliminate_boundary_dofs
    use mod_parameters, only : jorek_model
    use mpi_mod

    implicit none

    integer, intent(in) :: my_id

    integer             :: i, k, ierr
    character(len=120)  :: why

    if (.not. linear_compact) return

    ! --- Everything the compact layout does not handle stops the run here, on every rank
    why = ''
    if (jorek_model /= 183)                                   why = 'only model 183 is supported'
    if ((len_trim(why) == 0) .and. (.not. gmres))             why = 'it needs the iterative solver (gmres=.t.)'
    if ((len_trim(why) == 0) .and. freeboundary)              why = 'it does not support free boundary'
    if ((len_trim(why) == 0) .and. keep_n0_const)             why = 'it does not support keep_n0_const'
    if ((len_trim(why) == 0) .and. use_newton)                why = 'it does not support the Newton solver'
    if ((len_trim(why) == 0) .and. eliminate_boundary_dofs)   why = 'it does not support eliminate_boundary_dofs'
#ifdef USE_PETSC
    if (len_trim(why) == 0)                                   why = 'it does not support USE_PETSC (direct PETSc assembly)'
#endif
#ifdef JECCD
    if (len_trim(why) == 0)                                   why = 'it does not support JECCD'
#endif
#ifdef DIRECT_CONSTRUCTION
    if (len_trim(why) == 0)                                   why = 'it does not support USE_DIRECT_CONSTRUCTION=1'
#endif
#ifdef USE_BICGSTAB
    if (len_trim(why) == 0)                                   why = 'it does not support USE_BICGSTAB=1'
#endif
    if (freeze_harmonics) then
      if ((len_trim(why) == 0) .and. (.not. linear_pc_evolving)) why = 'a linear run needs linear_pc_evolving=.t. (families in slots)'
      if ((len_trim(why) == 0) .and. (.not. harmonic_frozen(1))) why = 'the harmonic n=0 (i_tor=1) must be frozen'
      do k = 1, (n_tor - 1)/2
        if ((len_trim(why) == 0) .and. (harmonic_frozen(2*k) .neqv. harmonic_frozen(2*k + 1))) then
          why = 'a cos/sin pair of harmonics is split by the freeze'
        endif
      enddo
    endif

    if (len_trim(why) /= 0) then
      if (my_id == 0) write(*,*) 'FATAL: linear_compact: ', trim(why)
      call MPI_Abort(MPI_COMM_WORLD, 6, ierr)
      stop
    endif

    ! --- Slot map: the active harmonics keep their order
    n_act = 0
    harm_active(1:n_tor) = 0
    do i = 1, n_tor
      if (freeze_harmonics .and. harmonic_frozen(i)) then
        harm_slot(i) = 0
      else
        n_act = n_act + 1
        harm_active(n_act) = i
        harm_slot(i)       = n_act
      endif
    enddo
    compact_on = .true.

    if (my_id == 0) then
      write(*,'(A,I0,A,I0,A)') ' linear_compact: the global system holds ', n_act, ' of ', n_tor, ' harmonics per node block'
      if (.not. freeze_harmonics) write(*,*) 'linear_compact: no linear run: every harmonic is active (identity map, smoke test only)'
    endif

  end subroutine setup_linear_compact

  !> Zero the increments (deltas) of the frozen harmonics on every node. The element matrix sums them over all harmonics (the time
  !! derivative terms), and a restart may hold stale values. In the plain code update_values forces them to zero at every step
  !! (through the zeroed right-hand side), with linear_compact the frozen harmonics are not in the solution vector and this does it.
  subroutine zero_frozen_deltas(node_list)

    use data_structure, only : type_node_list
    use phys_module,    only : harmonic_frozen, freeze_harmonics

    implicit none

    type(type_node_list), intent(inout) :: node_list

    integer :: i, in

    if (.not. freeze_harmonics) return

    do i = 1, node_list%n_nodes
      do in = 1, n_tor
        if (harmonic_frozen(in)) node_list%node(i)%deltas(in,:,:) = 0.d0
      enddo
    enddo

  end subroutine zero_frozen_deltas

end module mod_linear_compact
