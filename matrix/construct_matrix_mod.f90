module construct_matrix_mod

use mod_parameters, only : n_var, n_order, n_degrees_1d
#ifdef USE_PETSC
#include "petsc/finclude/petsc.h"
use petsc
#endif
implicit none

public :: construct_matrix

contains

!> Construct the main matrix from the contributions of the Bezier elements.
!!
!! The element contributions are determined by element_matrix(_fft). Additional
!! contributions from boundary conditions and the free boundary extension are
!! added by external routine calls.
subroutine construct_matrix(mhd_sim, local_elms, n_local_elms, a_mat, rhs_vec, harmonic_matrix)
  use tr_module
  use phys_module
  use data_structure, only: type_SP_MATRIX, type_RHS, type_element, type_node, thread_struct
  use mod_simulation_data, only: type_MHD_SIM
  use vacuum, only: sr
  use nodes_elements
  use mpi_mod
  use omp_lib
  use mod_axis_treatment, only: transform_basis_for_axis_element, penalize_dof_on_axis
  use mod_boundary_conditions, only: boundary_conditions
  use mod_fix_axis_nodes, only: fix_nodes_on_axis
  use mod_linear_compact, only: compact_on
  use vacuum_response, only: vacuum_boundary_integral
  use global_distributed_matrix, only: global_matrix_structure_vacuum
#ifdef USE_PETSC
  use mod_petsc, only: petsc_create_matrix
#endif
  implicit none

#include "r3_info.h"

  ! --- Routine parameters
  logical,              intent(in)    :: harmonic_matrix
  type(type_SP_MATRIX), intent(inout) :: a_mat !< Sparse matrix to be constructed
  type(type_RHS),       intent(inout) :: rhs_vec !< Right-hand side vector to be constructed
  type(type_MHD_SIM),   intent(in)    :: mhd_sim !< MHD simulation data structure containing a lot of information about the simulation
  integer, dimension(:), pointer      :: local_elms !< List of local element indices in the simulation
  integer                             :: n_local_elms ! Number of local elements in the simulation

  ! --- Internal variables
  type (type_element)               :: element
  type (type_node)                  :: nodes(n_vertex_max), aux_nodes(n_vertex_max)
  type (type_element)               :: element_father
  type (type_node)                  :: nodes_father(n_vertex_max)
  real*8,              allocatable  :: rhs_local(:) !< Local right-hand side vector for the current process

  integer                           :: ife, ielm
  integer                           :: i_v(n_var)
  integer, allocatable              :: i_harm(:)
  integer                           :: mpi_ierr
  integer                           :: omp_nthreads, omp_tid, n_tor_local
  integer                           :: my_ind_min, my_ind_max
  integer                           :: node_out(n_vertex_max)
  integer                           :: my_id
#ifdef USE_PETSC
  PetscErrorCode                    :: ierr
#endif
  integer                           :: xcase2
  real*8                            :: R_axis
  real*8                            :: Z_axis
  real*8                            :: psi_axis
  real*8                            :: psi_bnd
  real*8                            :: R_xpoint(2)
  real*8                            :: Z_xpoint(2)
  real*8                            :: psi_xpoint(2)
  logical                           :: xpoint2
  integer                           :: ierr_c
  real*8                            :: tw0, tw1                              ! wall-clock timers of one thread (s)
  real*8                            :: tw_start, tw_alloc, tw_loop, tw_bc    ! wall-clock timers of the call (s)
  logical                           :: lr_timers                             ! print the timers of the main matrix of a linear run
  logical, save                     :: skip_warned = .false.

  ! --- Timing call
  call r3_info_begin (r3_info_index_0, 'construct_matrix')
  tw_start = 0.d0
  !$ tw_start = omp_get_wtime()

  ! --- linear_compact is set up by setup_linear_compact (jorek2_main): any driver that does not call it must not run with the key on
  if (linear_compact .neqv. compact_on) then
    if (mhd_sim%my_id .eq. 0) write(*,*) 'FATAL: linear_compact is on but its slot map is not set up by this driver'
    call MPI_Abort(MPI_COMM_WORLD, 6, ierr_c)
    stop
  endif

  ! --- Same on every rank (the timers are reduced collectively): main matrix of a linear run by mode family
  lr_timers = freeze_harmonics .and. (.not. harmonic_matrix)

  ! --- Linear run by mode family, GMRES, PC from the evolving harmonics: the rows and columns of the frozen harmonics never reach
  ! --- the iterates (gather_solution and prec in mod_gmres2 keep every vector exactly 0 on them), so the FFT element matrix may
  ! --- leave them out. Set at every call, identical on every rank (all inputs are broadcast).
  skip_frozen_pairs = linear_skip_frozen .and. freeze_harmonics .and. linear_pc_evolving .and. gmres .and. (.not. use_newton) &
                      .and. (.not. harmonic_matrix)
  if (linear_skip_frozen .and. (.not. skip_frozen_pairs) .and. (.not. harmonic_matrix) .and. (.not. skip_warned)) then
    skip_warned = .true.
    if (mhd_sim%my_id .eq. 0) write(*,*) 'WARNING: linear_skip_frozen is ignored: it needs linear_run, linear_pc_evolving and gmres, and no Newton solver'
  endif

  ! --- Extract parameters from mhd_sim
  my_id           = mhd_sim%my_id
  xpoint2         = mhd_sim%es%xpoint
  xcase2          = mhd_sim%es%xcase
  R_axis          = mhd_sim%es%R_axis
  Z_axis          = mhd_sim%es%Z_axis
  psi_axis        = mhd_sim%es%psi_axis
  psi_bnd         = mhd_sim%es%psi_bnd
  R_xpoint(1:2)   = mhd_sim%es%R_xpoint(1:2)
  Z_xpoint(1:2)   = mhd_sim%es%Z_xpoint(1:2)
  psi_xpoint(1:2) = mhd_sim%es%psi_xpoint(1:2)

  ! --- Printout
  if (mhd_sim%my_id .eq. 0) then
    if(.NOT. harmonic_matrix) then
      write(*,*) '****************************************'
      write(*,*) '*       construct global matrix        *'
      write(*,*) '*       using CPU                      *'
      write(*,*) '****************************************'
    else
      write(*,*) '****************************************'
      write(*,*) '*        construct PC matrix           *'
      write(*,*) '****************************************'
    endif
  endif

  ! --- Memory allocation of temporary buffers
  call new_thread_buffers() 

  ! --- Allocation/Reallocation of the sparse matrix and right-hand side vector
#ifdef USE_PETSC
  if (.not. harmonic_matrix) then !TODO: Might be unnessecary ... constant protection by harmonic_matrix...
    ! Direct PETSc assembly: skip irn/jcn/val, use PETSc MPIBAIJ matrix
    if (.not. a_mat%petsc_assembled) then
      call petsc_create_matrix(a_mat%petsc_A, a_mat)
      a_mat%petsc_assembled = .true.
    else
      PetscCallA(MatZeroEntries(a_mat%petsc_A, ierr))
    endif
  else
#endif
  if (associated(a_mat%irn)) call tr_deallocatep(a_mat%irn, "irn", CAT_DMATRIX)
  if (associated(a_mat%jcn)) call tr_deallocatep(a_mat%jcn, "jcn", CAT_DMATRIX)
  if (associated(a_mat%val)) call tr_deallocatep(a_mat%val, "val", CAT_DMATRIX)

  call tr_allocatep(a_mat%irn, Int1, a_mat%nnz, "irn", CAT_DMATRIX)
  call tr_allocatep(a_mat%jcn, Int1, a_mat%nnz, "jcn", CAT_DMATRIX)
  call tr_allocatep(a_mat%val, Int1, a_mat%nnz, "val", CAT_DMATRIX)

  a_mat%irn(1:a_mat%nnz) = 0
  a_mat%jcn(1:a_mat%nnz) = 0
  a_mat%val(1:a_mat%nnz) = 0.0d0
#ifdef USE_PETSC
  endif
#endif

  if (associated(rhs_vec%val)) call tr_deallocatep(rhs_vec%val,"rhs",CAT_DMATRIX)
  call tr_allocatep(rhs_vec%val, Int1, a_mat%ng, "rhs", CAT_DMATRIX)
  rhs_vec%val(:) = 0.0d0 

  call tr_allocate(rhs_local, Int1, a_mat%ng, "rhs_local", CAT_DMATRIX)
  rhs_local  = 0.d0

  ! --- Element buffers (thread_struct) and matrix arrays are allocated here: the expected peak of the construction
  if (lr_timers) call tr_print_memsize("CMAfterAlloc")
  tw_alloc = 0.d0
  !$ tw_alloc = omp_get_wtime() - tw_start

  if (mhd_sim%freeboundary .and. (mhd_sim%sr_n_tor /= 0 ) ) then
    call global_matrix_structure_vacuum(mhd_sim%node_list, mhd_sim%bnd_node_list, a_mat, i_tor_min=1, i_tor_max=n_tor)
  endif

  ! --- Matrix construction
  my_ind_min = a_mat%index_min(my_id+1)
  my_ind_max = a_mat%index_max(my_id+1)

  n_tor_local = a_mat%i_tor_max - a_mat%i_tor_min + 1

  ! --- Start the OMP parallel region
    !$omp parallel default(none) &
    !$omp   shared(n_local_elms, local_elms, n_tor_local, element_list, node_list, aux_node_list, my_id, xpoint2, xcase2, &
    !$omp           R_axis, Z_axis, psi_axis, psi_bnd, R_xpoint, Z_xpoint, thread_struct, harmonic_matrix, &
    !$omp            a_mat, rhs_local, rhs_vec, &
    !$omp            fix_axis_nodes, treat_axis, my_ind_min, my_ind_max, bc_natural_open, bc_natural_flux, n_tor_fft_thresh, refinement) &
    !$omp   private(ife, ielm, element, element_father, nodes, nodes_father, aux_nodes, node_out, omp_nthreads, omp_tid, i_v, i_harm, tw0, tw1) 
    
    ! --- Fetch the OMP thread information
#ifdef _OPENMP
    omp_nthreads = omp_get_num_threads()
    omp_tid      = 1 + omp_get_thread_num()
#else
    omp_nthreads = 1
    omp_tid      = 1
#endif
    
    ! --- Initialise i_v and i_harm in case of treat_axis
    if (treat_axis) call initialize_v_and_harm(i_v, i_harm, a_mat%i_tor_min, a_mat%i_tor_max, n_tor_local)

    ! --- Element loop
    !$omp do schedule(runtime)
    do ife = 1, n_local_elms
      ielm = local_elms(ife)
      !$ tw0 = omp_get_wtime()

      ! --- Get element and nodes
      call get_element_and_nodes(ielm, element, element_father, nodes, nodes_father, aux_nodes, harmonic_matrix)

      ! --- Build element matrix
      !$ tw1 = omp_get_wtime()
      call build_element_matrix(element, nodes, xpoint2, xcase2, R_axis, Z_axis, psi_axis, psi_bnd, R_xpoint, Z_xpoint, omp_tid, ife, ielm, n_local_elms, node_list, a_mat%i_tor_min, a_mat%i_tor_max, aux_nodes, harmonic_matrix)
      !$ call tm_add(omp_tid, 10, omp_get_wtime() - tw1)

      ! --- Transform basis functions for axis nodes if required
      if(treat_axis .and. (nodes(1)%axis_node .or. nodes(2)%axis_node .or. nodes(3)%axis_node .or. nodes(4)%axis_node) ) then
        !$ tw1 = omp_get_wtime()
        call transform_basis_for_axis_element(nodes, ielm, thread_struct(omp_tid)%ELM, thread_struct(omp_tid)%RHS, i_v, n_var, i_harm, n_tor_local)
        !$ call tm_add(omp_tid, 11, omp_get_wtime() - tw1)
      endif

      ! TODO: Implement the diagnostic print_element_rhs subroutine

      ! --- Define element nodes
      call define_element_nodes(ielm, node_out, element, nodes, element_father, nodes_father, omp_tid, harmonic_matrix)

      ! --- Add element matrix to the global matrix
      !$ tw1 = omp_get_wtime()
      call add_to_a_mat(element, node_out, a_mat, rhs_local, my_ind_min, my_ind_max, omp_tid)
      !$ call tm_add(omp_tid, 12, omp_get_wtime() - tw1)

      !$ call tm_add(omp_tid, 16, omp_get_wtime() - tw0)
      thread_struct(omp_tid)%tm_n = thread_struct(omp_tid)%tm_n + 1

    enddo !< Element loop
    !$omp end do
    !$omp end parallel
    tw_loop = 0.d0
    !$ tw_loop = omp_get_wtime() - tw_start - tw_alloc

    ! --- Matrix and RHS are assembled in the element loop: memory after the scatter
    if (lr_timers) call tr_print_memsize("CMAfterScat")

#ifdef USE_PETSC
    if (.not. harmonic_matrix) then
      ! Flush element contributions (ADD_VALUES) before BCs (INSERT_VALUES)
      PetscCallA(MatAssemblyBegin(a_mat%petsc_A, MAT_FLUSH_ASSEMBLY, ierr))
      PetscCallA(MatAssemblyEnd(a_mat%petsc_A, MAT_FLUSH_ASSEMBLY, ierr))
    else
#endif
    ! --- Memory tracking
    call tr_vnorms("cm_A_bef_bc", a_mat%val, a_mat%nnz)
#ifdef USE_PETSC
    endif
#endif

    ! --- Apply boundary conditions.
    if (.not. eliminate_boundary_dofs) then !TODO: This is trivially not true in all cases and needs to be treated more precisely
      call boundary_conditions(my_id, node_list, element_list,  bnd_node_list,local_elms, n_local_elms,  &
                              my_ind_min, my_ind_max, rhs_local, xpoint2, xcase2, R_axis, Z_axis,        & 
                              psi_axis, psi_bnd, R_xpoint, Z_xpoint, psi_xpoint, a_mat)
    endif

    if (fix_axis_nodes) then
      call fix_nodes_on_axis(node_list, element_list, local_elms, n_local_elms, my_ind_min, my_ind_max, a_mat)
    elseif(treat_axis)then
      call penalize_dof_on_axis(node_list, 4, element_list, local_elms, n_local_elms, my_ind_min, my_ind_max, a_mat)
    endif

    tw_bc = 0.d0
    !$ tw_bc = omp_get_wtime() - tw_start - tw_alloc - tw_loop

#ifdef USE_PETSC
    if (.not. harmonic_matrix) then
      ! Flush BC contributions (INSERT_VALUES) before vacuum (ADD_VALUES)
      PetscCallA(MatAssemblyBegin(a_mat%petsc_A, MAT_FLUSH_ASSEMBLY, ierr))
      PetscCallA(MatAssemblyEnd(a_mat%petsc_A, MAT_FLUSH_ASSEMBLY, ierr))
    else
#endif
    ! --- Memory tracking
    call tr_vnorms("cm_A_aft_bc", a_mat%val, a_mat%nnz)
#ifdef USE_PETSC
    endif
#endif

    ! --- Add vacuum response (boundary integral) for free boundary computations
    if ( freeboundary .and. ( sr%n_tor /= 0 ) ) then
      call vacuum_boundary_integral(my_id, bnd_node_list, node_list, bnd_elm_list, freeboundary_equil, &
                                    resistive_wall, my_ind_min, my_ind_max, rhs_local, tstep, index_now, a_mat)
    endif

#ifdef USE_PETSC
    if (.not. harmonic_matrix) then
      ! Final assembly of PETSc matrix
      PetscCallA(MatAssemblyBegin(a_mat%petsc_A, MAT_FINAL_ASSEMBLY, ierr))
      PetscCallA(MatAssemblyEnd(a_mat%petsc_A, MAT_FINAL_ASSEMBLY, ierr))

      !call MatSetOption(a_mat%petsc_A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE, ierr)
    endif
#endif

    if ( .not. harmonic_matrix ) then
#ifdef COMPARE_ELEMENT_MATRIX
      ! TODO: Create a subroutine to compare the element matrix with the right-hand side vector
      !call summarise_element_matrix_comparison(a_mat, rhs_vec, my_id, index_now)
#endif
#ifdef NORMTRACE
      ! TODO: Reintroduce normtrace subroutine
      !call normtrace(...)
#endif
    endif
    
  ! --- Collect the right-hand side vector from all processes 
  call MPI_AllReduce(RHS_local,rhs_vec%val,a_mat%ng,MPI_DOUBLE_PRECISION,MPI_SUM,a_mat%comm,mpi_ierr)
  rhs_vec%n  = a_mat%ng

  ! --- Debug: sums over matrix and RHS, to compare two builds of the same case (all ranks take part in the reduction)
  if (linear_check_sums .and. (.not. harmonic_matrix)) call report_linear_check_sums(a_mat%comm, a_mat, rhs_vec)

  ! --- Timers of the element loop (all ranks take part in the reduction) and of the rest of the call
  if (lr_timers) call report_construct_timers(a_mat%comm, n_local_elms, tw_start, tw_alloc, tw_loop, tw_bc)

  ! --- Check if the matrix is distributed correctly
#ifdef USE_PETSC
  if (harmonic_matrix) call check_if_distributed(a_mat)
#else
  call check_if_distributed(a_mat)
#endif

  ! --- Handle all nessecarry cleanup operations
  call tr_locvnorms("cm_BCRhs",rhs_vec%val,a_mat%ng)
  call tr_debug_write("ndof",a_mat%ng)

  ! --- Timing
  call r3_info_end(r3_info_index_0)
  call tr_print_memsize("EndConstM")
end subroutine construct_matrix


!> Helps to interprete an element matrix index
subroutine  decrypt_index(ind, ivertex, iorder, ivar, itor)

  use mod_parameters,  only : n_tor, jorek_model, n_vertex_max, n_var, n_degrees

  integer, intent(in)  :: ind     !< Element matrix index
  integer, intent(out) :: ivertex !< Vertex index
  integer, intent(out) :: iorder  !< Degree of freedom
  integer, intent(out) :: ivar    !< Variable index
  integer, intent(out) :: itor    !< Toroidal mode index

  integer :: ind2

  ind2 = ind

  ivertex = ( ind2 - 1 ) / ( n_tor*n_var*n_degrees ) + 1
  ind2 = ind2 - ( ivertex - 1 ) * ( n_tor*n_var*n_degrees )

  iorder = ( ind2 - 1 ) / ( n_tor*n_var ) + 1
  ind2 = ind2 - ( iorder - 1 ) * ( n_tor*n_var )

  ivar = ( ind2 - 1 ) / ( n_tor ) + 1
  ind2 = ind2 - ( ivar - 1 ) * ( n_tor )

  itor = ind2

end subroutine decrypt_index

subroutine build_element_matrix(element, nodes, xpoint2, xcase2, R_axis,         &
      &                             Z_axis, psi_axis, psi_bnd, R_xpoint, Z_xpoint,   &
      &                             omp_tid, ife, ielm, n_local_elms, node_list, i_tor_min, i_tor_max, &
                                    aux_nodes, harmonic_matrix)

  use mod_parameters,           only : n_tor, jorek_model, n_vertex_max, n_degrees, unified_element_matrix, n_var
  use phys_module,              only : bc_natural_open, bc_natural_flux, n_tor_fft_thresh, grid_to_wall, n_wall_blocks, keep_n0_const
  use phys_module,              only : freeze_harmonics, harmonic_frozen
  use mod_linear_compact,       only : compact_on
  USE data_structure,           only : type_element, type_node, type_node_list, thread_struct
  use mod_boundary_matrix_open, only : boundary_matrix_open
  use mod_elt_matrix,           only : element_matrix
  use mod_elt_matrix_fft,   only : element_matrix_fft
  use mpi_mod
  implicit none

  type (type_element),              intent(inout)  :: element
  type (type_node),                 intent(inout)  :: nodes(n_vertex_max)
  logical,                          intent(in)     :: xpoint2
  integer,                          intent(in)     :: xcase2
  real*8,                           intent(in)     :: R_axis
  real*8,                           intent(in)     :: Z_axis
  real*8,                           intent(in)     :: psi_axis
  real*8,                           intent(in)     :: psi_bnd
  real*8,                           intent(in)     :: R_xpoint(2)
  real*8,                           intent(in)     :: Z_xpoint(2)
  integer,                          intent(in)     :: omp_tid
  integer,                          intent(in)     :: ife
  integer,                          intent(in)     :: ielm
  integer,                          intent(in)     :: n_local_elms
  integer,                          intent(in)     :: i_tor_min   
  integer,                          intent(in)     :: i_tor_max   
  TYPE (type_node_list),            intent(in)     :: node_list
  type (type_node), optional,       intent(inout)  :: aux_nodes(n_vertex_max)
  logical,                          intent(in)    :: harmonic_matrix
  
  ! -- internal parameters
  integer :: iv, iv2, iv3, iv4, inode1, inode2, inode3, inode4, i, j
  integer :: vertex(2), direction(n_degrees_1d), bnd1, bnd2, side1, side2
  type (type_node) :: bnd_nodes(n_vertex_max)   ! the edge's nodes, rotated for boundary_matrix_open
  integer :: i_max   
  integer :: n_tor_local

  ! --- Call element_matrix
  ! --- Call element_matrix
  if ( ( (i_tor_min .eq. 1) .and. (i_tor_max .eq. n_tor) .and. (n_tor .ge. n_tor_fft_thresh) )   &
    .or. (unified_element_matrix) ) then
    ! (use the FFT element matrix construction or the unified one in case it has been combined
    ! for the respective model)
    call element_matrix_fft(element,nodes, xpoint2, xcase2, R_axis, Z_axis, psi_axis, psi_bnd,   &
      R_xpoint, Z_xpoint, thread_struct(omp_tid)%ELM, thread_struct(omp_tid)%RHS, omp_tid,       &
      thread_struct(omp_tid)%ELM_p, thread_struct(omp_tid)%ELM_n, thread_struct(omp_tid)%ELM_k,  &
      thread_struct(omp_tid)%ELM_kn, thread_struct(omp_tid)%RHS_p, thread_struct(omp_tid)%RHS_k, &
      thread_struct(omp_tid)%eq_g, thread_struct(omp_tid)%eq_s, thread_struct(omp_tid)%eq_t,     &
      thread_struct(omp_tid)%eq_p, thread_struct(omp_tid)%eq_ss, thread_struct(omp_tid)%eq_st,   &
      thread_struct(omp_tid)%eq_tt, thread_struct(omp_tid)%delta_g,                              &
      thread_struct(omp_tid)%delta_s, thread_struct(omp_tid)%delta_t, i_tor_min, i_tor_max,      &
      aux_nodes, thread_struct(omp_tid)%ELM_pnn)
  else
    ! (use the element matrix by toroidal integration in case of very few harmonics or in case
    ! of direct construction of the harmonic matrices used in preconditioning)
    call element_matrix(element,nodes, xpoint2, xcase2, R_axis, Z_axis, psi_axis, psi_bnd,       &
      R_xpoint, Z_xpoint, thread_struct(omp_tid)%ELM, thread_struct(omp_tid)%RHS, omp_tid,       &
      i_tor_min, i_tor_max, aux_nodes)
  endif

  ! --- Apply sheath boundary conditions at the targets
  if (bc_natural_open) then
    ! --- Loop over the 4 nodes
    do iv = 1, n_vertex_max
      
      iv2  = mod(iv, n_vertex_max) + 1
      iv3 = mod(iv2, n_vertex_max) + 1
      iv4 = mod(iv3, n_vertex_max) + 1
      
      inode1 = element%vertex(iv)
      inode2 = element%vertex(iv2)
      inode3 = element%vertex(iv3)
      inode4 = element%vertex(iv4)
      
      bnd1 = node_list%node(inode1)%boundary
      bnd2 = node_list%node(inode2)%boundary
      
      ! --- carry on only if on boundary
      if ( (bnd1 .eq. 0) .or. (bnd2 .eq. 0)) cycle
      
      bnd_nodes(1) = node_list%node(inode1)
      bnd_nodes(2) = node_list%node(inode2)
      bnd_nodes(3) = node_list%node(inode3)
      bnd_nodes(4) = node_list%node(inode4)

      vertex    = (/ iv, iv2 /)
      
      if ( (grid_to_wall) .and. (n_wall_blocks .gt. 0) ) then
        side1 = 0                  ; side2 = 0
        if (bnd1 .eq. 1) side1 = 2 ; if (bnd2 .eq. 1) side2 = 2
        if (bnd1 .eq.11) side1 = 2 ; if (bnd2 .eq.11) side2 = 2
        if (bnd1 .eq. 5) side1 = 3 ; if (bnd2 .eq. 5) side2 = 3
        if (bnd1 .eq.15) side1 = 3 ; if (bnd2 .eq.15) side2 = 3
        if (bnd1 .eq. 2) side1 = 3 ; if (bnd2 .eq. 2) side2 = 3
        if (bnd1 .eq.12) side1 = 2 ; if (bnd2 .eq.12) side2 = 2
        if (bnd1 .eq. 4) side1 = 2 ; if (bnd2 .eq. 4) side2 = 2
        direction(1) = 1
        if     ( (side1 .eq. 2) .or. (side2 .eq. 2) ) then
          direction(2) = 2
          if (n_order .ge. 5) direction(3) = 5
        elseif ( (side1 .eq. 3) .or. (side2 .eq. 3) ) then
          direction(2) = 3
          if (n_order .ge. 5) direction(3) = 6
        endif
        ! --- This should never happen, but just in case...
        if (     ((side1 .eq. 2) .and. (side2 .eq. 3)) &
            .or. ((side1 .eq. 3) .and. (side2 .eq. 2)) ) then
          write(*,'(A,4i8)') 'WARNING: boundary_matrix_open, boundary element incoherent ',&
                             inode1,node_list%node(inode1)%boundary,inode2,node_list%node(inode2)%boundary  
          cycle
        endif
      else
        ! --- The target has boundary 1 or 3
        direction(1) = 1
        if ((bnd1 .eq. 9) .and. (bnd2 .eq. 9)) then ! for removed triangle, we get 9 - 9 boundary, where it could be either s=const or t=const
          ! so we use the fact that the nodes are consistently defining which sides are the s and t sides (identical to that in grids/mod_boundary.f90)
          if ( mod(iv,2) == 1 ) then
            direction(2) = 2
          else
            direction(2) = 3
          end if
            
        elseif ( (  ((bnd1 .eq. 1) .or. (bnd1 .eq. 3)) .and. ((bnd2 .eq. 1) .or. (bnd2 .eq. 3))  ) &
            .or. (  ((bnd1 .eq. 1) .or. (bnd1 .eq. 9)) .and. ((bnd2 .eq. 1) .or. (bnd2 .eq. 9))  ) &
            .or. (  ((bnd1 .eq. 4) .or. (bnd1 .eq. 9)) .and. ((bnd2 .eq. 4) .or. (bnd2 .eq. 9))  ) &
            .or. (  ((bnd1 .eq. 1) .or. (bnd1 .eq. 4)) .and. ((bnd2 .eq. 4) .or. (bnd2 .eq. 1))  ) ) then
          
          direction(2) = 2
          if (n_order .ge. 5) direction(3) = 5
          
        elseif (  ((bnd1 .eq. 5) .or. (bnd1 .eq. 9)) .and. ((bnd2 .eq. 5) .or. (bnd2 .eq. 9)) ) then
          
          direction(2) = 3
          if (n_order .ge. 5) direction(3) = 6
          
        elseif (  ((bnd1 .eq. 2) .or. (bnd1 .eq. 3)) .and. ((bnd2 .eq. 2) .or. (bnd2 .eq. 3)) ) then
          
          direction(2) = 3
          if (n_order .ge. 5) direction(3) = 6
          
        else
          write(*,'(A,4i8)') 'WARNING: boundary_matrix_open, boundary element not included ',&
                             inode1,node_list%node(inode1)%boundary,inode2,node_list%node(inode2)%boundary  
          cycle
        endif
      endif

      ! sanity check of the above code: the side between node 1 and 2 (iv=1), and nodes 3 and 4 (iv=3) should have direction(2)=2 (i.e. t=const.), while 
      ! the other two sides between nodes 2 and 3 (iv=2), and 4 and 1 (iv=4) should have direction(2)=3 (i.e. s=cont.)
      if (((direction(2)==2) .and. ( mod(iv,2) /= 1 )) .or. ((direction(2)==3) .and. (mod(iv,2) /= 0))) then
        !$omp critical
        write(*,"(A,3I6)") "ERROR: There seems to be an inconsistency in direction(2) in matrix/construct_matrix_mod.f90, (direction(2) / iv / node)=", direction(2),iv,inode1
        !$omp end critical
      end if
        

      ! --- Build matrix elements for boundary
#if JOREK_MODEL == 183
      call boundary_matrix_open(vertex, direction, element, bnd_nodes, & 
                                xpoint2, xcase2, R_axis, Z_axis, psi_axis, psi_bnd, R_xpoint, Z_xpoint, &
                                thread_struct(omp_tid)%ELM, thread_struct(omp_tid)%RHS, i_tor_min, i_tor_max, ielm)
#else
      call boundary_matrix_open(vertex, direction, element, bnd_nodes, &
                                xpoint2, xcase2, R_axis, Z_axis, psi_axis, psi_bnd, R_xpoint, Z_xpoint, &
                                thread_struct(omp_tid)%ELM, thread_struct(omp_tid)%RHS, i_tor_min, i_tor_max)
#endif

    enddo
  endif

  n_tor_local = i_tor_max - i_tor_min + 1
  if (keep_n0_const) then
    i_max = n_degrees * n_vertex_max * n_var * n_tor_local
#ifdef JECCD
    i_max = n_degrees * n_vertex_max * (n_var - 1) * n_tor_local
#endif
    do i = 1, i_max, n_tor_local
      thread_struct(omp_tid)%ELM(i, i) = 1.d15
    enddo
  endif

  ! --- Linear run by mode family: the same large diagonal entry on every frozen harmonic. The harmonic is the fastest
  ! --- index of ELM, and the block may start at i_tor_min > 1, so the row is mapped to its global harmonic.
  ! --- Not with linear_compact: the frozen harmonics are not in the system, the scatter leaves them out.
  if ( freeze_harmonics .and. (.not. compact_on) ) then
    i_max = n_degrees * n_vertex_max * n_var * n_tor_local
    do i = 1, i_max
      if ( harmonic_frozen(i_tor_min + mod(i-1, n_tor_local)) ) thread_struct(omp_tid)%ELM(i, i) = 1.d15
    enddo
  endif

  !--- Diagnostic comparison of element matrix constructed with and without FFT
  ! TODO: Reimplement compare_elm_rhs diagnostic

end subroutine build_element_matrix


subroutine get_element_and_nodes(ielm, element, element_father, nodes, nodes_father, aux_nodes, harmonic_matrix)
  use phys_module, only: refinement
  use mod_parameters, only: n_vertex_max
  use data_structure
  use nodes_elements
  implicit none

  integer, intent(in) :: ielm
  type(type_element), intent(out) :: element
  type(type_element), intent(out) :: element_father
  type(type_node)   , intent(out) :: nodes(n_vertex_max)
  type(type_node)   , intent(out) :: nodes_father(n_vertex_max)
  type(type_node)   , intent(out) :: aux_nodes(n_vertex_max)
  logical, intent(in) :: harmonic_matrix

  integer :: ifather, inode, iv, inode_father

  element = element_list%element(ielm)
  
  ! --- Define nodes (mhd_sim% depends on whether our element has been refined)
  if (refinement .and. .not. harmonic_matrix) then

    ifather = element_list%element(ielm)%father

    if (ifather .ne. 0) then
      element_father = element_list%element(ifather)
      do iv = 1, n_vertex_max
        inode_father=element_father%vertex(iv)
        nodes_father(iv) = node_list%node(inode_father)
      enddo
    endif

  else
  
    do iv = 1, n_vertex_max
        inode   = element%vertex(iv)
        nodes(iv) = node_list%node(inode)
        aux_nodes(iv) = aux_node_list%node(inode)
    enddo

  endif
end subroutine get_element_and_nodes


subroutine add_to_a_mat(element, node_out, a_mat, rhs_local, my_ind_min, my_ind_max, omp_tid)
  use phys_module, only: refinement, eliminate_boundary_dofs
  use mod_parameters, only: n_tor, n_var, n_degrees, n_vertex_max
  use data_structure, only: type_element, type_node, type_node_list, type_SP_MATRIX, thread_struct
  use nodes_elements
  use mod_linear_compact, only: compact_on, n_act, harm_active
  implicit none

  ! Arguments
  type(type_element), intent(in) :: element
  integer, intent(in) :: node_out(:)
  type(type_SP_MATRIX), intent(inout) :: a_mat
  real*8              :: rhs_local(:)
  integer, intent(in) :: my_ind_min, my_ind_max
  integer, intent(in) :: omp_tid

  ! Locals
  integer :: i, j
  integer :: i_order
  integer :: inode1
  integer :: index_node1
  integer :: index_large_i
  integer :: index_ij

  logical :: i_bnd
  integer :: i_bnd_type
  !> Per-variable diagonal of the eliminate_boundary_dofs boundary rows: the
  !! variables' diagonals differ by orders of magnitude, and one shared value
  !! would put some fields' boundary rows far off their bulk scale.
  real*8 :: elm_diagonal_average(n_var)
  integer :: nnz_counter(n_var)
  integer :: jvar

  integer :: n_tor_local
  integer :: cvj, csj, cj

  n_tor_local = a_mat%i_tor_max - a_mat%i_tor_min + 1

  if (eliminate_boundary_dofs) then
    elm_diagonal_average = 0.d0
    nnz_counter = 0
    do i = 1, n_vertex_max
      do i_order = 1, n_degrees
          do j = 1, n_var * n_tor_local
            index_ij = n_tor_local * n_var * n_degrees * (i-1) + n_tor_local * n_var * (i_order-1) + j
            ! j runs over (variable, harmonic), harmonic innermost -- the same
            ! layout the BC rows are built with elsewhere.
            jvar = (j - 1) / n_tor_local + 1
            if (abs(thread_struct(omp_tid)%ELM(index_ij,index_ij)) .ne. 0) &
              nnz_counter(jvar) = nnz_counter(jvar) + 1
            elm_diagonal_average(jvar) = elm_diagonal_average(jvar) &
                                       + abs(thread_struct(omp_tid)%ELM(index_ij,index_ij))
          enddo
      enddo
    enddo
    do jvar = 1, n_var
      if (nnz_counter(jvar) > 0) then
        elm_diagonal_average(jvar) = elm_diagonal_average(jvar) / nnz_counter(jvar)
      else
        elm_diagonal_average(jvar) = 1.d0
      endif
      ! NO absolute floor. max(...,1.d0) is what broke the small-scale variables,
      ! and it is meaningless once the value is per-variable -- the point is that
      ! each variable's own scale sets it. Guard only against a degenerate value.
      if (.not. (elm_diagonal_average(jvar) > 0.d0)) elm_diagonal_average(jvar) = 1.d0
      elm_diagonal_average(jvar) = min(elm_diagonal_average(jvar), 1.d12)
    enddo
  endif

  ! --- We only look at non-refined elements
  if ((.not. refinement) .or. (refinement .and. (element%n_sons .eq. 0))) then

    do i=1,n_vertex_max

      i_bnd = .false.

      inode1 = node_out(i)

      ! --- Get the boundary type of the node
      i_bnd_type = node_list%node(inode1)%boundary
      if (i_bnd_type .ne. 0) i_bnd = .true.

      do i_order = 1, n_degrees

        index_node1 = node_list%node(inode1)%index(i_order)

        index_large_i = n_tor_local * n_var * (index_node1 - 1)
        ! linear_compact: the system holds the active harmonics only, in blocks of n_act (ELM and RHS stay in the full order)
        if (compact_on) index_large_i = n_act * n_var * (index_node1 - 1)

        if ((index_node1 .ge. my_ind_min) .and. (index_node1 .le. my_ind_max)) then

          ! --- RHS assembly
          if (compact_on) then
            do cvj = 1, n_var
              do csj = 1, n_act
                cj       = (cvj-1) * n_act + csj                                        ! row in the block of the system
                index_ij = n_tor_local * n_var * n_degrees * (i-1) + n_tor_local * n_var * (i_order-1) &
                         + (cvj-1) * n_tor_local + harm_active(csj)                     ! index in the ELM matrix (full harmonic order)
                !$omp atomic
                rhs_local(index_large_i+cj) = rhs_local(index_large_i+cj) + thread_struct(omp_tid)%RHS(index_ij)
                !$omp end atomic
              enddo
            enddo
          else if (eliminate_boundary_dofs .and. i_bnd .and. &
                (     (i_bnd_type .eq. 1 .and. (i_order .eq. 1 .or. i_order .eq. 2)) &
                .or.  (i_bnd_type .eq. 2 .and. (i_order .eq. 1 .or. i_order .eq. 3)) &
                .or.  (i_bnd_type .eq. 3 .and. (i_order .eq. 1 .or. i_order .eq. 2 .or. i_order .eq. 3))   )) then

            do j = 1, n_var * n_tor_local
              index_ij = n_tor_local * n_var * n_degrees * (i-1) + n_tor_local * n_var * (i_order-1) + j
              !$omp critical
              rhs_local(index_large_i+j) = 0.d0
              !$omp end critical
            enddo
          else
            do j = 1, n_var * n_tor_local
              index_ij = n_tor_local * n_var * n_degrees * (i-1) + n_tor_local * n_var * (i_order-1) + j
              !$omp atomic
              rhs_local(index_large_i+j) = rhs_local(index_large_i+j) + thread_struct(omp_tid)%RHS(index_ij)
              !$omp end atomic
            enddo
          endif

          ! --- Matrix assembly: dispatch to backend-specific routine
#ifdef USE_PETSC
          if (a_mat%petsc_assembled) then
            call add_block_to_petsc(index_node1, i, i_order, i_bnd, i_bnd_type, &
                                    node_out, a_mat, omp_tid, n_tor_local, &
                                    eliminate_boundary_dofs, elm_diagonal_average)
          else
#endif
            call add_block_to_sp_matrix(index_node1, i, i_order, i_bnd, i_bnd_type, &
                                        node_out, a_mat, omp_tid, n_tor_local, &
                                        my_ind_min, my_ind_max, &
                                        eliminate_boundary_dofs, elm_diagonal_average)
#ifdef USE_PETSC
          endif
#endif

        endif ! my_ind_min < index < my_ind_max

      enddo ! n_degrees

    enddo ! n_vertex_max

  end if
end subroutine add_to_a_mat


!> Assemble element blocks into type_SP_MATRIX (irn/jcn/val arrays).
!! Called from add_to_a_mat for the non-PETSc path.
subroutine add_block_to_sp_matrix(index_node1, i, i_order, i_bnd, i_bnd_type, &
                                   node_out, a_mat, omp_tid, n_tor_local, &
                                   my_ind_min, my_ind_max, &
                                   eliminate_boundary_dofs, elm_diagonal_average)
  use mod_parameters, only: n_var, n_degrees, n_vertex_max
  use data_structure, only: type_node, type_node_list, type_SP_MATRIX, thread_struct
  use nodes_elements
  use mod_locate_irn_jcn
  use mod_linear_compact, only: compact_on, n_act, harm_active
  !$ use omp_lib
  implicit none

  integer, intent(in) :: index_node1, i, i_order
  logical, intent(in) :: i_bnd
  integer, intent(in) :: i_bnd_type
  integer, intent(in) :: node_out(:)
  type(type_SP_MATRIX), intent(inout) :: a_mat
  integer, intent(in) :: omp_tid, n_tor_local
  integer, intent(in) :: my_ind_min, my_ind_max
  logical, intent(in) :: eliminate_boundary_dofs
  real*8, intent(in) :: elm_diagonal_average(:)   !< one per variable

  integer :: j, k, l, k_order
  integer :: knode
  integer :: index_node2
  integer :: index_large_i, index_large_k
  integer :: index_ij, index_kl
  integer :: ijA_position, ilarge2
  logical :: interior, k_bnd
  integer :: k_bnd_type
  integer :: cj, cl, cvj, csj, cvl, csl
  real*8  :: tw2, tw3

  index_large_i = n_tor_local * n_var * (index_node1 - 1)
  if (compact_on) index_large_i = n_act * n_var * (index_node1 - 1)

  do k=1,n_vertex_max

    knode = node_out(k)
    k_bnd = .false.

    k_bnd_type = node_list%node(knode)%boundary
    if (k_bnd_type .ne. 0) k_bnd = .true.

    interior = .not. (i_bnd .or. k_bnd)

    do k_order = 1, n_degrees

      index_node2 = node_list%node(knode)%index(k_order)

      index_large_k = n_tor_local * n_var * (index_node2 - 1)
      if (compact_on) index_large_k = n_act * n_var * (index_node2 - 1)

      call locate_irn_jcn(index_node1,index_node2,my_ind_min,my_ind_max,ijA_position,a_mat)

      if (compact_on) then
        ! linear_compact: only the active harmonics enter the block (setup_linear_compact excludes eliminate_boundary_dofs)
        thread_struct(omp_tid)%synch_buff(1:n_var*n_act*n_var*n_act) = 0.d0
        do cvj = 1, n_var
          do csj = 1, n_act
            cj       = (cvj-1) * n_act + csj                                    ! row in the block of the system
            index_ij = n_tor_local * n_var * n_degrees * (i-1) + n_tor_local * n_var * (i_order-1) &
                     + (cvj-1) * n_tor_local + harm_active(csj)                 ! index in the ELM matrix (full harmonic order)
            do cvl = 1, n_var
              do csl = 1, n_act
                cl       = (cvl-1) * n_act + csl                                ! column in the block of the system
                index_kl = n_tor_local * n_var * n_degrees * (k-1) + n_tor_local * n_var * (k_order-1) &
                         + (cvl-1) * n_tor_local + harm_active(csl)             ! index in the ELM matrix (full harmonic order)
                ilarge2 = ijA_position - 1 + (cj-1) * n_var * n_act + cl
                a_mat%irn(ilarge2) = index_large_i + cj
                a_mat%jcn(ilarge2) = index_large_k + cl
                thread_struct(omp_tid)%synch_buff((cj-1)*n_var*n_act+cl) = &
                  thread_struct(omp_tid)%synch_buff((cj-1)*n_var*n_act+cl) + thread_struct(omp_tid)%ELM(index_ij,index_kl)
              enddo
            enddo
          enddo
        enddo
        !$ tw2 = omp_get_wtime()
        !$omp critical
        !$ tw3 = omp_get_wtime()
        !$ call tm_add(omp_tid, 14, tw3 - tw2)
        a_mat%val(ijA_position : ijA_position + n_var*n_act*n_var*n_act - 1) = &
          a_mat%val(ijA_position : ijA_position + n_var*n_act*n_var*n_act - 1) +  &
          thread_struct(omp_tid)%synch_buff(1:n_var*n_act*n_var*n_act)
        !$ call tm_add(omp_tid, 15, omp_get_wtime() - tw3)
        !$omp end critical
        cycle
      endif

      thread_struct(omp_tid)%synch_buff(1:n_var*n_tor_local*n_var*n_tor_local) = 0.d0

      do j = 1, n_var * n_tor_local
        index_ij = n_tor_local * n_var * n_degrees * (i-1) + n_tor_local * n_var * (i_order-1) + j

        do l = 1, n_var * n_tor_local

          index_kl = n_tor_local * n_var * n_degrees * (k-1) +  n_tor_local * n_var * (k_order-1) + l

          ilarge2 = ijA_position - 1 + (j-1) * n_var * n_tor_local + l

          a_mat%irn(ilarge2) = index_large_i + j
          a_mat%jcn(ilarge2) = index_large_k + l

          thread_struct(omp_tid)%synch_buff((j-1)*n_var*n_tor_local+l) = &
            thread_struct(omp_tid)%synch_buff((j-1)*n_var*n_tor_local+l) + thread_struct(omp_tid)%ELM(index_ij,index_kl)
        enddo

      enddo

      if (.not. eliminate_boundary_dofs) then
        !$omp critical
        a_mat%val(ijA_position : ijA_position + n_var*n_tor_local*n_var*n_tor_local - 1) = &
          a_mat%val(ijA_position : ijA_position + n_var*n_tor_local*n_var*n_tor_local - 1) +  &
          thread_struct(omp_tid)%synch_buff(1:n_var*n_tor_local*n_var*n_tor_local)
        !$omp end critical
      else

        if (interior) then
          !$omp critical
          a_mat%val(ijA_position : ijA_position + n_var*n_tor_local*n_var*n_tor_local - 1) = &
            a_mat%val(ijA_position : ijA_position + n_var*n_tor_local*n_var*n_tor_local - 1) +  &
            thread_struct(omp_tid)%synch_buff(1:n_var*n_tor_local*n_var*n_tor_local)
          !$omp end critical
        else
          if ((i .eq. k) .and. ((i_order .eq. 1 .and. k_order .eq. 1) &
                                    .or. ((i_bnd_type .eq. 2 .or. i_bnd_type .eq. 3) .and. (i_order .eq. 3 .and. k_order .eq. 3)) &
                                    .or. ((i_bnd_type .eq. 1 .or. i_bnd_type .eq. 3) .and. (i_order .eq. 2 .and. k_order .eq. 2)))) then

            !$omp critical
            do j = 1, n_var * n_tor_local
              do l = 1, n_var * n_tor_local
                ilarge2 = ijA_position - 1 + (j-1) * n_var * n_tor_local + l
                if (j .eq. l) then
                  a_mat%val(ilarge2) = a_mat%val(ilarge2) &
                                     + elm_diagonal_average((j-1)/n_tor_local + 1)
                else
                  a_mat%val(ilarge2) = 0.d0
                endif
              enddo
            enddo
            !$omp end critical
          else if ((i_bnd .and. (i_order .eq. 1 .or. (i_bnd_type .eq. 2 .and. i_order .eq. 3) .or. &
                                                    (i_bnd_type .eq. 1 .and. i_order .eq. 2) .or. &
                                                    (i_bnd_type .eq. 3 .and. (i_order .eq. 2 .or. i_order .eq. 3)))) &
                    .or. (k_bnd .and. (k_order .eq. 1 .or. (k_bnd_type .eq. 2 .and. k_order .eq. 3) .or. &
                                                    (k_bnd_type .eq. 1 .and. k_order .eq. 2) .or. &
                                                    (k_bnd_type .eq. 3 .and. (k_order .eq. 2 .or. k_order .eq. 3))))) then
            !$omp critical
            a_mat%val(ijA_position : ijA_position + n_var*n_tor_local*n_var*n_tor_local - 1) = 0.d0
            !$omp end critical
          else
            !$omp critical
            a_mat%val(ijA_position : ijA_position + n_var*n_tor_local*n_var*n_tor_local - 1) = &
              a_mat%val(ijA_position : ijA_position + n_var*n_tor_local*n_var*n_tor_local - 1) +  &
              thread_struct(omp_tid)%synch_buff(1:n_var*n_tor_local*n_var*n_tor_local)
            !$omp end critical
          endif
        endif

      endif ! eliminate_boundary_dofs

    enddo ! n_degrees
  enddo ! n_vertex_max

end subroutine add_block_to_sp_matrix


#ifdef USE_PETSC
!> Assemble element blocks directly into PETSc MPIBAIJ matrix.
!! Called from add_to_a_mat for the direct PETSc assembly path.
subroutine add_block_to_petsc(index_node1, i, i_order, i_bnd, i_bnd_type, &
                               node_out, a_mat, omp_tid, n_tor_local, &
                               eliminate_boundary_dofs, elm_diagonal_average)
#include "petsc/finclude/petsc.h"
  use petsc
  use mod_parameters, only: n_var, n_degrees, n_vertex_max
  use data_structure, only: type_node, type_node_list, type_SP_MATRIX, thread_struct
  use nodes_elements
  implicit none

  integer, intent(in) :: index_node1, i, i_order
  logical, intent(in) :: i_bnd
  integer, intent(in) :: i_bnd_type
  integer, intent(in) :: node_out(:)
  type(type_SP_MATRIX), intent(inout) :: a_mat
  integer, intent(in) :: omp_tid, n_tor_local
  logical, intent(in) :: eliminate_boundary_dofs
  real*8, intent(in) :: elm_diagonal_average(:)   !< one per variable

  integer :: j, k, l, k_order
  integer :: knode
  integer :: index_node2
  integer :: index_ij, index_kl
  logical :: interior, k_bnd
  integer :: k_bnd_type
  integer :: block_size

  PetscInt :: idxm_petsc(1), idxn_petsc(1)
  PetscErrorCode :: ierr

  block_size = n_var * n_tor_local

  do k=1,n_vertex_max

    knode = node_out(k)
    k_bnd = .false.

    k_bnd_type = node_list%node(knode)%boundary
    if (k_bnd_type .ne. 0) k_bnd = .true.

    interior = .not. (i_bnd .or. k_bnd)

    do k_order = 1, n_degrees

      index_node2 = node_list%node(knode)%index(k_order)

      ! --- Compute synch_buff from element matrix
      thread_struct(omp_tid)%synch_buff(1:block_size*block_size) = 0.d0
      do j = 1, block_size
        index_ij = n_tor_local * n_var * n_degrees * (i-1) + n_tor_local * n_var * (i_order-1) + j
        do l = 1, block_size
          index_kl = n_tor_local * n_var * n_degrees * (k-1) + n_tor_local * n_var * (k_order-1) + l
          thread_struct(omp_tid)%synch_buff((j-1)*block_size+l) = &
            thread_struct(omp_tid)%synch_buff((j-1)*block_size+l) + thread_struct(omp_tid)%ELM(index_ij,index_kl)
        enddo
      enddo

      ! --- Apply eliminate_boundary_dofs logic
      if (eliminate_boundary_dofs .and. .not. interior) then
        if ((i .eq. k) .and. ((i_order .eq. 1 .and. k_order .eq. 1) &
                                  .or. ((i_bnd_type .eq. 2 .or. i_bnd_type .eq. 3) .and. (i_order .eq. 3 .and. k_order .eq. 3)) &
                                  .or. ((i_bnd_type .eq. 1 .or. i_bnd_type .eq. 3) .and. (i_order .eq. 2 .and. k_order .eq. 2)))) then
          ! Diagonal boundary block: replace with diagonal of elm_diagonal_average
          thread_struct(omp_tid)%synch_buff(1:block_size*block_size) = 0.d0
          do j = 1, block_size
            thread_struct(omp_tid)%synch_buff((j-1)*block_size+j) = &
              elm_diagonal_average((j-1)/n_tor_local + 1)
          enddo
        else if ((i_bnd .and. (i_order .eq. 1 .or. (i_bnd_type .eq. 2 .and. i_order .eq. 3) .or. &
                                                  (i_bnd_type .eq. 1 .and. i_order .eq. 2) .or. &
                                                  (i_bnd_type .eq. 3 .and. (i_order .eq. 2 .or. i_order .eq. 3)))) &
                  .or. (k_bnd .and. (k_order .eq. 1 .or. (k_bnd_type .eq. 2 .and. k_order .eq. 3) .or. &
                                                  (k_bnd_type .eq. 1 .and. k_order .eq. 2) .or. &
                                                  (k_bnd_type .eq. 3 .and. (k_order .eq. 2 .or. k_order .eq. 3))))) then
          ! Off-diagonal boundary block: skip (matrix already zeroed by MatZeroEntries)
          cycle
        endif
      endif

      ! --- Insert block into PETSc matrix
      idxm_petsc(1) = index_node1 - 1  ! 0-based block row
      idxn_petsc(1) = index_node2 - 1  ! 0-based block col
      !$omp critical
      PetscCallA(MatSetValuesBlocked(a_mat%petsc_A, 1, idxm_petsc, 1, idxn_petsc, thread_struct(omp_tid)%synch_buff, ADD_VALUES, ierr))
      !$omp end critical

    enddo ! n_degrees
  enddo ! n_vertex_max

end subroutine add_block_to_petsc
#endif


subroutine define_element_nodes(ielm, node_out, element, nodes, element_father, nodes_father, omp_tid, harmonic_matrix)
  use phys_module, only: refinement
  use mod_parameters, only: n_vertex_max
  use data_structure, only: type_element, type_node, thread_struct
  use mod_ch_nod_rhs_elm
  implicit none

  integer,            intent(in)          :: ielm
  integer,            intent(inout)       :: node_out(n_vertex_max)
  type(type_element), intent(in)          :: element
  type(type_node)   , intent(inout)       :: nodes(n_vertex_max)
  type(type_element), intent(in)          :: element_father
  type(type_node)   , intent(inout)       :: nodes_father(n_vertex_max)
  integer,            intent(in)          :: omp_tid
  logical,            intent(in)          :: harmonic_matrix

  integer :: iv, inode

  if (refinement .and. .not. harmonic_matrix) then   
    call ch_nod_rhs_elm(ielm,element,nodes,element_father,nodes_father, &
            thread_struct(omp_tid)%ELM, thread_struct(omp_tid)%RHS,node_out)
  else
    do iv=1, n_vertex_max
      node_out(iv) = element%vertex(iv)   
    enddo 
  endif
end subroutine define_element_nodes


subroutine initialize_v_and_harm(i_v, i_harm, i_tor_min, i_tor_max, n_tor_local)
  use mod_parameters, only: n_var, n_degrees, n_tor
  implicit none

  integer, intent(out) :: i_v(n_var)
  integer, allocatable, intent(out) :: i_harm(:)
  integer, intent(in) :: i_tor_min, i_tor_max
  integer, intent(in) :: n_tor_local

  integer :: i

  ! --- Initialize i_v
  do i = 1, n_var
    i_v(i) = i
  enddo

  ! --- Initialize i_harm
    if (.not. allocated(i_harm)) allocate(i_harm(n_tor_local))
    do i = i_tor_min, i_tor_max
      i_harm(i) = i
    enddo
end subroutine initialize_v_and_harm


!> check if matrix is row distributed
subroutine check_if_distributed(a_mat)
  use mpi_mod
  use data_structure, only: type_SP_MATRIX
  use mod_integer_types
  implicit none

  type(type_SP_MATRIX)  :: a_mat
  integer(kind=int_all) :: nloc, nglob
  integer               :: ierr

  nloc = a_mat%irn(a_mat%nnz) - a_mat%irn(1) + 1
  call MPI_AllReduce(nloc,nglob,1,MPI_INTEGER_ALL,MPI_SUM,a_mat%comm,ierr)
  a_mat%row_distributed = (nglob.eq.a_mat%ng)

  nloc = a_mat%jcn(a_mat%nnz) - a_mat%jcn(1) + 1
  call MPI_AllReduce(nloc,nglob,1,MPI_INTEGER_ALL,MPI_SUM,a_mat%comm,ierr)
  a_mat%col_distributed = (nglob.eq.a_mat%ng)

end subroutine check_if_distributed


!> Add dt (s) to timer k of the buffers of thread tid. The list of timers is in report_construct_timers.
subroutine tm_add(tid, k, dt)
  use data_structure, only: thread_struct
  implicit none
  integer, intent(in) :: tid, k
  real*8,  intent(in) :: dt

  thread_struct(tid)%tm(k) = thread_struct(tid)%tm(k) + dt

end subroutine tm_add

!> Linear run by mode family: print where the time of construct_matrix goes. The accumulators of the threads (filled in
!! element_matrix_fft and in the element loop of construct_matrix) are reduced to their maximum over the threads of a rank, then over
!! the ranks of comm. Collective: every rank of comm has to call it, also ranks without local elements. Output on rank 0 of comm.
!! Timers (s): the first group is wall-clock of the whole call (alloc, loop, bc, end, total); the rest is the maximum over
!! threads and ranks, summed over the elements of the thread. fft is contained in blk_*; blk_* and rhs_asm are the harmonic
!! assembly; zero includes the final 0.5 scaling; loop is the whole element loop body, so loop - build - axis - scatter is the
!! rest (node copies etc.); scatter contains locate, crit_wait and crit_exec. rss and hwm are the current and peak resident set (GB).
subroutine report_construct_timers(comm, n_local_elms, tw_start, tw_alloc, tw_loop, tw_bc)
  use data_structure, only: thread_struct, nbthreads
  use tr_module,      only: get_memory_inkb
  use mpi_mod
  !$ use omp_lib
  implicit none

  integer, intent(in) :: comm, n_local_elms
  real*8,  intent(in) :: tw_start, tw_alloc, tw_loop, tw_bc

  integer, parameter  :: ntm = 16
  real*8              :: tmax(ntm+6), tglob(ntm+6), twall
  integer             :: i, ierr, crank, nloc(2), nglob(2)

  twall = 0.d0
  !$ twall = omp_get_wtime() - tw_start

  do i = 1, ntm
    tmax(i) = maxval(thread_struct(1:nbthreads)%tm(i))
  enddo
  tmax(ntm+1) = tw_alloc
  tmax(ntm+2) = tw_loop
  tmax(ntm+3) = tw_bc
  tmax(ntm+4) = twall - tw_alloc - tw_loop - tw_bc
  tmax(ntm+5) = real(get_memory_inkb("VmRSS"),8) / (1024.d0*1024.d0)
  tmax(ntm+6) = real(get_memory_inkb("VmHWM"),8) / (1024.d0*1024.d0)
  call MPI_Reduce(tmax, tglob, ntm+6, MPI_DOUBLE_PRECISION, MPI_MAX, 0, comm, ierr)

  nloc(1) = n_local_elms
  nloc(2) = maxval(thread_struct(1:nbthreads)%tm_n)
  call MPI_Reduce(nloc, nglob, 2, MPI_INTEGER, MPI_MAX, 0, comm, ierr)

  call MPI_Comm_rank(comm, crank, ierr)
  if (crank == 0) then
    write(*,'(A,*(1X,A,"=",F0.3))') ' linear_run timers:', 'alloc', tglob(ntm+1), 'loop', tglob(ntm+2), 'bc', tglob(ntm+3),     &
      'end', tglob(ntm+4), 'build', tglob(10), 'zero', tglob(1), 'setup', tglob(2), 'planes', tglob(3), 'fft', tglob(4),       &
      'blk_p', tglob(5), 'blk_n', tglob(6), 'blk_k', tglob(7), 'blk_kn', tglob(8), 'rhs_asm', tglob(9), 'axis', tglob(11),     &
      'scatter', tglob(12), 'locate', tglob(13), 'crit_wait', tglob(14), 'crit_exec', tglob(15), 'thread_loop', tglob(16),     &
      'rss_GB', tglob(ntm+5), 'hwm_GB', tglob(ntm+6), 'max_elem_rank', real(nglob(1),8), 'max_elem_thread', real(nglob(2),8), &
      'threads', real(nbthreads,8)
  endif

end subroutine report_construct_timers

!> Debug (linear_check_sums): sums over the global matrix and RHS that two builds of the same case can be compared by.
!! Matrix: sum of |a|, sum of a and number of entries whose row and column harmonic are both evolving (not harmonic_frozen). Without a
!! linear run every harmonic is evolving. RHS: the same over the evolving rows (taken on rank 0, where it is complete).
!! The harmonic of a row or column is mod(index-1, n_tor)+1, or the slot mod(index-1, n_act)+1 with linear_compact (then the whole
!! system is evolving). Collective: every rank of comm has to call it. Output on rank 0 of comm.
subroutine report_linear_check_sums(comm, a_mat, rhs_vec)
  use data_structure, only: type_SP_MATRIX, type_RHS
  use phys_module,    only: harmonic_frozen
  use mod_parameters, only: n_tor
  use mod_linear_compact, only: compact_on, n_act
  use mod_integer_types
  use mpi_mod
  implicit none

  integer,              intent(in) :: comm
  type(type_SP_MATRIX), intent(in) :: a_mat
  type(type_RHS),       intent(in) :: rhs_vec

  integer, parameter    :: ns = 6
  real*8                :: s(ns), sg(ns), sa, ss, na
  logical               :: ev(n_tor)
  integer(kind=int_all) :: i, nh, hr, hc
  integer               :: crank, ierr

  nh = n_tor
  ev = .not. harmonic_frozen
  if (compact_on) then
    nh = n_act
    ev = .true.
  endif

  sa = 0.d0; ss = 0.d0; na = 0.d0
  !$omp parallel do default(none) shared(a_mat, ev, nh) private(i, hr, hc) reduction(+:sa, ss, na) schedule(static)
  do i = 1, a_mat%nnz
    if ((a_mat%irn(i) .le. 0) .or. (a_mat%jcn(i) .le. 0)) cycle
    hr = mod(a_mat%irn(i) - Int1, nh) + 1
    hc = mod(a_mat%jcn(i) - Int1, nh) + 1
    if (ev(hr) .and. ev(hc)) then
      sa = sa + abs(a_mat%val(i))
      ss = ss + a_mat%val(i)
      na = na + 1.d0
    endif
  enddo
  !$omp end parallel do

  s = 0.d0
  s(1) = sa; s(2) = ss; s(3) = na
  call MPI_Comm_rank(comm, crank, ierr)
  if (crank .eq. 0) then
    do i = 1, rhs_vec%n
      hr = mod(i - Int1, nh) + 1
      if (ev(hr)) then
        s(4) = s(4) + abs(rhs_vec%val(i))
        s(5) = s(5) + rhs_vec%val(i)
        s(6) = s(6) + 1.d0
      endif
    enddo
  endif
  call MPI_Reduce(s, sg, ns, MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm, ierr)

  if (crank .eq. 0) then
    write(*,'(A,I0)') ' linear_check_sums: harmonics per block ', nh
    write(*,'(A,ES26.18,A,ES26.18,A,F0.0)') ' linear_check_sums: matrix sum|a| =', sg(1), '  sum a =', sg(2), '  entries =', sg(3)
    write(*,'(A,ES26.18,A,ES26.18,A,F0.0)') ' linear_check_sums: rhs    sum|b| =', sg(4), '  sum b =', sg(5), '  rows    =', sg(6)
  endif

end subroutine report_linear_check_sums

end module construct_matrix_mod

