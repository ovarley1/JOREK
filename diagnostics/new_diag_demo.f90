!> Demonstration of the diagnostic framework mod_position / mod_expression / mod_four_filter / mod_straight_field_line / ...
program demo
  
  use mod_parameters
  use data_structure
  use phys_module
  use mod_boundary
  use mod_new_diag
  use basis_at_gaussian
  use mod_import_restart
  use equil_info
  
  implicit none
  
  type(type_node_list),         pointer :: node_list
  type(type_element_list),      pointer :: element_list
  type (type_bnd_element_list), pointer :: bnd_elm_list
  type (type_bnd_node_list),    pointer :: bnd_node_list
  type(t_pol_pos_list) :: pol_pos_list
  type(t_tor_pos_list) :: tor_pos_list
  type(t_four_filter)  :: filter
  type(t_expr_list)    :: expr_list
  integer :: my_id, ierr, k_tor, i, j, k, n(4)
  real*8, allocatable :: result(:,:,:,:), res0d(:), res1d(:,:), res2d(:,:,:)
  complex*16, allocatable :: cp(:,:,:,:)

  ! --- Normal initialization
  allocate(node_list)
  allocate(element_list)
  allocate(bnd_elm_list)
  allocate(bnd_node_list)
  write(*,*) 'allocated necessary data structures'
  my_id = 0
  call initialise_parameters(my_id, "__NO_FILENAME__")
  call det_modes()
  call import_restart(node_list, element_list, 'jorek_restart',  rst_format, ierr, .true.)
  call initialise_basis()
  call boundary_from_grid(node_list, element_list, bnd_node_list, bnd_elm_list, .false.)
  
  ! --- Initialize the plasma equilibrium data structure
  call update_equil_state(my_id,node_list, element_list, bnd_elm_list, xpoint, xcase)
  call print_equil_state(.false.)
  
  ! --- Initialize the new_diag framework and print some information (.true.)
  call init_new_diag(.true.)

  if(my_id ==0 ) call export_nemec(node_list, element_list, xpoint, xcase)
  
  write(*,*) 'nemec finished.'

  stop
  ! **************** USING THE UNDERLYING ROUTINES TO PERFORM MORE FLEXIBLE TASKS ******************
  
  
  
  
  
  
!  ! --- How to select expressions, some examples:
!  expr_list = exprs_all
!  expr_list = exprs((/'Psi ', 'B_R ', 'xjac', 'T   ', 'rho ', 'zj  '/), 6)
!  
!  
!  
!  ! --- Evaluate several expressions at one single position.
!  call create_pol_pos(pol_pos_list, ierr, node_list, element_list, ES, R=3., Z=0.1)
!  call create_tor_pos(tor_pos_list, ierr, phi=0.)
!  expr_list = exprs((/'B_R ', 'xjac', 'T   ', 'rho ', 'zj  '/), 5)
!  call eval_expr(ES, JOREK_UNITS, expr_list, pol_pos_list, tor_pos_list, result, ierr)
!  
!  
!  
!  ! --- Print results in two different ways to the screen
!  call reduce_result_to_0d(ierr, result, res0d, 1, 1, 1)
!  call write_ascii_0d(ierr, ES, expr_list, res0d, FORM_TABLE, header=.true.)
!  call write_ascii_0d(ierr, ES, expr_list, res0d, FORM_LIST)
!  
!  
!  
!  ! --- Evaluate several expressions on the outboard midplane and write to file.
!  expr_list = exprs((/'R    ', 'Z    ', 'Psi_N', 'Psi  ', 'theta', 'x    ', 'y    ', &
!       'phi  ', 'B_R  ', 'xjac ', 'T    ', 'rho  ', 'zj   ', 'omega', 'u    '/), 15)
!  call create_pol_pos(pol_pos_list, ierr, node_list, element_list, ES,                    &
!      Rstart=ES%R_axis, Rend=ES%R_midpl(2)-1.d-3, Zstart=ES%Z_axis,     &
!      Zend=ES%Z_axis, n=500)
!  call eval_expr(ES, JOREK_UNITS, expr_list, pol_pos_list, tor_pos(phi=0.), result, ierr)
!  call reduce_result_to_1d(ierr, result, res1d, i1=1, i2=1)
!  call write_ascii_1d(ierr, ES, expr_list, res1d, FORM_TABLE, header=.true.,              &
!    filename='midplane_profiles.dat', append=.false., comment='Various profiles')
!  
!  
!  
!  ! --- Evaluate expressions along toroidal direction.
!  call eval_expr(ES, JOREK_UNITS, expr_list,                                              &
!    pol_pos(node_list,element_list,ES,R=3.,Z=0.1), tor_pos(nphi=128), result, ierr)
!  
!  
!  
!  ! --- Evaluate expressions on flux surfaces using straight field line theta.
!  call create_pol_pos(pol_pos_list, ierr, node_list, element_list, ES, PsiNmin=0.01,      &
!    PsiNmax=0.99, nPsiN=16, nTht=64)
!  call create_tor_pos(tor_pos_list, ierr, nphi=8)
!  call eval_expr(ES, JOREK_UNITS, expr_list, pol_pos_list, tor_pos_list, result, ierr)
!  
!  
!  
!  ! --- Output 2D VTK files
!  call reduce_result_to_2d(ierr, result, res2d, i1=1)
!  call write_vtk_2d(ierr, expr_list, res2d, 'testA.vtk', (/1,2/), close1=.true.) ! vs R,Z
!  res2d(:,:,5)=res2d(:,:,5)/(2.d0*PI)
!  call write_vtk_2d(ierr, expr_list, res2d, 'testB.vtk', (/3,5/)) ! vs PsiN,theta
!  call write_vtk_2d(ierr, expr_list, res2d, 'testC.vtk', (/6,2/)) ! vs x,z
!  call reduce_result_to_2d(ierr, result, res2d, i3=16)
!  call write_vtk_2d(ierr, expr_list, res2d, 'testC.vtk', (/5,8/)) ! vs theta,phi
!  
!  
!  
!  ! --- Check that 2D Fourier transform forward and backward recovers the original data
!  n(:) = (/ size(result,1), size(result,2), size(result,3), size(result,4) /)
!  call perform_four_trafo(result, POLTOR_TRAFO, FORWARD_TRAFO, ierr)
!  call perform_four_trafo(result, POLTOR_TRAFO, BACKWARD_TRAFO, ierr)
!  
!  
!  
!  ! --- Check that 1D toroidal Fourier transform forward and backward recovers the original data
!  call perform_four_trafo(result, TOROIDAL_TRAFO, FORWARD_TRAFO, ierr)
!  call perform_four_trafo(result, TOROIDAL_TRAFO, BACKWARD_TRAFO, ierr)
!  
!  
!  
!  ! --- Check that 1D poloidal Fourier transform forward and backward recovers the original data
!  call perform_pol_trafo(result, FORWARD_TRAFO, ierr)
!  call perform_pol_trafo(result, BACKWARD_TRAFO, ierr)
!  
!  
!  
!  ! --- Apply Fourier filter: Keep harmonics (m,n) = (0...1,0...5) and (m,n)=(<any>,16)
!  call perform_four_trafo(result, POLTOR_TRAFO, FORWARD_TRAFO, ierr)
!  call init_four_filter(filter)
!  call filter_add(filter, ierr, m_start=0, m_end=1, n_end=5)
!  call filter_add(filter, ierr, n=16)
!  call print_filter(filter)
!  call apply_four_filter(result, filter, ierr)
!  call perform_four_trafo(result, POLTOR_TRAFO, BACKWARD_TRAFO, ierr)
!  
!  
!  
!  ! --- Or simple filtering is possible as a single call:
!  call transform_and_filter(result, simple_filter(m=0,n=0), ierr) ! Keep only (m,n)=(0,0)
!  
!  
!  
!  ! --- Write result to a file = poloidally and toroidally averaged profiles
!  call reduce_result_to_1d(ierr, result, res1d, i1=1, i2=1)
!  call write_ascii_1d(ierr, ES, expr_list, res1d, FORM_TABLE, header=.true.,              &
!    filename='average_profiles_1.dat')
!  
!  
!  
!  ! --- Toroidally averaged expressions on the outboard midplane.
!  call create_pol_pos(pol_pos_list, ierr, node_list, element_list, ES,                    &
!      Rstart=ES%R_axis+1d-3, Rend=ES%R_midpl(2)-1d-3, Zstart=ES%Z_axis, &
!      Zend=ES%Z_axis, n=200)
!  call eval_expr(ES, JOREK_UNITS, expr_list, pol_pos_list, tor_pos(nphi=16), result, ierr)
!  call transform_and_filter(result, simple_filter(n=0), ierr)
!  call reduce_result_to_1d(ierr, result, res1d, i1=1, i2=1)
!  call write_ascii_1d(ierr, ES, expr_list, res1d, FORM_TABLE, header=.true.,              &
!    filename='toroidally_averaged_midplane_profiles.dat', append=.false.)
  
end program demo
