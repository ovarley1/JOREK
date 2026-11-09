module mod_assembly
implicit none
!> This module provides routines to factorize code in boundary conditions from each models

contains
  !>
  !! Subroutine: boundary_conditions_add_one_entry
  !!
  !! Add one entry to the product and/or harminic matrix if index_node is local.
  !!
  !! @param index_node  row node index
  !! @param k           row var index
  !! @param in          row tor index
  !! @param index_node2 col node index
  !! @param k2          col var index
  !! @param in2         col tor index
  !! @param zbig        value
  !! @param use_murge   Use murge interface.
  !! @param use_murge_element Murge interface with elementary matrices.
  !! @param index_min   Minimal local element index
  !! @param index_max   Maximal local element index
  !!
  subroutine boundary_conditions_add_one_entry( &
       &   index_node,  k,  in,                 &
       &   index_node2, k2, in2,                &
       &   zbig, index_min, index_max, a_mat)
    use mod_parameters
    use mod_locate_irn_jcn
    use mod_integer_types
    use mod_linear_compact, only: compact_on, n_act, harm_slot
    use data_structure, only: type_SP_MATRIX

    integer,               intent(in)                 :: k,  in
    integer,               intent(in)                 :: k2, in2
    integer,               intent(in)                 :: index_node
    integer,               intent(in)                 :: index_node2
    real*8,                intent(in)                 :: zbig
    integer,               intent(in)                 :: index_min, index_max
    type(type_SP_MATRIX)                              :: a_mat
    
    logical                                           :: is_local
    integer(kind=int_all)                             :: ija_position, ilarge_vp
    integer                                           :: n_tor_local, islot, islot2

    if ( (k==0) .or. (k2==0) ) return ! ignore calls for model family extensions not in use (variable number zero)

    n_tor_local = a_mat%i_tor_max - a_mat%i_tor_min +1
    islot  = in  - a_mat%i_tor_min + 1  ! position of the harmonic in the block of the system
    islot2 = in2 - a_mat%i_tor_min + 1
    if ( compact_on ) then              ! linear_compact: the block holds the active harmonics only, a frozen one has no row or column
      n_tor_local = n_act
      islot  = harm_slot(in)
      islot2 = harm_slot(in2)
      if ( (islot == 0) .or. (islot2 == 0) ) return
    endif
    if ((index_node .ge. index_min) .and. (index_node .le. index_max)) then

       call locate_irn_jcn(index_node,index_node2,index_min,index_max,ijA_position,a_mat)
                             
       !-------- index dans A_mat
       ilarge_vp  = ijA_position  - 1 + ((k-1)*n_tor_local + islot - 1) * n_var*n_tor_local + (k2-1)*n_tor_local + islot2
                               
                             
       a_mat%irn(ilarge_vp) =  n_tor_local * n_var * (index_node -1) + (k -1)*n_tor_local + islot
       a_mat%jcn(ilarge_vp) =  n_tor_local * n_var * (index_node2-1) + (k2-1)*n_tor_local + islot2
       a_mat%val(ilarge_vp) = ZBIG
    endif
  end subroutine boundary_conditions_add_one_entry

  !>
  !! Subroutine: boundary_conditions_add_RHS
  !!
  !! Add one entry to the product and/or harminic matrix if index_node is local.
  !!
  !! @param index_node  row node index
  !! @param k           row var index
  !! @param in          row tor index
  !! @param index_node2 col node index
  !! @param k2          col var index
  !! @param in2         col tor index
  !! @param zbig        value
  !! @param index_min   Minimal local element index
  !! @param index_max   Maximal local element index
  !!
  subroutine boundary_conditions_add_RHS( &
       &   index_node, k, in,             &
       &   index_min, index_max,          &
       &   rhs_loc,  val,                 &
       &   i_tor_min, i_tor_max)
    use mod_parameters
    use mod_integer_types
    use mod_linear_compact, only: compact_on, n_act, harm_slot
    integer,               intent(in)    :: k,  in
    integer,               intent(in)    :: index_node
    integer,               intent(in)    :: index_min, index_max
    integer,               intent(in)    :: i_tor_min, i_tor_max 
    real*8,                intent(in)    :: val
    real*8,                intent(inOUT) :: rhs_loc(*) 
    integer                              :: n_tor_local, islot
    logical                              :: is_local

    if ( (k==0) ) return ! ignore calls for model family extensions not in use (variable number zero)

    n_tor_local = i_tor_max - i_tor_min +1
    islot       = in - i_tor_min + 1
    if ( compact_on ) then ! linear_compact: the block of the system holds the active harmonics only (not reached by model 183)
      n_tor_local = n_act
      islot       = harm_slot(in)
      if ( islot == 0 ) return
    endif
    if ((index_node .ge. index_min) .and. (index_node .le. index_max)) then
       RHS_loc(n_tor_local*n_var * (index_node-1) + (k-1)*n_tor_local + islot) = val
    endif
  end subroutine boundary_conditions_add_RHS
end module mod_assembly
