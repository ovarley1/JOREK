subroutine update_deltas(node_list, deltas)
!---------------------------------------------------------------------
! subroutine to create a local list of delta values
!---------------------------------------------------------------------
use tr_module
use data_structure, only: type_node_list, type_RHS
use mod_parameters, only: n_tor, n_var, n_degrees
use mod_linear_compact, only: compact_on, n_act, harm_slot
!use global_distributed_matrix

implicit none
integer               :: i, j, k, in, index, index_node, n_sys, index_total
type(type_node_list)  :: node_list
type(type_RHS)        :: deltas

if (.not. associated(deltas%val)) then
  n_sys = node_list%n_dof
  if (compact_on) then ! linear_compact: n_act harmonics per node block in the vectors of the system (n_dof of the nodes is unchanged)
    index_total = 0
    do i = 1, node_list%n_nodes
      index_total = max(index_total, maxval(node_list%node(i)%index))
    enddo
    n_sys = index_total * n_act * n_var
  endif
  call tr_allocatep(deltas%val,1,n_sys,"deltas",CAT_FEM)
  deltas%n = n_sys
  deltas%val(1:deltas%n) = 0.d0
endif

do i = 1, node_list%n_nodes
 if ((.not. node_list%node(i)%constrained)) then
  do j=1,n_degrees

    index_node = node_list%node(i)%index(j)

    do k=1,n_var

      do in=1,n_tor

        if (compact_on) then
          if (harm_slot(in) == 0) cycle ! frozen harmonic: not in the system
          index = n_act*n_var * (index_node - 1) + n_act*(k-1) + harm_slot(in)
        else
          index = n_tor*n_var * (index_node - 1) + n_tor*(k-1) + in
        endif

        deltas%val(index) = node_list%node(i)%deltas(in,j,k)

      enddo

    enddo

  enddo
 endif
enddo

return
end
