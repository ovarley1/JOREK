subroutine initial_conditions(my_id,node_list,element_list,bnd_node_list, bnd_elm_list, xpoint2, xcase2)

use constants
use data_structure
use phys_module
use mod_poiss
use equil_info
use mod_sources

implicit none

type (type_node_list)    :: node_list
type (type_element_list) :: element_list
type (type_surface_list) :: surface_list
type (type_bnd_node_list)    :: bnd_node_list
type (type_bnd_element_list) :: bnd_elm_list

integer    :: my_id, i, in, mm, i_elm, ifail, xcase2
real*8     :: amplitude, psi, psi_n, theta
real*8     :: zn,  dn_dpsi,  dn_dpsi2,  dn_dz,  dn_dz2,  dn_dpsi_dz,  dn_dpsi3,  dn_dpsi2_dz,  dn_dpsi_dz2
real*8     :: zrn, drn_dpsi, drn_dpsi2, drn_dz, drn_dz2, drn_dpsi_dz, drn_dpsi3, drn_dpsi2_dz, drn_dpsi_dz2
real*8     :: zT,  dT_dpsi,  dT_dpsi2,  dT_dz,  dT_dz2,  dT_dpsi_dz,  dT_dpsi3,  dT_dpsi2_dz,  dT_dpsi_dz2
real*8     :: zTi, dTi_dpsi, dTi_dpsi2, dTi_dz, dTi_dz2, dTi_dpsi_dz, dTi_dpsi3, dTi_dpsi2_dz, dTi_dpsi_dz2
real*8     :: zTe, dTe_dpsi, dTe_dpsi2, dTe_dz, dTe_dz2, dTe_dpsi_dz, dTe_dpsi3, dTe_dpsi2_dz, dTe_dpsi_dz2
real*8     :: zFFprime,dFFprime_dpsi,dFFprime_dz, dFFprime_dpsi_dz, dFFprime_dz2, dFFprime_dpsi2
real*8     :: R, Z, BigR, T0, BigR_s, T0_s
real*8     :: zjz, dj_dpsi, dj_dR, dj_dZ, dj_dR_dZ, dj_dR_DR, dj_dZ_dZ, dj_dpsi2, dj_dR_dpsi, dj_dZ_dpsi
real*8     :: zp, dp_dpsi, dp_dpsi2, dp_dz, dp_dz2, dp_dpsi_dz, P_ss, P_st, P_tt, R_out,Z_out,s_out,t_out
real*8     :: ps0_s, ps0_t, p_s, p_t, zj0_s, zj0_t,R_s, R_t, ps0_x, ps0_y, Z_s, Z_t, xjac, direction, Btot
logical    :: xpoint2
real*8     :: zV, dV_dpsi, dV_dpsi2, dV_dz, dV_dz2, dV_dpsi_dz, dV_dpsi3, dV_dpsi2_dz, dV_dpsi_dz2
real*8     :: Omega, dOmega_dpsi, dOmega_dz, dOmega_dpsi2, dOmega_dz2, dOmega_dpsi_dz
! --- u0rho patch (021_castor_comparison, 2026-10-20): table of G(psi) = int_{psi_bnd}^{psi} (dp/dpsi')/rho dpsi'
integer, parameter   :: n_u0tab = 20001
integer              :: k_u0
real*8, allocatable  :: u0tab_psi(:), u0tab_f(:), u0tab_G(:)
real*8               :: u0_psi_lo, u0_psi_hi, u0_G_bnd, u0_wgt, u0_zG, u0_zf, u0_df
if (my_id .eq. 0) then
  write(*,*) '***************************************'
  write(*,*) '*      initial conditions  (600)      *'
  write(*,*) '***************************************'
endif

! --- This old projection method should be replaced by a direct solve on the elements
! --- It is much cleaner and more generic
if ( (my_id .eq. 0) .and. (n_order .le. 3) ) then

  ! --- u0rho patch (021_castor_comparison, 2026-10-20). Equilibrium flow with the ion perpendicular flow exactly zero.
  ! --- Model 600's ion perpendicular flow is rho*grad(u) + 2*tauIC*grad(Pi_i) with Pi_i = rho*T/2 (mod_elt_matrix_fft.f90
  ! --- l.1632, 1762, 2455), so it vanishes for du/dpsi = -tauIC*(dp/dpsi)/rho. The stock initialisation, u0 = -tauIC*p,
  ! --- is exact only where rho = 1 (it leaves R = 1 - rho). Here u0(psi) = -tauIC*G(psi), G(psi) = int_{psi_bnd}^{psi}
  ! --- (dp/dpsi')/rho dpsi' (zero on the boundary), tabulated once in psi from the same profile routines; w0 then comes
  ! --- from the stock inverse Poisson solve further down. Valid only for profiles of psi alone (no X-point).
  if (tauIC .ne. 0.d0) then
    if (xpoint2) then
      write(*,*) 'u0rho patch: profiles depend on Z with an X-point, G(psi) is not defined; stopping'
      stop
    endif
    allocate(u0tab_psi(n_u0tab), u0tab_f(n_u0tab), u0tab_G(n_u0tab))
    u0_psi_lo = min(ES%psi_axis, ES%psi_bnd)
    u0_psi_hi = max(ES%psi_axis, ES%psi_bnd)
    do i = 1, node_list%n_nodes
      u0_psi_lo = min(u0_psi_lo, node_list%node(i)%values(1,1,var_psi))
      u0_psi_hi = max(u0_psi_hi, node_list%node(i)%values(1,1,var_psi))
    enddo
    do k_u0 = 1, n_u0tab
      u0tab_psi(k_u0) = u0_psi_lo + (u0_psi_hi - u0_psi_lo) * dble(k_u0 - 1) / dble(n_u0tab - 1)
      psi = u0tab_psi(k_u0)
      call density(xpoint2, xcase2, 0.d0, ES%Z_xpoint, psi, ES%psi_axis, ES%psi_bnd, zn, dn_dpsi, dn_dz, dn_dpsi2, dn_dz2, &
                   dn_dpsi_dz, dn_dpsi3, dn_dpsi_dz2, dn_dpsi2_dz)
      if (zn .le. 0.d0) then
        write(*,*) 'u0rho patch: density <= 0 at psi = ', psi, '; stopping'
        stop
      endif
      if ( with_TiTe ) then
        call temperature_i(xpoint2, xcase2, 0.d0, ES%Z_xpoint, psi, ES%psi_axis, ES%psi_bnd, zTi, dTi_dpsi, dTi_dz, dTi_dpsi2, &
                           dTi_dz2, dTi_dpsi_dz, dTi_dpsi3, dTi_dpsi_dz2, dTi_dpsi2_dz)
        call temperature_e(xpoint2, xcase2, 0.d0, ES%Z_xpoint, psi, ES%psi_axis, ES%psi_bnd, zTe, dTe_dpsi, dTe_dz, dTe_dpsi2, &
                           dTe_dz2, dTe_dpsi_dz, dTe_dpsi3, dTe_dpsi_dz2, dTe_dpsi2_dz)
        u0tab_f(k_u0) = (zn * (dTi_dpsi + dTe_dpsi) + dn_dpsi * (zTi + zTe)) / zn
      else
        call temperature(xpoint2, xcase2, 0.d0, ES%Z_xpoint, psi, ES%psi_axis, ES%psi_bnd, zT, dT_dpsi, dT_dz, dT_dpsi2, &
                         dT_dz2, dT_dpsi_dz, dT_dpsi3, dT_dpsi_dz2, dT_dpsi2_dz)
        u0tab_f(k_u0) = (zn * dT_dpsi + dn_dpsi * zT) / zn          ! p = rho*T, as zp below (Ti = Te = T/2)
      endif
    enddo
    u0tab_G(1) = 0.d0
    do k_u0 = 2, n_u0tab
      u0tab_G(k_u0) = u0tab_G(k_u0-1) + 0.5d0 * (u0tab_f(k_u0) + u0tab_f(k_u0-1)) * (u0tab_psi(k_u0) - u0tab_psi(k_u0-1))
    enddo
    u0_wgt   = (ES%psi_bnd - u0_psi_lo) / (u0_psi_hi - u0_psi_lo) * dble(n_u0tab - 1)
    k_u0     = min(max(int(u0_wgt) + 1, 1), n_u0tab - 1)
    u0_wgt   = u0_wgt - dble(k_u0 - 1)
    u0_G_bnd = (1.d0 - u0_wgt) * u0tab_G(k_u0) + u0_wgt * u0tab_G(k_u0+1)
    u0tab_G  = u0tab_G - u0_G_bnd
    write(*,'(A,ES13.5,A,2ES13.5,A,I6,A)') ' u0rho patch: u0 = -tauIC*int(dp/dpsi/rho), tauIC =', tauIC, ', psi from', &
                                           u0_psi_lo, u0_psi_hi, ',', n_u0tab, ' points'
  endif

  do i=1,node_list%n_nodes

    psi = node_list%node(i)%values(1,1,var_psi)
    R   = node_list%node(i)%x(1,1,1)
    Z   = node_list%node(i)%x(1,1,2)
   

    call density(  xpoint2, xcase2, Z, ES%Z_xpoint, psi,ES%psi_axis,ES%psi_bnd,zn,dn_dpsi,dn_dz,dn_dpsi2,dn_dz2,          &
                                                               dn_dpsi_dz,dn_dpsi3,dn_dpsi_dz2, dn_dpsi2_dz)

    if ( with_TiTe ) then
      call temperature_i(xpoint2, xcase2, Z, ES%Z_xpoint, psi,ES%psi_axis,ES%psi_bnd,zTi,dTi_dpsi,dTi_dz,dTi_dpsi2,dTi_dz2, &
                                                                 dTi_dpsi_dz,dTi_dpsi3,dTi_dpsi_dz2, dTi_dpsi2_dz)

      call temperature_e(xpoint2, xcase2, Z, ES%Z_xpoint, psi,ES%psi_axis,ES%psi_bnd,zTe,dTe_dpsi,dTe_dz,dTe_dpsi2,dTe_dz2, &
                                                                 dTe_dpsi_dz,dTe_dpsi3,dTe_dpsi_dz2, dTe_dpsi2_dz)
    else
      call temperature(xpoint2,  xcase2, Z, ES%Z_xpoint, psi,ES%psi_axis,ES%psi_bnd,zT, dT_dpsi, dT_dz, dT_dpsi2, dT_dz2, &
                                                                 dT_dpsi_dz, dT_dpsi3, dT_dpsi_dz2,  dT_dpsi2_dz)
      zTi          = zT          / 2.d0
      dTi_dpsi     = dT_dpsi     / 2.d0
      dTi_dz       = dT_dz       / 2.d0
      dTi_dpsi2    = dT_dpsi2    / 2.d0
      dTi_dz2      = dT_dz2      / 2.d0
      dTi_dpsi_dz  = dT_dpsi_dz  / 2.d0
      dTi_dpsi3    = dT_dpsi3    / 2.d0
      dTi_dpsi_dz2 = dT_dpsi_dz2 / 2.d0
      dTi_dpsi2_dz = dT_dpsi2_dz / 2.d0

      zTe          = zTi
      dTe_dpsi     = dTi_dpsi
      dTe_dz       = dTi_dz
      dTe_dpsi2    = dTi_dpsi2
      dTe_dz2      = dTi_dz2
      dTe_dpsi_dz  = dTi_dpsi_dz
      dTe_dpsi3    = dTi_dpsi3
      dTe_dpsi_dz2 = dTi_dpsi_dz2
      dTe_dpsi2_dz = dTi_dpsi2_dz
    end if

    call FFprime(   xpoint2, xcase2, Z, ES%Z_xpoint, psi,ES%psi_axis,ES%psi_bnd,zFFprime,dFFprime_dpsi,dFFprime_dz, &
                                                               dFFprime_dpsi2,dFFprime_dz2, dFFprime_dpsi_dz)
    if ( (abs(V_0) .ge. 1.d-19) .or. (num_rot) ) then
       if (normalized_velocity_profile) then
          call velocity(xpoint2, xcase2, Z, ES%Z_xpoint, psi,ES%psi_axis,ES%psi_bnd,zV,dV_dpsi,dV_dz,dV_dpsi2,dV_dz2, &
               dV_dpsi_dz,dV_dpsi3,dV_dpsi_dz2, dV_dpsi2_dz)
       else
          call velocity(xpoint2, xcase2, Z, ES%Z_xpoint, psi,ES%psi_axis,ES%psi_bnd,Omega,dOmega_dpsi,dOmega_dz,dOmega_dpsi2,dOmega_dz2, &
               dOmega_dpsi_dz,dV_dpsi3,dV_dpsi_dz2, dV_dpsi2_dz)
       endif
    endif

    zp         = zn * (zTi + zTe)
    dp_dpsi    = zn * (dTi_dpsi  + dTe_dpsi)  + dn_dpsi * (zTi + zTe)
    dp_dpsi2   = zn * (dTi_dpsi2 + dTe_dpsi2) + 2.d0 * dn_dpsi * (dTi_dpsi + dTe_dpsi) + dn_dpsi2 * (zTi + zTe)
    dp_dz      = zn * (dTi_dz    + dTe_dz)    + dn_dz * (zTi + zTe)
    dp_dz2     = zn * (dTi_dz2   + dTe_dz2)   + 2.d0 * dn_dz * (dTi_dz + dTe_dz) + dn_dz2 * (zTi + zTe)
    dp_dpsi_dz = zn * (dTi_dpsi_dz + dTe_dpsi_dz) + dn_dz * (dTi_dpsi + dTe_dpsi) + dn_dpsi * (dTi_dz + dTe_dz) + dn_dpsi_dz * (zTi + zTe)

    if (with_rho) then
      node_list%node(i)%values(1,1,var_rho) = zn
      node_list%node(i)%values(1,2,var_rho) = dn_dpsi    * node_list%node(i)%values(1,2,var_psi) + dn_dz * node_list%node(i)%x(1,2,2)
      node_list%node(i)%values(1,3,var_rho) = dn_dpsi    * node_list%node(i)%values(1,3,var_psi) + dn_dz * node_list%node(i)%x(1,3,2)
      node_list%node(i)%values(1,4,var_rho) = dn_dpsi    * node_list%node(i)%values(1,4,var_psi) + dn_dz * node_list%node(i)%x(1,4,2) &
                                      + dn_dpsi2   * node_list%node(i)%values(1,2,var_psi) * node_list%node(i)%values(1,3,var_psi)  &
                                      + dn_dz2     * node_list%node(i)%x(1,2,2)              * node_list%node(i)%x(1,3,2)         &
                                      + dn_dpsi_dz * node_list%node(i)%values(1,3,var_psi) * node_list%node(i)%x(1,2,2)         &
                                      + dn_dpsi_dz * node_list%node(i)%values(1,2,var_psi) * node_list%node(i)%x(1,3,2)
    endif

    if (with_neutrals) then

      call neutral_density(xpoint2, xcase2, Z, ES%Z_xpoint, psi,ES%psi_axis,ES%psi_bnd,zrn,drn_dpsi,drn_dz,drn_dpsi2,drn_dz2,             &
                                                                 drn_dpsi_dz,drn_dpsi3,drn_dpsi_dz2, drn_dpsi2_dz)

      node_list%node(i)%values(1,1,var_rhon) = zrn
      node_list%node(i)%values(1,2,var_rhon) = drn_dpsi    * node_list%node(i)%values(1,2,var_psi) + drn_dz * node_list%node(i)%x(1,2,2)
      node_list%node(i)%values(1,3,var_rhon) = drn_dpsi    * node_list%node(i)%values(1,3,var_psi) + drn_dz * node_list%node(i)%x(1,3,2)
      node_list%node(i)%values(1,4,var_rhon) = drn_dpsi    * node_list%node(i)%values(1,4,var_psi) + drn_dz * node_list%node(i)%x(1,4,2) &
                                      + drn_dpsi2   * node_list%node(i)%values(1,2,var_psi) * node_list%node(i)%values(1,3,var_psi)      &
                                      + drn_dz2     * node_list%node(i)%x(1,2,2)            * node_list%node(i)%x(1,3,2)                 &
                                      + drn_dpsi_dz * node_list%node(i)%values(1,3,var_psi) * node_list%node(i)%x(1,2,2)                 &
                                      + drn_dpsi_dz * node_list%node(i)%values(1,2,var_psi) * node_list%node(i)%x(1,3,2)

    endif ! with_neutrals


    if ( with_TiTe ) then
      node_list%node(i)%values(1,1,var_Ti) = zTi
      node_list%node(i)%values(1,2,var_Ti) = dTi_dpsi    * node_list%node(i)%values(1,2,var_psi) + dTi_dz * node_list%node(i)%x(1,2,2)
      node_list%node(i)%values(1,3,var_Ti) = dTi_dpsi    * node_list%node(i)%values(1,3,var_psi) + dTi_dz * node_list%node(i)%x(1,3,2)
      node_list%node(i)%values(1,4,var_Ti) = dTi_dpsi    * node_list%node(i)%values(1,4,var_psi) + dTi_dz * node_list%node(i)%x(1,4,2) &
                                      + dTi_dpsi2   * node_list%node(i)%values(1,2,var_psi) * node_list%node(i)%values(1,3,var_psi) &
                                      + dTi_dz2     * node_list%node(i)%x(1,2,2)              * node_list%node(i)%x(1,3,2)         &
                                      + dTi_dpsi_dz * node_list%node(i)%values(1,3,var_psi) * node_list%node(i)%x(1,2,2)         &
                                      + dTi_dpsi_dz * node_list%node(i)%values(1,2,var_psi) * node_list%node(i)%x(1,3,2)

      node_list%node(i)%values(1,1,var_Te) = zTe
      node_list%node(i)%values(1,2,var_Te) = dTe_dpsi    * node_list%node(i)%values(1,2,var_psi) + dTe_dz * node_list%node(i)%x(1,2,2)
      node_list%node(i)%values(1,3,var_Te) = dTe_dpsi    * node_list%node(i)%values(1,3,var_psi) + dTe_dz * node_list%node(i)%x(1,3,2)
      node_list%node(i)%values(1,4,var_Te) = dTe_dpsi    * node_list%node(i)%values(1,4,var_psi) + dTe_dz * node_list%node(i)%x(1,4,2) &
                                      + dTe_dpsi2   * node_list%node(i)%values(1,2,var_psi) * node_list%node(i)%values(1,3,var_psi) &
                                      + dTe_dz2     * node_list%node(i)%x(1,2,2)              * node_list%node(i)%x(1,3,2)         &
                                      + dTe_dpsi_dz * node_list%node(i)%values(1,3,var_psi) * node_list%node(i)%x(1,2,2)         &
                                      + dTe_dpsi_dz * node_list%node(i)%values(1,2,var_psi) * node_list%node(i)%x(1,3,2)
    else
      node_list%node(i)%values(1,1,var_T) = zT
      node_list%node(i)%values(1,2,var_T) = dT_dpsi    * node_list%node(i)%values(1,2,var_psi) + dT_dz * node_list%node(i)%x(1,2,2)
      node_list%node(i)%values(1,3,var_T) = dT_dpsi    * node_list%node(i)%values(1,3,var_psi) + dT_dz * node_list%node(i)%x(1,3,2)
      node_list%node(i)%values(1,4,var_T) = dT_dpsi    * node_list%node(i)%values(1,4,var_psi) + dT_dz * node_list%node(i)%x(1,4,2) &
                                      + dT_dpsi2   * node_list%node(i)%values(1,2,var_psi) * node_list%node(i)%values(1,3,var_psi) &
                                      + dT_dz2     * node_list%node(i)%x(1,2,2)            * node_list%node(i)%x(1,3,2)         &
                                      + dT_dpsi_dz * node_list%node(i)%values(1,3,var_psi) * node_list%node(i)%x(1,2,2)         &
                                      + dT_dpsi_dz * node_list%node(i)%values(1,2,var_psi) * node_list%node(i)%x(1,3,2)

    end if

    if (tauIC .ne. 0.d0) then
      ! --- u0rho patch: u0 = -tauIC*G(psi), dG/dpsi = (dp/dpsi)/rho; chain rule as for the profiles above (no Z dependence)
      u0_wgt  = (psi - u0_psi_lo) / (u0_psi_hi - u0_psi_lo) * dble(n_u0tab - 1)
      k_u0    = min(max(int(u0_wgt) + 1, 1), n_u0tab - 1)
      u0_wgt  = u0_wgt - dble(k_u0 - 1)
      u0_zG   = (1.d0 - u0_wgt) * u0tab_G(k_u0) + u0_wgt * u0tab_G(k_u0+1)
      u0_zf   = dp_dpsi / zn
      u0_df   = dp_dpsi2 / zn - dp_dpsi * dn_dpsi / zn**2
      node_list%node(i)%values(1,1,var_u) = - tauIC * u0_zG
      node_list%node(i)%values(1,2,var_u) = - tauIC * u0_zf * node_list%node(i)%values(1,2,var_psi)
      node_list%node(i)%values(1,3,var_u) = - tauIC * u0_zf * node_list%node(i)%values(1,3,var_psi)
      node_list%node(i)%values(1,4,var_u) = - tauIC * ( u0_zf * node_list%node(i)%values(1,4,var_psi)                          &
                                          + u0_df * node_list%node(i)%values(1,2,var_psi) * node_list%node(i)%values(1,3,var_psi) )
    else
      ! --- stock initialisation (u0 = -tauIC*p, i.e. 0 for tauIC = 0)
      node_list%node(i)%values(1,1,var_u) = - tauIC * zp 
      node_list%node(i)%values(1,2,var_u) = - tauIC * (dp_dpsi  * node_list%node(i)%values(1,2,var_psi) + dp_dz * node_list%node(i)%x(1,2,2))
      node_list%node(i)%values(1,3,var_u) = - tauIC * (dp_dpsi  * node_list%node(i)%values(1,3,var_psi) + dp_dz * node_list%node(i)%x(1,3,2))
      node_list%node(i)%values(1,4,var_u) = - tauIC * (dP_dpsi  * node_list%node(i)%values(1,4,var_psi) + dP_dz * node_list%node(i)%x(1,4,2) &
                                      + dP_dpsi2   * node_list%node(i)%values(1,2,var_psi) * node_list%node(i)%values(1,3,var_psi)  &
                                      + dP_dz2     * node_list%node(i)%x(1,2,2)        * node_list%node(i)%x(1,3,2)               &
                                      + dP_dpsi_dz * node_list%node(i)%values(1,3,var_psi) * node_list%node(i)%x(1,2,2)         &
                                      + dP_dpsi_dz * node_list%node(i)%values(1,2,var_psi) * node_list%node(i)%x(1,3,2) )
    endif
 
    node_list%node(i)%values(1,:,var_w) = 0.d0        ! vorticity (will be filled just below with inverse Poisson)

    if ( with_vpar ) then
      node_list%node(i)%values(1,:,var_Vpar) = 0.d0        ! parallel velocity
      if ( (abs(V_0) .ge. 1.d-19) .or. (num_rot) ) then
        ! set Vpar,0 (JOREK normalized, ie without unit)= parallel velocity given as input profile 
        if (normalized_velocity_profile) then
          node_list%node(i)%values(1,1,var_Vpar) = zV
          node_list%node(i)%values(1,2,var_Vpar) = dV_dpsi    * node_list%node(i)%values(1,2,1) + dV_dz * node_list%node(i)%x(1,2,2)
          node_list%node(i)%values(1,3,var_Vpar) = dV_dpsi    * node_list%node(i)%values(1,3,1) + dV_dz * node_list%node(i)%x(1,3,2)
          node_list%node(i)%values(1,4,var_Vpar) = dV_dpsi    * node_list%node(i)%values(1,4,1) + dV_dz * node_list%node(i)%x(1,4,2) &
                                                 + dV_dpsi2   * node_list%node(i)%values(1,2,1) * node_list%node(i)%values(1,3,1)    &
                                                 + dV_dz2     * node_list%node(i)%x(1,2,2)      * node_list%node(i)%x(1,3,2)         &
                                                 + dV_dpsi_dz * node_list%node(i)%values(1,3,1) * node_list%node(i)%x(1,2,2)         &
                                                 + dV_dpsi_dz * node_list%node(i)%values(1,2,1) * node_list%node(i)%x(1,3,2) 
        ! set Vpar,0 = 2piR^2/F0 * Omega_tor,0
        ! where Omega_tor,0 is the angular toroidal (~parallel) rotation given as input profile.
        else
          node_list%node(i)%values(1,1,var_Vpar) = R**2 * Omega
          node_list%node(i)%values(1,2,var_Vpar) = 2.d0 * R * Omega       * node_list%node(i)%x(1,2,1)        &
                                                   + R**2   * dOmega_dpsi * node_list%node(i)%values(1,2,1)   &
                                                   + R**2   * dOmega_dz   * node_list%node(i)%x(1,2,2)
          node_list%node(i)%values(1,3,var_Vpar) = 2.d0 * R * Omega       * node_list%node(i)%x(1,3,1)        &
                                                   + R**2   * dOmega_dpsi * node_list%node(i)%values(1,3,1)   &
                                                   + R**2   * dOmega_dz   * node_list%node(i)%x(1,3,2)
        
          node_list%node(i)%values(1,4,var_Vpar) =   2.d0 *     node_list%node(i)%x(1,2,1)**2  * Omega &
                                                   + 2.d0 * R * node_list%node(i)%x(1,4,1)     * Omega &
                                                   + 2.d0 * R * node_list%node(i)%x(1,2,1)     * dOmega_dpsi  * node_list%node(i)%values(1,3,1) &
                                                   + 2.d0 * R * node_list%node(i)%x(1,2,1)     * dOmega_dz    * node_list%node(i)%x(1,3,2) 
          node_list%node(i)%values(1,4,var_Vpar) = node_list%node(i)%values(1,4,var_Vpar) &
                                                   + 2.d0 * R * dOmega_dpsi    * node_list%node(i)%x(1,3,1)      * node_list%node(i)%values(1,2,1) &
                                                   + R**2     * dOmega_dpsi2   * node_list%node(i)%values(1,3,1) * node_list%node(i)%values(1,2,1) &
                                                   + R**2     * dOmega_dpsi_dz * node_list%node(i)%x(1,3,2)      * node_list%node(i)%values(1,2,1) &
                                                   + R**2     * dOmega_dpsi    * node_list%node(i)%values(1,4,1)
          node_list%node(i)%values(1,4,var_Vpar) = node_list%node(i)%values(1,4,var_Vpar) &
                                                   + 2.d0 * R * dOmega_dz      * node_list%node(i)%x(1,3,1)      * node_list%node(i)%x(1,2,2) &
                                                   + R**2     * dOmega_dpsi_dz * node_list%node(i)%values(1,3,1) * node_list%node(i)%x(1,3,2) &
                                                   + R**2     * dOmega_dz2     * node_list%node(i)%x(1,3,2)      * node_list%node(i)%x(1,3,2) &
                                                   + R**2     * dOmega_dz      * node_list%node(i)%x(1,4,2)
        
          node_list%node(i)%values(1,1,var_Vpar) = 2.d0 * PI / F0 * node_list%node(i)%values(1,1,var_Vpar)
          node_list%node(i)%values(1,2,var_Vpar) = 2.d0 * PI / F0 * node_list%node(i)%values(1,2,var_Vpar)
          node_list%node(i)%values(1,3,var_Vpar) = 2.d0 * PI / F0 * node_list%node(i)%values(1,3,var_Vpar)
          node_list%node(i)%values(1,4,var_Vpar) = 2.d0 * PI / F0 * node_list%node(i)%values(1,4,var_Vpar)
        endif
      endif
    end if 
    
    node_list%node(i)%deltas = 0.d0

  enddo

  if (allocated(u0tab_psi)) deallocate(u0tab_psi, u0tab_f, u0tab_G)

endif ! if n_order<=3

! --- Variable projection is better at higher order...
! --- (by the way, we could use this for n_order=3 and remove all the above as well, 
! --- and remove all derivatives from profiles functions, which are not really needed, 
! --- except dn_dpsi and dT_dpsi for current profile...)
if (n_order .ge. 5) then
  if (with_rho) then
    call Poisson(my_id,0,node_list,element_list,bnd_node_list,bnd_elm_list, &
                 var_psi,var_rho,1, ES%psi_axis,ES%psi_bnd,xpoint2,xcase2,ES%Z_xpoint,freeboundary_equil,refinement,1)
  endif
  if ( with_TiTe ) then
    call Poisson(my_id,0,node_list,element_list,bnd_node_list,bnd_elm_list, &
                 var_psi,var_Ti,1, ES%psi_axis,ES%psi_bnd,xpoint2,xcase2,ES%Z_xpoint,freeboundary_equil,refinement,1)
    call Poisson(my_id,0,node_list,element_list,bnd_node_list,bnd_elm_list, &
                 var_psi,var_Te,1, ES%psi_axis,ES%psi_bnd,xpoint2,xcase2,ES%Z_xpoint,freeboundary_equil,refinement,1)
  else
    call Poisson(my_id,0,node_list,element_list,bnd_node_list,bnd_elm_list, &
                 var_psi,var_T,1, ES%psi_axis,ES%psi_bnd,xpoint2,xcase2,ES%Z_xpoint,freeboundary_equil,refinement,1)
  endif
  if ( with_vpar ) then
    call Poisson(my_id,0,node_list,element_list,bnd_node_list,bnd_elm_list, &
                 var_psi,var_Vpar,1, ES%psi_axis,ES%psi_bnd,xpoint2,xcase2,ES%Z_xpoint,freeboundary_equil,refinement,1)
  endif
  if (with_neutrals) then
    call Poisson(my_id,0,node_list,element_list,bnd_node_list,bnd_elm_list, &
                 var_psi,var_rhon,1, ES%psi_axis,ES%psi_bnd,xpoint2,xcase2,ES%Z_xpoint,freeboundary_equil,refinement,1)
  endif
endif



if (tauIC .ne. 0.d0) then
  call Poisson(my_id,2,node_list,element_list,bnd_node_list,bnd_elm_list, &
               var_u,var_w,1, ES%psi_axis,ES%psi_bnd,xpoint2, xcase2,ES%Z_xpoint,freeboundary_equil,refinement,1)      ! inverse Poisson
endif

    
!---------------------------- initialise perturbations
amplitude = 1.d-12
mm = 2

do in=2,n_tor

  if (my_id .eq. 0) then

    do i=1,node_list%n_nodes

      node_list%node(i)%values(in,:,:) = 0.d0

      psi = node_list%node(i)%values(1,1,var_psi)
      Z   = node_list%node(i)%x(1,1,2)
      psi_n = (psi - ES%psi_axis)/(ES%psi_bnd - ES%psi_axis)

      node_list%node(i)%values(in,1,var_w) = amplitude * psi_n * (1.d0 -psi_n)
      node_list%node(i)%values(in,2,var_w) = amplitude * (1.d0 - 2.d0 * psi_n)/(ES%psi_bnd - ES%psi_axis) * node_list%node(i)%values(1,2,var_psi)
      node_list%node(i)%values(in,3,var_w) = amplitude * (1.d0 - 2.d0 * psi_n)/(ES%psi_bnd - ES%psi_axis) * node_list%node(i)%values(1,3,var_psi)
      node_list%node(i)%values(in,4,var_w) = amplitude * (1.d0 - 2.d0 * psi_n)/(ES%psi_bnd - ES%psi_axis) * node_list%node(i)%values(1,4,var_psi)
      
      if (xpoint2 .and. ((psi_n .gt. 1.d0) .or. ((Z .lt. ES%Z_xpoint(1)) .and. (xcase2 .ne. UPPER_XPOINT)) ) ) then
        node_list%node(i)%values(in,1:4,var_w) = 0.d0
      endif
      if (xpoint2 .and. ((psi_n .gt. 1.d0) .or. ((Z .gt. ES%Z_xpoint(2)) .and. (xcase2 .ne. LOWER_XPOINT)) ) ) then
        node_list%node(i)%values(in,1:4,var_w) = 0.d0
      endif

      node_list%node(i)%deltas = 0.d0

    enddo

  endif

  call Poisson(my_id,1,node_list,element_list,bnd_node_list,bnd_elm_list, &
               var_w,var_u,1, ES%psi_axis,ES%psi_bnd,xpoint2, xcase2,ES%Z_xpoint,freeboundary_equil,refinement,1)
enddo

return
end
