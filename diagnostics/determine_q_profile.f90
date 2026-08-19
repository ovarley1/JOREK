!> Calculates the q-profile from the flux surface representation
!! (adapted from helena20)
subroutine determine_q_profile(node_list,element_list,surface_list,psi_axis,psi_xpoint,Z_xpoint,q,rad)

use constants
use tr_module 
use data_structure
use phys_module
use equil_info
use mod_interp


implicit none

! --- Gaussian points between (-1.,1.) for Gauss-integration
real*8, parameter :: xgs(4) = (/-0.861136311594053, -0.339981043584856, 0.339981043584856,  0.861136311594053 /)
real*8, parameter :: wgs(4) = (/ 0.347854845137454,  0.652145154862546, 0.652145154862546,  0.347854845137454 /)

! --- Input parameters.
type (type_node_list),    intent(in)    :: node_list
type (type_element_list), intent(in)    :: element_list
type (type_surface_list), intent(in)    :: surface_list
real*8,                   intent(in)    :: psi_axis
real*8,                   intent(in)    :: psi_xpoint(2)
real*8,                   intent(in)    :: Z_xpoint(2)
real*8,                   intent(inout) :: q(surface_list%n_psi)
real*8,                   intent(inout) :: rad(surface_list%n_psi)

! --- Local variables
integer :: n_int
integer :: i_elm, j, k, n1, n2, n3
real*8  :: t,rr1, rr2, drr1, drr2, ss1, ss2, dss1, dss2, ri, si, dri, dsi, dl
real*8  :: RRgi, dRRgi_dr, dRRgi_ds, ZZgi, dZZgi_dr, dZZgi_ds, dRRgi_dt, dZZgi_dt
real*8  :: PSgi, dPSgi_dr, dPSgi_ds, PSI_R, PSI_Z, RZJAC, grad_psi, psi_n
real*8  :: Fgi,dFgi_dr,dFgi_ds,dFgi_drs,dFgi_drr,dFgi_dss
real*8  :: sum_dl, B_tot2
real*8  :: dRRgi_drs,dRRgi_drr,dRRgi_dss, dZZgi_drs,dZZgi_drr,dZZgi_dss, dPSgi_drs,dPSgi_drr,dPSgi_dss
integer :: i,m, ig, ip

!write(*,*) '   i     psi           q          sum_dl'

  rad(:) = 0.d0
  q(:)   = 0.d0
  Fgi    = 0.d0

do i=2, surface_list%n_psi
  rad(i) = 0.d0
  q(i)   = 0.d0
  sum_dl = 0.d0
  do k=1, surface_list%flux_surfaces(i)%n_pieces
    do ig = 1, 4
      t = xgs(ig)

      rr1  = surface_list%flux_surfaces(i)%s(1,k)
      drr1 = surface_list%flux_surfaces(i)%s(2,k)
      rr2  = surface_list%flux_surfaces(i)%s(3,k)
      drr2 = surface_list%flux_surfaces(i)%s(4,k)

      ss1  = surface_list%flux_surfaces(i)%t(1,k)
      dss1 = surface_list%flux_surfaces(i)%t(2,k)
      ss2  = surface_list%flux_surfaces(i)%t(3,k)
      dss2 = surface_list%flux_surfaces(i)%t(4,k)

      call CUB1D(rr1, drr1, rr2, drr2, t, ri, dri)
      call CUB1D(ss1, dss1, ss2, dss2, t, si, dsi)

      i_elm = surface_list%flux_surfaces(i)%elm(k)

      call interp(node_list,element_list,i_elm,1,1,ri,si,PSgi,dPSgi_dr,dPSgi_ds,dPSgi_drs,dPSgi_drr,dPSgi_dss)

      call interp_RZ(node_list,element_list,i_elm,ri,si,RRgi,dRRgi_dr,dRRgi_ds,dRRgi_drs,dRRgi_drr,dRRgi_dss, &
                                                        ZZgi,dZZgi_dr,dZZgi_ds,dZZgi_drs,dZZgi_drr,dZZgi_dss)
                                                        
      ! --- Ignore open and private field line regions
      if ( get_psi_n(PSgi, ZZgi) > 1.d0 ) cycle

      dRRgi_dt = dRRgi_dr * dri + dRRgi_ds * dsi
      dZZgi_dt = dZZgi_dr * dri + dZZgi_ds * dsi

      dl = sqrt(dRRgi_dt**2 + dZZgi_dt**2)

      RZjac  = DRRgi_dr * dZZgi_ds - dRRgi_ds * dZZgi_dr

      PSI_R = (   dPSgi_dr * dZZgi_ds - dPSgi_ds * dZZgi_dr ) / RZjac
      PSI_Z = ( - dPSgi_dr * dRRgi_ds + dPSgi_ds * dRRgi_dr ) / RZjac

      grad_psi = sqrt(PSI_R * PSI_R + PSI_Z * PSI_Z)
      grad_psi = sign(1d0, ES%psi_axis - ES%psi_bnd ) * grad_psi ! to get q profile with the correct sign 

#ifdef fullmhd
      call interp(node_list,element_list,i_elm, 710 ,1,ri,si, Fgi, dFgi_dr,dFgi_ds,dFgi_drs,dFgi_drr,dFgi_dss)  ! ivar = 710 for Fprof_eq
#else
      Fgi = F0
#endif

      B_tot2 =  (Fgi / RRgi)**2 + (grad_psi/ RRgi)**2

      sum_dl = sum_dl +  wgs(ig) * dl

      q(i) = q(i) +  wgs(ig) / (RRgi * grad_psi) * dl
      rad(i)= rad(i) + sqrt( (RRgi-R_geo)**2.+(ZZgi-Z_geo)**2.)
    end do
  end do

  q(i) = Fgi * q(i) / (2.d0 * PI)
  if ( surface_list%flux_surfaces(i)%n_pieces /= 0 ) then
    rad(i)=rad(i)/(4.d0*surface_list%flux_surfaces(i)%n_pieces)
  end if
!  write(*,'(i5,3es13.5)') i, surface_list%psi_values(i), q(i), sum_dl
end do

!----------------------------------- values on axis
!q(1) =  PI / sqrt(CRR_axis*CZZ_Axis)

end subroutine determine_q_profile

subroutine determine_Itor_profile(node_list,element_list,surface_list,psi_axis,psi_xpoint,Z_xpoint,Itor,Area)

use constants
use mod_parameters
use tr_module 
use data_structure
use gauss
use basis_at_gaussian
use phys_module
use equil_info

implicit none

! --- Input parameters.
type (type_node_list),    intent(in)    :: node_list
type (type_element_list), intent(in)    :: element_list
type (type_surface_list), intent(in)    :: surface_list
real*8,                   intent(in)    :: psi_axis
real*8,                   intent(in)    :: psi_xpoint(2)
real*8,                   intent(in)    :: Z_xpoint(2)
real*8,                   intent(inout) :: Itor(surface_list%n_psi)
real*8,                   intent(inout) :: Area(surface_list%n_psi)

! --- Local variables
type (type_element)      :: element
type (type_node)         :: nodes(n_vertex_max)
real*8     :: x_g(n_gauss,n_gauss), x_s(n_gauss,n_gauss), x_t(n_gauss,n_gauss)
real*8     :: y_g(n_gauss,n_gauss), y_s(n_gauss,n_gauss), y_t(n_gauss,n_gauss)
real*8     :: eq_g(0:n_var,n_gauss,n_gauss)
integer    :: i, j, k, in, ms, mt, iv, inode, ife, n_elements
real*8     :: xjac, BigR, wst, ZJ_0, PS_0, psin, psin_tmp

Itor(:)     = 0.d0
Area(:)     = 0.d0

do ife =1, element_list%n_elements

  element = element_list%element(ife)

  do iv = 1, n_vertex_max
    inode     = element%vertex(iv)
    nodes(iv) = node_list%node(inode)
  enddo

  x_g(:,:)    = 0.d0; x_s(:,:)    = 0.d0; x_t(:,:)    = 0.d0;
  y_g(:,:)    = 0.d0; y_s(:,:)    = 0.d0; y_t(:,:)    = 0.d0;
  eq_g(:,:,:) = 0.d0

  do i=1,n_vertex_max
    do j=1,n_degrees
      do ms=1, n_gauss
        do mt=1, n_gauss

          do in=1, n_coord_tor
            x_g(ms,mt) = x_g(ms,mt) + nodes(i)%x(in,j,1) * element%size(i,j) * H(i,j,ms,mt) * HZ_coord(in,i_plane_rtree)
            y_g(ms,mt) = y_g(ms,mt) + nodes(i)%x(in,j,2) * element%size(i,j) * H(i,j,ms,mt) * HZ_coord(in,i_plane_rtree)
            
            x_s(ms,mt) = x_s(ms,mt) + nodes(i)%x(in,j,1) * element%size(i,j) * H_s(i,j,ms,mt) * HZ_coord(in,i_plane_rtree)
            x_t(ms,mt) = x_t(ms,mt) + nodes(i)%x(in,j,1) * element%size(i,j) * H_t(i,j,ms,mt) * HZ_coord(in,i_plane_rtree)
            y_s(ms,mt) = y_s(ms,mt) + nodes(i)%x(in,j,2) * element%size(i,j) * H_s(i,j,ms,mt) * HZ_coord(in,i_plane_rtree)
            y_t(ms,mt) = y_t(ms,mt) + nodes(i)%x(in,j,2) * element%size(i,j) * H_t(i,j,ms,mt) * HZ_coord(in,i_plane_rtree)
          enddo
        enddo
      enddo
    enddo
  enddo

  eq_g(:,:,:) = 0.d0

  do i=1,n_vertex_max
    do j=1,n_degrees
      do ms=1, n_gauss
        do mt=1, n_gauss

          do k=1,n_var
            eq_g(k,ms,mt)  = eq_g(k,ms,mt)  + nodes(i)%values(1,j,k) * element%size(i,j) * H(i,j,ms,mt)
          enddo

        enddo
      enddo
    enddo
  enddo
!--------------------------------------------------- sum over the Gaussian integration points

  do ms=1, n_gauss

    do mt=1, n_gauss

      wst = wgauss(ms)*wgauss(mt)

      xjac = x_s(ms,mt)*y_t(ms,mt) - x_t(ms,mt)*y_s(ms,mt)
      BigR = x_g(ms,mt)

      !! rho_00 = eq_g(var_rho,ms,mt)
      !! if (with_TiTe) then
      !!   Ti_00 = eq_g(var_Ti,ms,mt)
      !!   Te_00 = eq_g(var_Te,ms,mt)
      !!   T_00  = Ti_00 + Te_00
      !! else
      !!   T_00  = eq_g(var_T,ms,mt)
      !! endif
      ZJ_0  = eq_g(var_zj,ms,mt)
      PS_0  = eq_g(var_psi,ms,mt) 

      psin = (PS_0-psi_axis)/(psi_xpoint(1)-psi_axis)

      ! --- Ignore open and private field line regions
      if ( get_psi_n(PS_0, y_g(ms,mt)) > 1.d0 ) cycle

      
      do k=1, surface_list%n_psi
        psin_tmp = (surface_list%psi_values(k)-psi_axis)/(psi_xpoint(1)-psi_axis)

        if ( psin .le. psin_tmp ) then
          Itor(k)     = Itor(k)     - ZJ_0 /BigR    * xjac * wst
          Area(k)     = Area(k)     +                 xjac * wst
        endif
      enddo
      
    enddo
  enddo
enddo

!----------------------------------- values on axis
!q(1) =  PI / sqrt(CRR_axis*CZZ_Axis)

end subroutine determine_Itor_profile


subroutine determine_heatflux_profile(node_list,element_list,surface_list,psi_axis,psi_xpoint,Z_xpoint,Qheat,Area)

use constants
use mod_parameters
use tr_module 
use data_structure
use gauss
use basis_at_gaussian
use phys_module
use equil_info
use diffusivities

implicit none

! --- Input parameters.
type (type_node_list),    intent(in)    :: node_list
type (type_element_list), intent(in)    :: element_list
type (type_surface_list), intent(in)    :: surface_list
real*8,                   intent(in)    :: psi_axis
real*8,                   intent(in)    :: psi_xpoint(2)
real*8,                   intent(in)    :: Z_xpoint(2)
real*8,                   intent(inout) :: Qheat(surface_list%n_psi)
real*8,                   intent(inout) :: Area(surface_list%n_psi)

! --- Local variables
type (type_element)      :: element
type (type_node)         :: nodes(n_vertex_max), node_k
real*8     :: x_g(n_gauss,n_gauss), x_s(n_gauss,n_gauss), x_t(n_gauss,n_gauss)
real*8     :: y_g(n_gauss,n_gauss), y_s(n_gauss,n_gauss), y_t(n_gauss,n_gauss)
real*8     :: eq_g(0:n_var,n_gauss,n_gauss), eq_s(0:n_var,n_gauss,n_gauss)
real*8     :: eq_t(0:n_var,n_gauss,n_gauss)
integer    :: i, j, k, in, ms, mt, iv, inode, ife, n_elements
real*8     :: xjac, BigR, wst, ZJ_0, PS_0, rho0, T0, psin, psin_tmp
real*8     :: T0_s, T0_t, ps0_s, ps0_t
real*8     :: dTdx, dTdy, dpsidx, dpsidy, dpsidp, grad_psi
real*8     :: ZK_prof, cond_perp
real*8     :: fact_rho, fact_resistiv, fact_T, fact_flux
!real*8     :: R_c, Z_c, vec_inside(2), grad_t(2)

Qheat(:)    = 0.d0
Area(:)     = 0.d0

cond_perp   = 0.

do ife =1, element_list%n_elements

  element = element_list%element(ife)

  do iv = 1, n_vertex_max
    inode     = element%vertex(iv)
    nodes(iv) = node_list%node(inode)
  enddo

  x_g(:,:)    = 0.d0; x_s(:,:)    = 0.d0; x_t(:,:)    = 0.d0;
  y_g(:,:)    = 0.d0; y_s(:,:)    = 0.d0; y_t(:,:)    = 0.d0;
  eq_g(:,:,:) = 0.d0

  cond_perp   = 0.

  do i=1,n_vertex_max
    do j=1,n_degrees
      do ms=1, n_gauss
        do mt=1, n_gauss

          ! 1D only
          x_g(ms,mt) = x_g(ms,mt) + nodes(i)%x(1,j,1) * element%size(i,j) * H(i,j,ms,mt)
          y_g(ms,mt) = y_g(ms,mt) + nodes(i)%x(1,j,2) * element%size(i,j) * H(i,j,ms,mt)
          
          x_s(ms,mt) = x_s(ms,mt) + nodes(i)%x(1,j,1) * element%size(i,j) * H_s(i,j,ms,mt)
          x_t(ms,mt) = x_t(ms,mt) + nodes(i)%x(1,j,1) * element%size(i,j) * H_t(i,j,ms,mt)
          y_s(ms,mt) = y_s(ms,mt) + nodes(i)%x(1,j,2) * element%size(i,j) * H_s(i,j,ms,mt)
          y_t(ms,mt) = y_t(ms,mt) + nodes(i)%x(1,j,2) * element%size(i,j) * H_t(i,j,ms,mt)

          do k=1,n_var
            eq_g(k,ms,mt) = eq_g(k,ms,mt)  + nodes(i)%values(1,j,k) * element%size(i,j) * H(i,j,ms,mt)
          enddo

        enddo
      enddo
    enddo
  enddo

!--------------------------------------------------- sum over the Gaussian integration points

  do ms=1, n_gauss
    do mt=1, n_gauss

      wst = wgauss(ms)*wgauss(mt)

      xjac = x_s(ms,mt)*y_t(ms,mt) - x_t(ms,mt)*y_s(ms,mt)
      BigR = x_g(ms,mt)

      PS_0  = eq_g(var_psi,ms,mt) 
      ps0_s = eq_s(var_psi,ms,mt) 
      ps0_t = eq_t(var_psi,ms,mt)
      psin = (PS_0-psi_axis)/(psi_xpoint(1)-psi_axis)

      ZJ_0  = eq_g(var_zj,ms,mt)
      rho0 = eq_g(var_rho,ms,mt)
      T0   = eq_g(var_T,ms,mt)
      T0_s = eq_s(var_T,ms,mt) 
      T0_t = eq_t(var_T,ms,mt)

      dpsidx   = (   y_t(ms,mt) * ps0_s - y_s(ms,mt) * ps0_t ) / xjac
      dpsidy   = ( - x_t(ms,mt) * ps0_s + x_s(ms,mt) * ps0_t ) / xjac
      grad_psi = sqrt(dpsidx*dpsidx + dpsidy*dpsidy)

      dTdx   = (   y_t(ms,mt) * T0_s - y_s(ms,mt) * T0_t ) / xjac
      dTdy   = ( - x_t(ms,mt) * T0_s + x_s(ms,mt) * T0_t ) / xjac

      fact_rho      = central_density * 1.d20 * central_mass*MASS_PROTON
      fact_resistiv = sqrt ( MU_zero / fact_rho )
      fact_T        = 1.d0 / ( MU_zero * central_density * 1.d20 * EL_CHG )
      fact_flux     = 1. / ( MU_zero * sqrt ( MU_zero*fact_rho ) )  

      if (set_chi_perp_constant) then
        ZK_prof = get_zkperp(psin, y_s(ms,mt)) * max(rho0,1.d-2) / (gamma - 1) / fact_resistiv / fact_rho
        if (T0 .lt. ZK_prof_neg_thresh) then
          ZK_prof = ZK_prof_neg * max(rho0,1.d-2)                / (gamma - 1) / fact_resistiv / fact_rho
        endif
      else
        ZK_prof = get_zkperp(psin, y_s(ms,mt))                   / (gamma - 1) / fact_resistiv / fact_rho / rho0
        if (T0 .lt. ZK_prof_neg_thresh) then
          ZK_prof = ZK_prof_neg                                  / (gamma - 1) / fact_resistiv / fact_rho / rho0
        endif
      end if

      if (grad_psi .gt. 1.d-6) then
        cond_perp   = - ZK_prof * (dTdx*dpsidx + dTdy*dpsidy) / grad_psi / fact_T
      endif

      ! --- Ignore open and private field line regions
      if ( get_psi_n(PS_0, y_g(ms,mt)) > 1.d0 ) cycle

      
      do k=1, surface_list%n_psi
        psin_tmp = (surface_list%psi_values(k)-psi_axis)/(psi_xpoint(1)-psi_axis)

        if ( psin .le. psin_tmp ) then
          Qheat(k)    = Qheat(k)    + cond_perp * xjac * wst
          Area(k)     = Area(k)     +             xjac * wst
        endif
      enddo
      
    enddo
  enddo
enddo

!----------------------------------- values on axis
!q(1) =  PI / sqrt(CRR_axis*CZZ_Axis)

end subroutine determine_heatflux_profile
