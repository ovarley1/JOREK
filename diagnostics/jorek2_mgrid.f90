!> Create a file with pressure, magnetic field, plasma beta to compare with HINT.
program jorek2_fields_orin
  
use data_structure
use phys_module
use basis_at_gaussian
use mod_parameters
use mod_log_params
use mod_import_restart
use mod_interp
use mod_chi
use mpi_mod
implicit none

type(type_node_list)    :: node_list
type(type_element_list) :: element_list
integer :: my_id, ierr, ifail, checked_elms
integer :: nR, nZ, nfp, n_planes_out
real*8 :: R_min, R_max, Z_min, Z_max
real*8, allocatable :: pressure(:,:,:), magnetic_pressure(:,:,:)
real*8, allocatable :: BR(:,:,:), BZ(:,:,:), Bphi(:,:,:), B_vac(:,:,:), B_pert(:,:,:)  ! (nR, nZ, n_planes_out)
integer :: i_R, i_Z, i_plane, i_elm
real*8 :: R, Z, phi, R_out, R_s_out, R_t_out, R_p_out, Z_out, Z_s_out, Z_t_out, Z_p_out, s, t
real*8 :: rho, temp, P, P_s, P_t, P_st, P_ss, P_tt, dummy
real*8 :: psi, psi_s, psi_t, psi_p, dpsidx, dpsidy, dpsidp
real*8 :: xjac, Bv2
real*8 :: dR, dZ, dphi, fact_mu0
character(len=256) :: dummy_s
logical :: in_domain
real*8 :: max_pres = 0.d0
real*8, dimension(0:n_order-1,0:n_order-1,0:n_order-1) :: chi

! required = MPI_THREAD_FUNNELED
! call MPI_Init_thread(required, provided, StatInfo)
! call init_threads()  ! on some systems init_threads needs to come after mpi_init_thread
! call MPI_COMM_SIZE(MPI_COMM_WORLD, comm_size, ierr)
! n_cpu = comm_size

! ! --- Determine ID of each MPI proc
! call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
! my_id = rank
my_id = 0

if ( my_id == 0 ) then
  write(*,*) '***************************************'
  write(*,*) '* JOREK2_fields_orin                  *'
  write(*,*) '***************************************'
endif

! Initialize MPI
! call MPI_Init(ierr)
! call MPI_Comm_rank(MPI_COMM_WORLD, my_id, ierr)

! --- Initialize mode and mode_type arrays
call det_modes()
call initialise_basis
call init_chi_basis

call initialise_parameters(my_id,  "__NO_FILENAME__")
call log_parameters(my_id)

call import_restart(node_list,element_list, 'jorek_restart', rst_format, ierr, .true.)


! ===== READ GRID CONFIGURATION =====
! Read configuration file for grid parameters
open(21, file='grid_config.dat', status='old', action='read', iostat=ierr)
if ( ierr == 0 ) then ! config file exists, use it.
  
  read(21, '(a)') dummy_s ! read comment line (ignored)
  read(21,*) R_min, R_max, nR
  read(21,*) Z_min, Z_max, nZ
  read(21,*) n_planes_out, nfp
  close(21)
  
  if (my_id == 0) then
    write(*,*) 'Grid configuration read from grid_config.dat:'
    write(*,*) '  R: ', R_min, ' to ', R_max, ' with ', nR, ' points'
    write(*,*) '  Z: ', Z_min, ' to ', Z_max, ' with ', nZ, ' points'
    write(*,*) '  Planes: ', n_planes_out, ', Field periods: ', nfp
  endif

else ! if no config file exists, use defaults
  R_min   = 1.0d0
  R_max   = 3.0d0
  nR      = 121
  Z_min   = -1.5d0
  Z_max   =  1.5d0
  nZ      = 121
  n_planes_out = 128
  nfp     = 4
  
  if (my_id == 0) then
    write(*,*) 'Warning: grid_config.dat not found, using default grid:'
    write(*,*) '  R: ', R_min, ' to ', R_max, ' with ', nR, ' points'
    write(*,*) '  Z: ', Z_min, ' to ', Z_max, ' with ', nZ, ' points'
    write(*,*) '  Planes: ', n_planes_out, ', Field periods: ', nfp
  endif
end if

fact_mu0 = 1/MU_ZERO

! ===== CREATE GRID ARRAYS =====
allocate(pressure(nR, nZ, n_planes_out))
allocate(magnetic_pressure(nR, nZ, n_planes_out))
allocate(BR(nR, nZ, n_planes_out))
allocate(BZ(nR, nZ, n_planes_out))
allocate(Bphi(nR, nZ, n_planes_out))
allocate(B_vac(nR, nZ, n_planes_out))
allocate(B_pert(nR, nZ, n_planes_out))

! Initialize variables to NaN
pressure = 0.d0
pressure = pressure / pressure
magnetic_pressure = 0.d0
magnetic_pressure = magnetic_pressure / magnetic_pressure
BR = 0.d0
BR = BR / BR
BZ = 0.d0
BZ = BZ / BZ
Bphi = 0.d0
Bphi = Bphi / Bphi
B_vac = 0.d0
B_vac = B_vac / B_vac
B_pert = 0.d0
B_pert = B_pert / B_pert

if (my_id == 0) then
  write(*,*) ''
  write(*,*) '===== Computing pressure on regular grid ====='
  write(*,*) 'Total grid points: ', nR * nZ * n_planes_out
endif

! ===== LOOP THROUGH GRID AND COMPUTE physical values =====
!$omp parallel default(none) &
!$omp   shared(node_list, element_list, nR, nZ, n_planes_out, nfp, R_min, R_max, Z_min, Z_max, pressure,      &
!$omp          magnetic_pressure, BR, BZ, Bphi, B_vac, B_pert, fact_mu0, F0)                                                 &
!$omp   private(i_R, i_Z, i_plane, R, Z, phi, R_out, R_s_out, R_t_out, R_p_out, Z_out, Z_s_out, Z_t_out,      &
!$omp           Z_p_out, i_elm, s, t, ifail, rho, temp, dummy, checked_elms, chi, Bv2,                        &
!$omp           psi, psi_s, psi_t, psi_p, dpsidx, dpsidy, dpsidp, xjac)
!$omp do
do i_plane = 1, n_planes_out
  phi = 2.d0 * PI / (nfp) * (i_plane - 1) / (n_planes_out)
  do i_Z = 1, nZ
    Z = Z_min + (Z_max - Z_min) * (i_Z - 1) / (nZ - 1)
    do i_R = 1, nR
      R = R_min + (R_max - R_min) * (i_R - 1) / (nR - 1)
      
      ! Check if point is in domain
      call find_RZP(node_list,element_list,R,Z,phi,R_out,Z_out,i_elm,s,t,ifail,checked_elms)

      if (ifail == 0) then ! Point is in domain - compute values
        ! Check the output R,Z is close to input R,Z
        if (R - R_out > 1e-6 .or. Z - Z_out > 1e-6) then
          write(*,*) 'Error in find_RZP: point outside element despite ifail=0'
          write(*,*) '  R,Z,phi = ', R, Z, phi
          write(*,*) '  R_out,Z_out = ', R_out, Z_out
          stop
        end if

        ! Get density and temperature then calculate pressure at this location
        call var_value(node_list,element_list,i_elm,var_rho,s,t,phi,rho, dummy, dummy, dummy)
        call var_value(node_list,element_list,i_elm,var_T,s,t,phi,temp, dummy, dummy, dummy)

        pressure(i_R, i_Z, i_plane) = rho * temp * fact_mu0

        ! Calculate magnetic pressure.
        chi = get_chi(R,Z,phi,node_list,element_list,i_elm,s,t)   ! Vacuum field
        Bv2 = chi(1,0,0)**2 + chi(0,1,0)**2 + chi(0,0,1)**2/R**2  ! Vacuum B squared

        ! Calculate gradient of psi: dpsi/dx, dpsi/dy, dpsi/dphi
        ! First get psi and its derivatives at the location (s,t,phi)
        call var_value(node_list,element_list,i_elm,var_psi,s,t,phi,psi,psi_s,psi_t,psi_p)
        ! Now get the derivatives of R,Z with respect to s,t,phi
        call interp_RZP(node_list,element_list,i_elm,s,t,phi,R_out,R_s_out,R_t_out,R_p_out,Z_out,Z_s_out,Z_t_out,Z_p_out)
        xjac = R_s_out * Z_t_out - R_t_out * Z_s_out

        dpsidx = (Z_t_out * psi_s - Z_s_out * psi_t) / xjac
        dpsidy = (-Z_t_out * psi_s + R_s_out * psi_t) / xjac
        dpsidp = psi_p - dpsidx * R_p_out - dpsidy * Z_p_out
  
        magnetic_pressure(i_R, i_Z, i_plane) = (Bv2*(1.d0 + (dpsidx**2 + dpsidy**2 + dpsidp**2/R**2)/F0**2) - &
                      ((chi(1,0,0)*dpsidx + chi(0,1,0)*dpsidy + chi(0,0,1)*dpsidp/R**2)**2)/F0**2) * fact_mu0 / 2.d0
        if (magnetic_pressure(i_R, i_Z, i_plane) .ne. magnetic_pressure(i_R, i_Z, i_plane)) then
          write(*,*) 'Error: magnetic pressure is NaN at R,Z,phi = ', R, Z, phi
          stop
        end if

        ! Calculate magnetic field components
        BR(i_R, i_Z, i_plane)   = chi(1,0,0)   + (dpsidy * chi(0,0,1) - dpsidp * chi(0,1,0))/(R*F0)
        BZ(i_R, i_Z, i_plane)   = chi(0,1,0)   - (dpsidx * chi(0,0,1) - dpsidp * chi(1,0,0))/(R*F0)
        Bphi(i_R, i_Z, i_plane) = chi(0,0,1)/R + (dpsidx * chi(0,1,0) - dpsidy * chi(1,0,0))/F0

        ! B_vac and B_pert components
        B_vac(i_R, i_Z, i_plane)  = sqrt(Bv2)
        B_pert(i_R, i_Z, i_plane) = sqrt((BR(i_R, i_Z, i_plane)   - chi(1,0,0))**2 + &
                                         (BZ(i_R, i_Z, i_plane)   - chi(0,1,0))**2 + &
                                         (Bphi(i_R, i_Z, i_plane) - chi(0,0,1)/R)**2)

      else if (R .ge. 0) then
        BR(i_R, i_Z, i_plane) = 1.d3
        BZ(i_R, i_Z, i_plane) = 0.d0
        Bphi(i_R, i_Z, i_plane) = 0.d0
      else
        BR(i_R, i_Z, i_plane) = 1.d3
        BZ(i_R, i_Z, i_plane) = 0.d0
        Bphi(i_R, i_Z, i_plane) = 0.d0
      endif ! ifail == 0
    enddo ! i_R
  enddo ! i_Z
enddo ! i_plane
!$omp end do
!$omp end parallel

if (my_id == 0) then
  write(*,*) ''
  write(*,*) '===== Grid calculation complete ====='
endif

! ===== SAVE TO NETCDF FILE =====
call save_to_netcdf(R_min, R_max, nR, Z_min, Z_max, nZ, n_planes_out, nfp, pressure, magnetic_pressure, BR, BZ, Bphi, B_vac, B_pert)
! call save_to_ascii(R_min, R_max, nR, Z_min, Z_max, nZ, n_planes_out, nfp, pressure)


! ===== CLEANUP =====
deallocate(pressure)

if (my_id == 0) then
  write(*,*) ''
  write(*,*) '===== Diagnostic complete ====='
endif





! Finalize
! call MPI_Finalize(ierr)



contains 

!> Get variable and derivatives at a specific location
subroutine var_value(node_list, element_list, i_elm, i_var, s_in, t_in, phi_in, &
                      value_out, value_s_out, value_t_out, value_p_out)
  type (type_node_list),    intent(in)  :: node_list
  type (type_element_list), intent(in)  :: element_list
  integer,                  intent(in)  :: i_elm, i_var
  real*8,                   intent(in)  :: s_in, t_in, phi_in
  real*8,                   intent(out) :: value_out, value_s_out, value_t_out, value_p_out
  
  ! Local variables
  integer :: i_tor, i_harm
  real*8 :: V0, V0_s, V0_t, V0_st, V0_ss, V0_tt, V0_p
  real*8 :: Vcos, Vcos_s, Vcos_t, Vcos_st, Vcos_ss, Vcos_tt, Vcos_p
  real*8 :: Vsin, Vsin_s, Vsin_t, Vsin_st, Vsin_ss, Vsin_tt, Vsin_p
  
  ! Get n=0 mode
  call interp(node_list, element_list, i_elm, i_var, 1, s_in, t_in, &
              V0, V0_s, V0_t, V0_st, V0_ss, V0_tt)
  value_out   = V0
  value_s_out = V0_s
  value_t_out = V0_t
  value_p_out = 0.d0
  
  ! Add Fourier modes
  do i_tor = 1, (n_tor-1)/2
    i_harm = 2*i_tor
    
    ! Cosine mode
    call interp(node_list, element_list, i_elm, i_var, i_harm, s_in, t_in, &
                Vcos, Vcos_s, Vcos_t, Vcos_st, Vcos_ss, Vcos_tt)
    value_out   = value_out   + Vcos   * cos(mode(i_harm)*phi_in)
    value_s_out = value_s_out + Vcos_s * cos(mode(i_harm)*phi_in)
    value_t_out = value_t_out + Vcos_t * cos(mode(i_harm)*phi_in)
    value_p_out = value_p_out - Vcos   * sin(mode(i_harm)*phi_in) * mode(i_harm)
    
    ! Sine mode
    call interp(node_list, element_list, i_elm, i_var, i_harm+1, s_in, t_in, &
                Vsin, Vsin_s, Vsin_t, Vsin_st, Vsin_ss, Vsin_tt)
    value_out   = value_out   + Vsin   * sin(mode(i_harm+1)*phi_in)
    value_s_out = value_s_out + Vsin_s * sin(mode(i_harm+1)*phi_in)
    value_t_out = value_t_out + Vsin_t * sin(mode(i_harm+1)*phi_in)
    value_p_out = value_p_out + Vsin   * cos(mode(i_harm+1)*phi_in) * mode(i_harm+1)
  enddo
  
end subroutine var_value



!> Fallback: Save to ASCII format if NetCDF not available
subroutine save_to_ascii(R_min, R_max, nR, Z_min, Z_max, nZ, n_planes_out, nfp, pressure, magnetic_pressure)
  
  implicit none
  
  integer, intent(in) :: nR, nZ, n_planes_out, nfp
  real*8, intent(in)  :: R_min, R_max, Z_min, Z_max
  real*8, intent(in)  :: pressure(nR, nZ, n_planes_out)
  real*8, intent(in)  :: magnetic_pressure(nR, nZ, n_planes_out)
  integer :: i_R, i_Z, i_plane
  
  if (my_id /= 0) return
  
  open(30, file='jorek_data.dat', status='replace')
  
  ! Write header
  write(30,*) '# JOREK values on regular grid'
  write(30,*) '# Columns: R, Z, phi, pressure, magnetic_pressure'
  write(30,*) '# Radial grid: R_min=', R_min, ' R_max=', R_max, ' nR=', nR
  write(30,*) '# Vertical grid: Z_min=', Z_min, ' Z_max=', Z_max, ' nZ=', nZ
  write(30,*) '# Toroidal grid: n_planes_out=', n_planes_out, ' field periods=', nfp
  write(30,*) '#'
  
  ! Write data
  do i_plane = 1, n_planes_out
    phi = 2.d0 * PI / nfp * (i_plane - 1) / (n_planes_out - 1)
    do i_Z = 1, nZ
      Z = Z_min + (Z_max - Z_min) * (i_Z - 1) / (nZ - 1)
      do i_R = 1, nR
        R = R_min + (R_max - R_min) * (i_R - 1) / (nR - 1)

        write(30,'(4ES16.8)') R, Z, phi, pressure(i_R, i_Z, i_plane), magnetic_pressure(i_R, i_Z, i_plane)
      enddo
    enddo
  enddo


  close(30)

  write(*,*) 'ASCII file written: jorek_data.dat'
  
end subroutine save_to_ascii



!> Save pressure data to NetCDF file in mgrid-compatible format
subroutine save_to_netcdf(R_min, R_max, nR, Z_min, Z_max, nZ, n_planes_out, nfp, pressure, magnetic_pressure, BR, BZ, Bphi, B_vac, B_pert)
  use netcdf
  implicit none  
  ! Arguments
  integer, intent(in) :: nR, nZ, n_planes_out, nfp
  real*8, intent(in)  :: R_min, R_max, Z_min, Z_max
  real*8, intent(in)  :: pressure(nR, nZ, n_planes_out), magnetic_pressure(nR, nZ, n_planes_out)
  real*8, intent(in)  :: BR(nR, nZ, n_planes_out), BZ(nR, nZ, n_planes_out), Bphi(nR, nZ, n_planes_out)
  real*8, intent(in)  :: B_vac(nR, nZ, n_planes_out), B_pert(nR, nZ, n_planes_out)
  
  ! NetCDF IDs
  integer :: ncid
  integer :: dim_rad_id, dim_zee_id, dim_phi_id
  integer :: var_pres_id, var_mag_pres_id, var_br_id, var_bz_id, var_bp_id, var_bvac_id, var_bpert_id
  integer :: var_rmax_id, var_rmin_id, var_zmax_id, var_zmin_id
  integer :: var_ir_id, var_jz_id, var_kp_id, var_nfp_id
  character(len=256) :: filename
  
  ! Only rank 0 writes the file
  if (my_id /= 0) return
  
  filename = 'jorek_data.nc'
  
  write(*,*) ''
  write(*,*) '===== Writing NetCDF file: ', trim(filename), ' ====='

! #ifdef USE_NETCDF
  
  ! Create the NetCDF file
  call check(nf90_create(trim(filename), NF90_CLOBBER, ncid))

  ! ----------------------------------
  ! Define dimensions
  ! ----------------------------------
  call check(nf90_def_dim(ncid, "rad", nR, dim_rad_id))
  call check(nf90_def_dim(ncid, "zee", nZ, dim_zee_id))
  call check(nf90_def_dim(ncid, "phi", n_planes_out, dim_phi_id))
  
  ! ----------------------------------
  ! Define variables 
  ! ----------------------------------
  ! 3D Field arrays (dimension order: R, Z, phi)
  call check(nf90_def_var(ncid, "pres",   NF90_DOUBLE, (/ dim_rad_id, dim_zee_id, dim_phi_id /), var_pres_id))
  call check(nf90_def_var(ncid, "mag_pres",   NF90_DOUBLE, (/ dim_rad_id, dim_zee_id, dim_phi_id /), var_mag_pres_id))
  call check(nf90_def_var(ncid, "br_001",   NF90_DOUBLE, (/ dim_rad_id, dim_zee_id, dim_phi_id /), var_br_id))
  call check(nf90_def_var(ncid, "bz_001",   NF90_DOUBLE, (/ dim_rad_id, dim_zee_id, dim_phi_id /), var_bz_id))
  call check(nf90_def_var(ncid, "bp_001",   NF90_DOUBLE, (/ dim_rad_id, dim_zee_id, dim_phi_id /), var_bp_id))
  call check(nf90_def_var(ncid, "bvac_001",   NF90_DOUBLE, (/ dim_rad_id, dim_zee_id, dim_phi_id /), var_bvac_id))
  call check(nf90_def_var(ncid, "bpert_001",   NF90_DOUBLE, (/ dim_rad_id, dim_zee_id, dim_phi_id /), var_bpert_id))
  ! Note this uses Fortran ordering so it corresponds to a python/C array pres[R,Z,phi]

  ! Geometry Scalars (Doubles)
  call check(nf90_def_var(ncid, "rmax", NF90_DOUBLE, var_rmax_id))
  call check(nf90_def_var(ncid, "rmin", NF90_DOUBLE, var_rmin_id))
  call check(nf90_def_var(ncid, "zmax", NF90_DOUBLE, var_zmax_id))
  call check(nf90_def_var(ncid, "zmin", NF90_DOUBLE, var_zmin_id))

  ! Integer Parameters
  call check(nf90_def_var(ncid, "ir", NF90_INT, var_ir_id))
  call check(nf90_def_var(ncid, "jz", NF90_INT, var_jz_id))
  call check(nf90_def_var(ncid, "kp", NF90_INT, var_kp_id))
  call check(nf90_def_var(ncid, "nfp", NF90_INT, var_nfp_id))

  call check(nf90_enddef(ncid))  ! End define mode to write scalar variables
  
  ! ----------------------------------
  ! 7. Write Data
  ! ----------------------------------
  
  ! Write scalars
  call check(nf90_put_var(ncid, var_rmax_id, R_max))
  call check(nf90_put_var(ncid, var_rmin_id, R_min))
  call check(nf90_put_var(ncid, var_zmax_id, Z_max))
  call check(nf90_put_var(ncid, var_zmin_id, Z_min))

  ! Write integer parameters
  call check(nf90_put_var(ncid, var_ir_id, nR))
  call check(nf90_put_var(ncid, var_jz_id, nZ))
  call check(nf90_put_var(ncid, var_kp_id, n_planes_out))
  call check(nf90_put_var(ncid, var_nfp_id, nfp))
  
  ! Write pressure array
  ! Fortran data is (nR, nZ, n_planes_out). 
  ! NetCDF (via the def_var IDs) expects (nR, nZ, n_planes_out).
  ! This results in a file structure of (phi, zee, rad) when viewed in ncdump/python.
  call check(nf90_put_var(ncid, var_pres_id, pressure))
  call check(nf90_put_var(ncid, var_mag_pres_id, magnetic_pressure))
  call check(nf90_put_var(ncid, var_br_id, BR))
  call check(nf90_put_var(ncid, var_bz_id, BZ))
  call check(nf90_put_var(ncid, var_bp_id, Bphi))
  call check(nf90_put_var(ncid, var_bvac_id, B_vac))
  call check(nf90_put_var(ncid, var_bpert_id, B_pert))

  ! ----------------------------------
  ! Close file
  ! ----------------------------------
  call check(nf90_close(ncid))
  
  write(*,*) 'NetCDF file written successfully!'
  write(*,*) '  Dimensions: R=', nR, ', Z=', nZ, ', phi=', n_planes_out
  write(*,*) '  Field periods: ', nfp
  
! #else
!   ! NetCDF not available - write ASCII file instead
!   write(*,*) 'Warning: NetCDF not available, writing ASCII file instead'
! call save_to_ascii(R_min, R_max, nR, Z_min, Z_max, nZ, n_planes_out, nfp, pressure)
! #endif
  
end subroutine save_to_netcdf


subroutine check(status)
use netcdf
  implicit none
  integer, intent(in) :: status
  if (status /= nf90_noerr) then
    ! print *, "NetCDF Error: ", trim(nf90_strerror(status))
    print *, "NetCDF Error"
    stop
  end if
end subroutine check

 
end program jorek2_fields_orin

