!> Toroidal FFT scatter of the preconditioner element matrices: the transform
!! of one (row, column) basis pair, assembled on n_plane toroidal planes, into
!! the cos/sin-packed element matrix. The channel says where d/dphi acts:
!! nowhere (scatter_fft_to_elm), on the trial function (_n), on the test
!! function (_k) or on both (_kn); on a Fourier mode d/dphi = i n.
module mod_pc_fft_scatter

  implicit none
  private

  public :: scatter_fft_to_elm
  public :: scatter_fft_to_elm_n
  public :: scatter_fft_to_elm_k
  public :: scatter_fft_to_elm_kn
  public :: pc_my_fft

contains

!-----------------------------------------------------------------
!> Scatter FFT output for row i, col j into ELM (no toroidal derivative).
!-----------------------------------------------------------------
subroutine scatter_fft_to_elm(out_fft, i, j, ELM, ndim)

  use mod_parameters, only: n_tor, n_plane

  implicit none

  complex*16, intent(in)    :: out_fft(1:n_plane)
  integer,    intent(in)    :: i, j, ndim
  real*8,     intent(inout) :: ELM(ndim, ndim)

  integer :: k, m, l, index_k, index_m

  do k = 1, (n_tor+1)/2

    index_k = n_tor*(i-1) + max(2*(k-1), 1)

    do m = 1, (n_tor+1)/2

      index_m = n_tor*(j-1) + max(2*(m-1), 1)

      l = (k-1) + (m-1)
      if (l .ge. 0 .and. l .le. n_plane/2) then
        ELM(index_k,   index_m  ) = ELM(index_k,   index_m  ) +  real(out_fft(l+1))
        ELM(index_k+1, index_m  ) = ELM(index_k+1, index_m  ) - imag(out_fft(l+1))
        ELM(index_k,   index_m+1) = ELM(index_k,   index_m+1) - imag(out_fft(l+1))
        ELM(index_k+1, index_m+1) = ELM(index_k+1, index_m+1) -  real(out_fft(l+1))
      elseif (l .lt. 0 .and. abs(l) .le. n_plane/2) then
        ELM(index_k,   index_m  ) = ELM(index_k,   index_m  ) +  real(out_fft(abs(l)+1))
        ELM(index_k+1, index_m  ) = ELM(index_k+1, index_m  ) + imag(out_fft(abs(l)+1))
        ELM(index_k,   index_m+1) = ELM(index_k,   index_m+1) + imag(out_fft(abs(l)+1))
        ELM(index_k+1, index_m+1) = ELM(index_k+1, index_m+1) -  real(out_fft(abs(l)+1))
      endif

      l = (k-1) - (m-1)
      if (l .ge. 0 .and. l .le. n_plane/2) then
        ELM(index_k,   index_m  ) = ELM(index_k,   index_m  ) +  real(out_fft(l+1))
        ELM(index_k+1, index_m  ) = ELM(index_k+1, index_m  ) - imag(out_fft(l+1))
        ELM(index_k,   index_m+1) = ELM(index_k,   index_m+1) + imag(out_fft(l+1))
        ELM(index_k+1, index_m+1) = ELM(index_k+1, index_m+1) +  real(out_fft(l+1))
      elseif (l .lt. 0 .and. abs(l) .le. n_plane/2) then
        ELM(index_k,   index_m  ) = ELM(index_k,   index_m  ) +  real(out_fft(abs(l)+1))
        ELM(index_k+1, index_m  ) = ELM(index_k+1, index_m  ) + imag(out_fft(abs(l)+1))
        ELM(index_k,   index_m+1) = ELM(index_k,   index_m+1) - imag(out_fft(abs(l)+1))
        ELM(index_k+1, index_m+1) = ELM(index_k+1, index_m+1) +  real(out_fft(abs(l)+1))
      endif

    enddo  ! m
  enddo    ! k

end subroutine scatter_fft_to_elm

!-----------------------------------------------------------------
!> ELM_n channel: d/dphi acts on the trial (column) function.
!! Weight float(mode(im)) with im the toroidal index of the column.
!-----------------------------------------------------------------
subroutine scatter_fft_to_elm_n(out_fft, i, j, ELM, ndim)

  use mod_parameters, only: n_tor, n_plane
  use phys_module,    only: mode

  implicit none

  complex*16, intent(in)    :: out_fft(1:n_plane)
  integer,    intent(in)    :: i, j, ndim
  real*8,     intent(inout) :: ELM(ndim, ndim)

  integer :: k, m, im, l, index_k, index_m

  do k = 1, (n_tor+1)/2

    index_k = n_tor*(i-1) + max(2*(k-1), 1)

    do m = 1, (n_tor+1)/2

      im      = max(2*(m-1), 1)
      index_m = n_tor*(j-1) + im

      l = (k-1) + (m-1)
      if (l .ge. 0 .and. l .le. n_plane/2) then
        ELM(index_k,   index_m  ) = ELM(index_k,   index_m  ) + imag(out_fft(l+1))        * float(mode(im))
        ELM(index_k+1, index_m  ) = ELM(index_k+1, index_m  ) + real(out_fft(l+1))        * float(mode(im))
        ELM(index_k,   index_m+1) = ELM(index_k,   index_m+1) + real(out_fft(l+1))        * float(mode(im))
        ELM(index_k+1, index_m+1) = ELM(index_k+1, index_m+1) - imag(out_fft(l+1))        * float(mode(im))
      elseif (l .lt. 0 .and. abs(l) .le. n_plane/2) then
        ELM(index_k,   index_m  ) = ELM(index_k,   index_m  ) - imag(out_fft(abs(l)+1))   * float(mode(im))
        ELM(index_k+1, index_m  ) = ELM(index_k+1, index_m  ) + real(out_fft(abs(l)+1))   * float(mode(im))
        ELM(index_k,   index_m+1) = ELM(index_k,   index_m+1) + real(out_fft(abs(l)+1))   * float(mode(im))
        ELM(index_k+1, index_m+1) = ELM(index_k+1, index_m+1) + imag(out_fft(abs(l)+1))   * float(mode(im))
      endif

      l = (k-1) - (m-1)
      if (l .ge. 0 .and. l .le. n_plane/2) then
        ELM(index_k,   index_m  ) = ELM(index_k,   index_m  ) - imag(out_fft(l+1))        * float(mode(im))
        ELM(index_k+1, index_m  ) = ELM(index_k+1, index_m  ) - real(out_fft(l+1))        * float(mode(im))
        ELM(index_k,   index_m+1) = ELM(index_k,   index_m+1) + real(out_fft(l+1))        * float(mode(im))
        ELM(index_k+1, index_m+1) = ELM(index_k+1, index_m+1) - imag(out_fft(l+1))        * float(mode(im))
      elseif (l .lt. 0 .and. abs(l) .le. n_plane/2) then
        ELM(index_k,   index_m  ) = ELM(index_k,   index_m  ) + imag(out_fft(abs(l)+1))   * float(mode(im))
        ELM(index_k+1, index_m  ) = ELM(index_k+1, index_m  ) - real(out_fft(abs(l)+1))   * float(mode(im))
        ELM(index_k,   index_m+1) = ELM(index_k,   index_m+1) + real(out_fft(abs(l)+1))   * float(mode(im))
        ELM(index_k+1, index_m+1) = ELM(index_k+1, index_m+1) + imag(out_fft(abs(l)+1))   * float(mode(im))
      endif

    enddo  ! m
  enddo    ! k

end subroutine scatter_fft_to_elm_n

!-----------------------------------------------------------------
!> ELM_k channel: d/dphi acts on the test (row) function.
!! Weight float(mode(ik)) with ik the toroidal index of the row.
!-----------------------------------------------------------------
subroutine scatter_fft_to_elm_k(out_fft, i, j, ELM, ndim)

  use mod_parameters, only: n_tor, n_plane
  use phys_module,    only: mode

  implicit none

  complex*16, intent(in)    :: out_fft(1:n_plane)
  integer,    intent(in)    :: i, j, ndim
  real*8,     intent(inout) :: ELM(ndim, ndim)

  integer :: k, m, ik, l, index_k, index_m

  do k = 1, (n_tor+1)/2

    ik      = max(2*(k-1), 1)
    index_k = n_tor*(i-1) + ik

    do m = 1, (n_tor+1)/2

      index_m = n_tor*(j-1) + max(2*(m-1), 1)

      l = (k-1) + (m-1)
      if (l .ge. 0 .and. l .le. n_plane/2) then
        ELM(index_k,   index_m  ) = ELM(index_k,   index_m  ) + imag(out_fft(l+1))        * float(mode(ik))
        ELM(index_k+1, index_m  ) = ELM(index_k+1, index_m  ) + real(out_fft(l+1))        * float(mode(ik))
        ELM(index_k,   index_m+1) = ELM(index_k,   index_m+1) + real(out_fft(l+1))        * float(mode(ik))
        ELM(index_k+1, index_m+1) = ELM(index_k+1, index_m+1) - imag(out_fft(l+1))        * float(mode(ik))
      elseif (l .lt. 0 .and. abs(l) .le. n_plane/2) then
        ELM(index_k,   index_m  ) = ELM(index_k,   index_m  ) - imag(out_fft(abs(l)+1))   * float(mode(ik))
        ELM(index_k+1, index_m  ) = ELM(index_k+1, index_m  ) + real(out_fft(abs(l)+1))   * float(mode(ik))
        ELM(index_k,   index_m+1) = ELM(index_k,   index_m+1) + real(out_fft(abs(l)+1))   * float(mode(ik))
        ELM(index_k+1, index_m+1) = ELM(index_k+1, index_m+1) + imag(out_fft(abs(l)+1))   * float(mode(ik))
      endif

      l = (k-1) - (m-1)
      if (l .ge. 0 .and. l .le. n_plane/2) then
        ELM(index_k,   index_m  ) = ELM(index_k,   index_m  ) + imag(out_fft(l+1))        * float(mode(ik))
        ELM(index_k+1, index_m  ) = ELM(index_k+1, index_m  ) + real(out_fft(l+1))        * float(mode(ik))
        ELM(index_k,   index_m+1) = ELM(index_k,   index_m+1) - real(out_fft(l+1))        * float(mode(ik))
        ELM(index_k+1, index_m+1) = ELM(index_k+1, index_m+1) + imag(out_fft(l+1))        * float(mode(ik))
      elseif (l .lt. 0 .and. abs(l) .le. n_plane/2) then
        ELM(index_k,   index_m  ) = ELM(index_k,   index_m  ) - imag(out_fft(abs(l)+1))   * float(mode(ik))
        ELM(index_k+1, index_m  ) = ELM(index_k+1, index_m  ) + real(out_fft(abs(l)+1))   * float(mode(ik))
        ELM(index_k,   index_m+1) = ELM(index_k,   index_m+1) - real(out_fft(abs(l)+1))   * float(mode(ik))
        ELM(index_k+1, index_m+1) = ELM(index_k+1, index_m+1) - imag(out_fft(abs(l)+1))   * float(mode(ik))
      endif

    enddo  ! m
  enddo    ! k

end subroutine scatter_fft_to_elm_k

!-----------------------------------------------------------------
!> ELM_kn channel: d/dphi acts on BOTH test and trial functions.
!! Weight float(mode(ik))*float(mode(im)); the extra sign is i^2 = -1.
!-----------------------------------------------------------------
subroutine scatter_fft_to_elm_kn(out_fft, i, j, ELM, ndim)

  use mod_parameters, only: n_tor, n_plane
  use phys_module,    only: mode

  implicit none

  complex*16, intent(in)    :: out_fft(1:n_plane)
  integer,    intent(in)    :: i, j, ndim
  real*8,     intent(inout) :: ELM(ndim, ndim)

  integer :: k, m, ik, im, l, index_k, index_m

  do k = 1, (n_tor+1)/2

    ik      = max(2*(k-1), 1)
    index_k = n_tor*(i-1) + ik

    do m = 1, (n_tor+1)/2

      im      = max(2*(m-1), 1)
      index_m = n_tor*(j-1) + im

      l = (k-1) + (m-1)
      if (l .ge. 0 .and. l .le. n_plane/2) then
        ELM(index_k,   index_m  ) = ELM(index_k,   index_m  ) - real(out_fft(l+1))        * float(mode(im)) * float(mode(ik))
        ELM(index_k+1, index_m  ) = ELM(index_k+1, index_m  ) + imag(out_fft(l+1))        * float(mode(im)) * float(mode(ik))
        ELM(index_k,   index_m+1) = ELM(index_k,   index_m+1) + imag(out_fft(l+1))        * float(mode(im)) * float(mode(ik))
        ELM(index_k+1, index_m+1) = ELM(index_k+1, index_m+1) + real(out_fft(l+1))        * float(mode(im)) * float(mode(ik))
      elseif (l .lt. 0 .and. abs(l) .le. n_plane/2) then
        ELM(index_k,   index_m  ) = ELM(index_k,   index_m  ) - real(out_fft(abs(l)+1))   * float(mode(im)) * float(mode(ik))
        ELM(index_k+1, index_m  ) = ELM(index_k+1, index_m  ) - imag(out_fft(abs(l)+1))   * float(mode(im)) * float(mode(ik))
        ELM(index_k,   index_m+1) = ELM(index_k,   index_m+1) - imag(out_fft(abs(l)+1))   * float(mode(im)) * float(mode(ik))
        ELM(index_k+1, index_m+1) = ELM(index_k+1, index_m+1) + real(out_fft(abs(l)+1))   * float(mode(im)) * float(mode(ik))
      endif

      l = (k-1) - (m-1)
      if (l .ge. 0 .and. l .le. n_plane/2) then
        ELM(index_k,   index_m  ) = ELM(index_k,   index_m  ) + real(out_fft(l+1))        * float(mode(im)) * float(mode(ik))
        ELM(index_k+1, index_m  ) = ELM(index_k+1, index_m  ) - imag(out_fft(l+1))        * float(mode(im)) * float(mode(ik))
        ELM(index_k,   index_m+1) = ELM(index_k,   index_m+1) + imag(out_fft(l+1))        * float(mode(im)) * float(mode(ik))
        ELM(index_k+1, index_m+1) = ELM(index_k+1, index_m+1) + real(out_fft(l+1))        * float(mode(im)) * float(mode(ik))
      elseif (l .lt. 0 .and. abs(l) .le. n_plane/2) then
        ELM(index_k,   index_m  ) = ELM(index_k,   index_m  ) + real(out_fft(abs(l)+1))   * float(mode(im)) * float(mode(ik))
        ELM(index_k+1, index_m  ) = ELM(index_k+1, index_m  ) + imag(out_fft(abs(l)+1))   * float(mode(im)) * float(mode(ik))
        ELM(index_k,   index_m+1) = ELM(index_k,   index_m+1) - imag(out_fft(abs(l)+1))   * float(mode(im)) * float(mode(ik))
        ELM(index_k+1, index_m+1) = ELM(index_k+1, index_m+1) + real(out_fft(abs(l)+1))   * float(mode(im)) * float(mode(ik))
      endif

    enddo  ! m
  enddo    ! k

end subroutine scatter_fft_to_elm_kn

!-----------------------------------------------------------------
!> Real-to-complex transform of one n-plane channel without FFTW (RFT2), the
!! fallback the element matrices take when JOREK is built without USE_FFTW.
!! Same output layout as dfftw_execute_dft_r2c: out_fft(k) holds harmonic k-1.
!-----------------------------------------------------------------
subroutine pc_my_fft(in_fft, out_fft, n)

  implicit none

  integer,    intent(in)  :: n
  real*8,     intent(in)  :: in_fft(n)
  complex*16, intent(out) :: out_fft(*)

  real*8  :: tmp_fft(2*n+2)
  integer :: i

  tmp_fft(1:n) = in_fft(1:n)
  call RFT2(tmp_fft, n, 1)
  do i = 1, n/2 + 1
    out_fft(i) = cmplx(tmp_fft(2*i-1), tmp_fft(2*i), kind=8)
  enddo

end subroutine pc_my_fft

end module mod_pc_fft_scatter
