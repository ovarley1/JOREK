module mod_petsc_raw_csr
#ifdef USE_PETSC
  use iso_c_binding
#include "petsc/finclude/petsc.h"
  use petsc
  implicit none
  private
#if defined(PETSC_USE_COMPLEX) || !defined(PETSC_USE_REAL_DOUBLE)
#error "mod_petsc_raw_csr: the value maps are real(c_double)"
#endif

  !> Raw CSR access to PETSc (MPI)AIJ/BAIJ matrices and the C kernels, via iso_c_binding (the Fortran
  !! MatMPIBAIJGetSeqBAIJ mishandles its colmap). Any failure ends the run on all ranks.

  interface
    integer(c_int) function c_mpibaij_seq(A, Ad, Ao, cmap) bind(C, name="MatMPIBAIJGetSeqBAIJ")
      use, intrinsic :: iso_c_binding, only: c_int, c_intptr_t, c_ptr
      integer(c_intptr_t), value :: A
      integer(c_intptr_t)        :: Ad, Ao
      type(c_ptr)                :: cmap
    end function c_mpibaij_seq
    integer(c_int) function c_mpiaij_seq(A, Ad, Ao, cmap) bind(C, name="MatMPIAIJGetSeqAIJ")
      use, intrinsic :: iso_c_binding, only: c_int, c_intptr_t, c_ptr
      integer(c_intptr_t), value :: A
      integer(c_intptr_t)        :: Ad, Ao
      type(c_ptr)                :: cmap
    end function c_mpiaij_seq
    integer(c_int) function c_baij_get(M, arr) bind(C, name="MatSeqBAIJGetArray")
      use, intrinsic :: iso_c_binding, only: c_int, c_intptr_t, c_ptr
      integer(c_intptr_t), value :: M
      type(c_ptr)                :: arr
    end function c_baij_get
    integer(c_int) function c_baij_restore(M, arr) bind(C, name="MatSeqBAIJRestoreArray")
      use, intrinsic :: iso_c_binding, only: c_int, c_intptr_t, c_ptr
      integer(c_intptr_t), value :: M
      type(c_ptr)                :: arr
    end function c_baij_restore
    integer(c_int) function c_aij_get(M, arr) bind(C, name="MatSeqAIJGetArray")
      use, intrinsic :: iso_c_binding, only: c_int, c_intptr_t, c_ptr
      integer(c_intptr_t), value :: M
      type(c_ptr)                :: arr
    end function c_aij_get
    integer(c_int) function c_aij_restore(M, arr) bind(C, name="MatSeqAIJRestoreArray")
      use, intrinsic :: iso_c_binding, only: c_int, c_intptr_t, c_ptr
      integer(c_intptr_t), value :: M
      type(c_ptr)                :: arr
    end function c_aij_restore
    ! read-only access: does not advance the matrix's object state
    integer(c_int) function c_aij_get_read(M, arr) bind(C, name="MatSeqAIJGetArrayRead")
      use, intrinsic :: iso_c_binding, only: c_int, c_intptr_t, c_ptr
      integer(c_intptr_t), value :: M
      type(c_ptr)                :: arr
    end function c_aij_get_read
    integer(c_int) function c_aij_restore_read(M, arr) bind(C, name="MatSeqAIJRestoreArrayRead")
      use, intrinsic :: iso_c_binding, only: c_int, c_intptr_t, c_ptr
      integer(c_intptr_t), value :: M
      type(c_ptr)                :: arr
    end function c_aij_restore_read
    ! the threaded block matvec of jorek_blockmv_attach.c (this declaration's
    ! form is what util/makedepend links the C file by)
    function jorek_blockmv_attach(M, bs, ok) bind(C)
      use, intrinsic :: iso_c_binding, only: c_int, c_intptr_t
      integer(c_int)             :: jorek_blockmv_attach
      integer(c_intptr_t), value :: M
      integer(c_int), value      :: bs
      integer(c_int)             :: ok
    end function jorek_blockmv_attach
    ! the Galerkin product on dense slot blocks of jorek_galerkin_blk.c
    function jorek_galerkin_blk(A, Ps, bs, reuse, C, ok) bind(C)
      use, intrinsic :: iso_c_binding, only: c_int, c_intptr_t
      integer(c_int)             :: jorek_galerkin_blk
      integer(c_intptr_t), value :: A, Ps
      integer(c_int), value      :: bs, reuse
      integer(c_intptr_t)        :: C
      integer(c_int)             :: ok
    end function jorek_galerkin_blk
  end interface

  public :: c_baij_get, c_baij_restore, c_aij_get, c_aij_restore
  public :: split_parts, get_ij, put_ij, aij_parts, aij_vals_read, aij_vals_done
  public :: blockmv_attach, galerkin_blk

contains

  !> Ends the run on all ranks (a PETSc error may come from one rank only, inside collective work).
  subroutine fatal(what)
    character(len=*), intent(in) :: what
    integer :: ierr
    write(*,'(A,A)') "[Physics PC]   FATAL: ", what
    call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
  end subroutine fatal

  !> fatal() unless rc (a PETSc error code) is 0; PETSc has printed the error on this rank.
  subroutine chk(rc, what)
    integer(c_int), intent(in) :: rc
    character(len=*), intent(in) :: what
    character(len=16) :: num
    if (rc == 0) return
    write(num, '(I0)') rc
    call fatal(what//": PETSc error "//trim(num))
  end subroutine chk

  !> M's MatMult/MatMultAdd on the OpenMP kernel of jorek_blockmv_attach.c (.false.: PETSc's kept); collective, idempotent.
  logical function blockmv_attach(M, bs)
    Mat, intent(in)     :: M
    integer, intent(in) :: bs
    integer(c_int) :: ok
    ok = 0
    call chk(jorek_blockmv_attach(transfer(M%v, 0_c_intptr_t), int(bs, c_int), ok), "blockmv_attach")
    blockmv_attach = (ok /= 0)
  end function blockmv_attach

  !> C = P^T A P on dense slot blocks (jorek_galerkin_blk.c), reuse refills; .false. on all ranks if A / Ps do not fit.
  logical function galerkin_blk(A, Ps, bs, reuse, C)
    Mat, intent(in)     :: A, Ps
    integer, intent(in) :: bs
    logical, intent(in) :: reuse
    Mat, intent(inout)  :: C
    integer(c_int) :: ok
    integer(c_intptr_t) :: cv
    cv = 0
    if (reuse) cv = transfer(C%v, 0_c_intptr_t)
    ok = 0
    call chk(jorek_galerkin_blk(transfer(A%v, 0_c_intptr_t), transfer(Ps%v, 0_c_intptr_t), int(bs, c_int), &
                                merge(1_c_int, 0_c_int, reuse), cv, ok), "galerkin_blk")
    galerkin_blk = (ok /= 0)
    if (galerkin_blk .and. .not. reuse) C%v = transfer(cv, C%v)
  end function galerkin_blk

  !> Diagonal/off-diagonal sequential parts of an MPI(B)AIJ matrix, and the
  !! off-diagonal part's global (block) column of each local column.
  subroutine split_parts(M, is_baij, Md, Mo, garr)
    Mat, intent(in)  :: M
    logical, intent(in) :: is_baij
    Mat, intent(out) :: Md, Mo
    PetscInt, pointer, intent(out) :: garr(:)
    integer(c_intptr_t) :: pd, po
    type(c_ptr) :: cm
    PetscInt :: nr, nc, bs
    PetscErrorCode :: ierr
    if (is_baij) then
      call chk(c_mpibaij_seq(transfer(M%v, 0_c_intptr_t), pd, po, cm), "split_parts: MatMPIBAIJGetSeqBAIJ")
    else
      call chk(c_mpiaij_seq(transfer(M%v, 0_c_intptr_t), pd, po, cm), "split_parts: MatMPIAIJGetSeqAIJ")
    endif
    Md%v = transfer(pd, Md%v); Mo%v = transfer(po, Mo%v)
    call MatGetSize(Mo, nr, nc, ierr)
    call chk(ierr, "split_parts: MatGetSize")
    bs = 1
    if (is_baij) then
      call MatGetBlockSize(Mo, bs, ierr)
      call chk(ierr, "split_parts: MatGetBlockSize")
    endif
    if (nc / bs > 0) then
      call c_f_pointer(cm, garr, [nc / bs])
    else
      garr => null()
    endif
  end subroutine split_parts

  !> split_parts for an AIJ matrix that may also be sequential (a one-rank
  !! communicator can give either): has_o = .false. means M is its own
  !! diagonal part and there is no off-diagonal part. The AIJ families as
  !! MatMPIAIJGetSeqAIJ takes them (types beginning mpiaij / seqaij).
  subroutine aij_parts(M, Md, Mo, garr, has_o)
    Mat, intent(in)  :: M
    Mat, intent(out) :: Md, Mo
    PetscInt, pointer, intent(out) :: garr(:)
    logical, intent(out) :: has_o
    MatType :: mt
    PetscErrorCode :: ierr
    call MatGetType(M, mt, ierr)
    call chk(ierr, "aij_parts: MatGetType")
    has_o = (index(mt, "mpiaij") == 1)
    if (has_o) then
      call split_parts(M, .false., Md, Mo, garr)
    else if (index(mt, "seqaij") == 1) then
      Md = M; garr => null()
    else
      call fatal("raw CSR access needs an AIJ matrix, got "//trim(mt))
    endif
  end subroutine aij_parts

  !> Block CSR (BAIJ) or CSR (AIJ) of a sequential part, 0-based.
  subroutine get_ij(M, blocked, n, ia, ja)
    Mat, intent(in) :: M
    logical, intent(in) :: blocked
    PetscInt, intent(out) :: n
    PetscInt, pointer :: ia(:), ja(:)
    PetscBool :: done, bc
    PetscErrorCode :: ierr
    PetscInt, parameter :: zero = 0
    bc = PETSC_FALSE
    if (blocked) bc = PETSC_TRUE
    call MatGetRowIJ(M, zero, PETSC_FALSE, bc, n, ia, ja, done, ierr)
    call chk(ierr, "get_ij: MatGetRowIJ")
    if (.not. done) call fatal("get_ij: MatGetRowIJ gave no row IJ")
  end subroutine get_ij

  subroutine put_ij(M, blocked, n, ia, ja)
    Mat, intent(in) :: M
    logical, intent(in) :: blocked
    PetscInt, intent(inout) :: n
    PetscInt, pointer :: ia(:), ja(:)
    PetscBool :: done, bc
    PetscErrorCode :: ierr
    PetscInt, parameter :: zero = 0
    bc = PETSC_FALSE
    if (blocked) bc = PETSC_TRUE
    call MatRestoreRowIJ(M, zero, PETSC_FALSE, bc, n, ia, ja, done, ierr)
    call chk(ierr, "put_ij: MatRestoreRowIJ")
    if (.not. done) call fatal("put_ij: MatRestoreRowIJ failed")
  end subroutine put_ij

  !> Read-only value array (length n) of a sequential AIJ matrix; p is the
  !! handle aij_vals_done needs.
  subroutine aij_vals_read(M, n, p, v)
    Mat, intent(in) :: M
    integer(8), intent(in) :: n
    type(c_ptr), intent(out) :: p
    real(c_double), pointer, intent(out) :: v(:)
    call chk(c_aij_get_read(transfer(M%v, 0_c_intptr_t), p), "aij_vals_read: MatSeqAIJGetArrayRead")
    call c_f_pointer(p, v, [max(n, 1_8)])
  end subroutine aij_vals_read

  subroutine aij_vals_done(M, p)
    Mat, intent(in) :: M
    type(c_ptr), intent(inout) :: p
    call chk(c_aij_restore_read(transfer(M%v, 0_c_intptr_t), p), "aij_vals_done: MatSeqAIJRestoreArrayRead")
  end subroutine aij_vals_done

#endif
end module mod_petsc_raw_csr
