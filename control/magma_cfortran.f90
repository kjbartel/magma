!
!   -- MAGMA (version 1.3.0) --
!      Univ. of Tennessee, Knoxville
!      Univ. of California, Berkeley
!      Univ. of Colorado, Denver
!      November 2012
!
!   @generated c Wed Nov 14 22:52:30 2012
!

#define PRECISION_c

module magma_cfortran

  use magma_param, only: sizeof_complex

  implicit none

  !---- Fortran interfaces to MAGMA subroutines ----
  interface

     subroutine magmaf_cgetptr( m, n, A, lda, d, e,tauq, taup, work, lwork, info)
       integer       :: m
       integer       :: n
       complex    :: A(*)
       integer       :: lda
       real:: d(*)
       real:: e(*)
       complex    :: tauq(*)
       complex    :: taup(*)
       complex    :: work(*)
       integer       :: lwork
       integer       :: info
     end subroutine magmaf_cgetptr


     subroutine magmaf_cgebrd( m, n, A, lda, d, e,tauq, taup, work, lwork, info)
       integer       :: m
       integer       :: n
       complex    :: A(*)
       integer       :: lda
       real:: d(*)
       real:: e(*)
       complex    :: tauq(*)
       complex    :: taup(*)
       complex    :: work(*)
       integer       :: lwork
       integer       :: info
     end subroutine magmaf_cgebrd

     subroutine magmaf_cgehrd2(n, ilo, ihi,A, lda, tau, work, lwork, info)
       integer       :: n
       integer       :: ilo
       integer       :: ihi
       complex    :: A(*)
       integer       :: lda
       complex    :: tau(*)
       complex    :: work(*)
       integer       :: lwork
       integer       :: info
     end subroutine magmaf_cgehrd2

     subroutine magmaf_cgehrd(n, ilo, ihi,A, lda, tau, work, lwork, d_T, info)
       integer       :: n
       integer       :: ilo
       integer       :: ihi
       complex    :: A(*)
       integer       :: lda
       complex    :: tau(*)
       complex    :: work(*)
       integer       :: lwork
       complex    :: d_T(*)
       integer       :: info
     end subroutine magmaf_cgehrd

     subroutine magmaf_cgelqf( m, n, A,    lda,   tau, work, lwork, info)
       integer       :: m
       integer       :: n
       complex    :: A(*)
       integer       :: lda
       complex    :: tau(*)
       complex    :: work(*)
       integer       :: lwork
       integer       :: info
     end subroutine magmaf_cgelqf

     subroutine magmaf_cgeqlf( m, n, A,    lda,   tau, work, lwork, info)
       integer       :: m
       integer       :: n
       complex    :: A(*)
       integer       :: lda
       complex    :: tau(*)
       complex    :: work(*)
       integer       :: lwork
       integer       :: info
     end subroutine magmaf_cgeqlf

     subroutine magmaf_cgeqrf( m, n, A, lda, tau, work, lwork, info)
       integer       :: m
       integer       :: n
       complex    :: A(*)
       integer       :: lda
       complex    :: tau(*)
       complex    :: work(*)
       integer       :: lwork
       integer       :: info
     end subroutine magmaf_cgeqrf

     subroutine magmaf_cgesv(  n, nrhs, A, lda, ipiv, B, ldb, info)
       integer       :: n
       integer       :: nrhs
       complex    :: A
       integer       :: lda
       integer       :: ipiv(*)
       complex    :: B
       integer       :: ldb
       integer       :: info
     end subroutine magmaf_cgesv

     subroutine magmaf_cgetrf( m, n, A, lda, ipiv, info)
       integer       :: m
       integer       :: n
       complex    :: A(*)
       integer       :: lda
       integer       :: ipiv(*)
       integer       :: info
     end subroutine magmaf_cgetrf

     subroutine magmaf_chegst( itype, uplo, n, A, lda, B, ldb, info)
       integer       :: itype
       character     :: uplo
       integer       :: n
       complex    :: A(*)
       integer       :: lda
       complex    :: B(*)
       integer       :: ldb
       integer       :: info
     end subroutine magmaf_chegst

     subroutine magmaf_cposv(  uplo, n, nrhs, dA, ldda, dB, lddb, info)
       character     :: uplo
       integer       :: n
       integer       :: nrhs
       magma_devptr_t:: dA
       integer       :: ldda
       magma_devptr_t:: dB
       integer       :: lddb
       integer       :: info
     end subroutine magmaf_cposv
     
     subroutine magmaf_cpotrf( uplo, n, A, lda, info)
       character          :: uplo
       integer       :: n
       complex    :: A(*)
       integer       :: lda
       integer       :: info
     end subroutine magmaf_cpotrf

     subroutine magmaf_chetrd( uplo, n, A, lda, d, e, tau, work, lwork, info)
       character          :: uplo
       integer       :: n
       complex    :: A(*)
       integer       :: lda
       real:: d(*)
       real:: e(*)
       complex    :: tau(*)
       complex    :: work(*)
       integer       :: lwork
       integer       :: info
     end subroutine magmaf_chetrd

     subroutine magmaf_cunmqr( side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info)
       character          :: side
       character          :: trans
       integer       :: m
       integer       :: n
       integer       :: k
       complex    :: a(*)
       integer       :: lda
       complex    :: tau(*)
       complex    :: c(*)
       integer       :: ldc
       complex    :: work(*)
       integer       :: lwork
       integer       :: info
     end subroutine magmaf_cunmqr

     subroutine magmaf_cunmtr( side, uplo, trans, m, n, a, lda,tau,c,    ldc,work, lwork,info)
       character          :: side
       character          :: uplo
       character          :: trans
       integer       :: m
       integer       :: n
       complex    :: a(*)
       integer       :: lda
       complex    :: tau(*)
       complex    :: c(*)
       integer       :: ldc
       complex    :: work(*)
       integer       :: lwork
       integer       :: info
     end subroutine magmaf_cunmtr
#if defined(PRECISION_z) || defined(PRECISION_c)

     subroutine magmaf_cgeev( jobvl, jobvr, n, a, lda, w, vl, ldvl, vr, ldvr, work, lwork, rwork, info)
       character          :: jobvl
       character          :: jobvr
       integer       :: n
       complex    :: a(*)
       integer       :: lda
       complex    :: w(*)
       complex    :: vl(*)
       integer       :: ldvl
       complex    :: vr(*)
       integer       :: ldvr
       complex    :: work(*)
       integer       :: lwork
       real:: rwork(*)
       integer       :: info
     end subroutine magmaf_cgeev

     subroutine magmaf_cgesvd( jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, rwork, info)
       character          :: jobu
       character          :: jobvt
       integer       :: m
       integer       :: n
       complex    :: a(*)
       integer       :: lda
       real:: s(*)
       complex    :: u(*)
       integer       :: ldu
       complex    :: vt(*)
       integer       :: ldvt
       complex    :: work(*)
       integer       :: lwork
       real:: rwork(*)
       integer       :: info
     end subroutine magmaf_cgesvd

     subroutine magmaf_cheevd( jobz, uplo, n, a, lda, w, work, lwork, rwork, lrwork, iwork, liwork, info)
       character     :: jobz
       character     :: uplo
       integer       :: n
       complex    :: a(*)
       integer       :: lda
       real:: w(*)
       complex    :: work(*)
       integer       :: lwork
       real:: rwork(*)
       integer       :: lrwork
       integer       :: iwork(*)
       integer       :: liwork
       integer       :: info
     end subroutine magmaf_cheevd

     subroutine magmaf_chegvd( itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, rwork, lrwork, iwork, liwork, info)
       integer       :: itype
       character     :: jobz
       character     :: uplo
       integer       :: n
       complex    :: a(*)
       integer       :: lda
       complex    :: b(*)
       integer       :: ldb
       real:: w(*)
       complex    :: work(*)
       integer       :: lwork
       real:: rwork(*)
       integer       :: lrwork
       integer       :: iwork(*)
       integer       :: liwork
       integer       :: info
     end subroutine magmaf_chegvd

#else
     subroutine magmaf_cgeev( jobvl, jobvr, n, a, lda, wr, wi, vl, ldvl, vr, ldvr, work, lwork, info)
       character          :: jobvl
       character          :: jobvr
       integer       :: n
       complex    :: a(*)
       integer       :: lda
       complex    :: wr(*)
       complex    :: wi(*)
       complex    :: vl(*)
       integer       :: ldvl
       complex    :: vr(*)
       integer       :: ldvr
       complex    :: work(*)
       integer       :: lwork
       integer       :: info
     end subroutine magmaf_cgeev

     subroutine magmaf_cgesvd( jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, info)
       character          :: jobu
       character          :: jobvt
       integer       :: m
       integer       :: n
       complex    :: a(*)
       integer       :: lda
       real:: s(*)
       complex    :: u(*)
       integer       :: ldu
       complex    :: vt(*)
       integer       :: ldvt
       complex    :: work(*)
       integer       :: lwork
       integer       :: info
     end subroutine magmaf_cgesvd

     subroutine magmaf_cheevd( jobz, uplo, n, a, lda, w, work, lwork, iwork, liwork, info)
       character          :: jobz
       character          :: uplo
       integer       :: n
       complex    :: a(*)
       integer       :: lda
       real:: w(*)
       complex    :: work(*)
       integer       :: lwork
       integer       :: iwork(*)
       integer       :: liwork
       integer       :: info
     end subroutine magmaf_cheevd

     subroutine magmaf_chegvd( itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, iwork, liwork, info)
       integer       :: itype
       character     :: jobz
       character     :: uplo
       integer       :: n
       complex    :: a(*)
       integer       :: lda
       complex    :: b(*)
       integer       :: ldb
       real:: w(*)
       complex    :: work(*)
       integer       :: lwork
       integer       :: iwork(*)
       integer       :: liwork
       integer       :: info
     end subroutine magmaf_chegvd
#endif

     subroutine magmaf_cgels_gpu(  trans, m, n, nrhs, dA, ldda, dB, lddb, hwork, lwork, info)
       character          :: trans
       integer       :: m
       integer       :: n
       integer       :: nrhs
       magma_devptr_t:: dA
       integer       :: ldda
       magma_devptr_t:: dB
       integer       :: lddb
       complex    :: hwork(*)
       integer       :: lwork
       integer       :: info
     end subroutine magmaf_cgels_gpu

     subroutine magmaf_cgeqrf_gpu( m, n, dA, ldda, tau, dT, info)
       integer       :: m
       integer       :: n
       magma_devptr_t:: dA
       integer       :: ldda
       complex    :: tau(*)
       magma_devptr_t:: dT
       integer       :: info
     end subroutine magmaf_cgeqrf_gpu

     subroutine magmaf_cgeqrf2_gpu(m, n, dA, ldda, tau, info)
       integer       :: m
       integer       :: n
       magma_devptr_t:: dA
       integer       :: ldda
       complex    :: tau(*)
       integer       :: info
     end subroutine magmaf_cgeqrf2_gpu

     subroutine magmaf_cgeqrf3_gpu(m, n, dA, ldda, tau, dT, info)
       integer       :: m
       integer       :: n
       magma_devptr_t:: dA
       integer       :: ldda
       complex    :: tau(*)
       magma_devptr_t:: dT
       integer       :: info
     end subroutine magmaf_cgeqrf3_gpu

     subroutine magmaf_cgeqrs_gpu( m, n, nrhs, dA, ldda, tau, dT, dB, lddb, hwork, lhwork, info)
       integer       :: m
       integer       :: n
       integer       :: nrhs
       magma_devptr_t:: dA
       integer       :: ldda
       complex    :: tau
       magma_devptr_t:: dT
       magma_devptr_t:: dB
       integer       :: lddb
       complex    :: hwork(*)
       integer       :: lhwork
       integer       :: info
     end subroutine magmaf_cgeqrs_gpu

     subroutine magmaf_cgeqrs3_gpu( m, n, nrhs, dA, ldda, tau, dT, dB, lddb, hwork, lhwork, info)
       integer       :: m
       integer       :: n
       integer       :: nrhs
       magma_devptr_t:: dA
       integer       :: ldda
       complex    :: tau
       magma_devptr_t:: dT
       magma_devptr_t:: dB
       integer       :: lddb
       complex    :: hwork(*)
       integer       :: lhwork
       integer       :: info
     end subroutine magmaf_cgeqrs3_gpu

     subroutine magmaf_cgessm_gpu( storev, m, n, k, ib, ipiv, dL1, lddl1, dL,  lddl, dA,  ldda, info)
       character          :: storev
       integer       :: m
       integer       :: n
       integer       :: k
       integer       :: ib
       integer       :: ipiv(*)
       magma_devptr_t:: dL1
       integer       :: lddl1
       magma_devptr_t:: dL
       integer       :: lddl
       magma_devptr_t:: dA
       integer       :: ldda
       integer       :: info
     end subroutine magmaf_cgessm_gpu

     subroutine magmaf_cgesv_gpu(  n, nrhs, dA, ldda, ipiv, dB, lddb, info)
       integer       :: n
       integer       :: nrhs
       magma_devptr_t:: dA
       integer       :: ldda
       integer       :: ipiv(*)
       magma_devptr_t:: dB
       integer       :: lddb
       integer       :: info
     end subroutine magmaf_cgesv_gpu

     subroutine magmaf_cgetrf_gpu( m, n, dA, ldda, ipiv, info)
       integer       :: m
       integer       :: n
       magma_devptr_t:: dA
       integer       :: ldda
       integer       :: ipiv(*)
       integer       :: info
     end subroutine magmaf_cgetrf_gpu

     subroutine magmaf_cgetrs_gpu( trans, n, nrhs, dA, ldda, ipiv, dB, lddb, info)
       character          :: trans
       integer       :: n
       integer       :: nrhs
       magma_devptr_t:: dA
       integer       :: ldda
       integer       :: ipiv(*)
       magma_devptr_t:: dB
       integer       :: lddb
       integer       :: info
     end subroutine magmaf_cgetrs_gpu

     subroutine magmaf_chegst_gpu( itype, uplo, n, dA, ldda, dB, lddb, info)
       integer       :: itype
       character     :: uplo
       integer       :: n
       magma_devptr_t:: dA(*)
       integer       :: ldda
       magma_devptr_t:: dB(*)
       integer       :: lddb
       integer       :: info
     end subroutine magmaf_chegst_gpu

     subroutine magmaf_clabrd_gpu( m, n, nb, a, lda, da, ldda, d, e, tauq, taup, x, ldx, dx, lddx, y, ldy, dy, lddy)
       integer       :: m
       integer       :: n
       integer       :: nb
       complex    :: a(*)
       integer       :: lda
       magma_devptr_t:: da
       integer       :: ldda
       real:: d(*)
       real:: e(*)
       complex    :: tauq(*)
       complex    :: taup(*)
       complex    :: x(*)
       integer       :: ldx
       magma_devptr_t:: dx
       integer       :: lddx
       complex    :: y(*)
       integer       :: ldy
       magma_devptr_t:: dy
       integer       :: lddy
     end subroutine magmaf_clabrd_gpu

     subroutine magmaf_clarfb_gpu( side, trans, direct, storev, m, n, k, dv, ldv, dt, ldt, dc, ldc, dowrk, ldwork)
       character          :: side
       character          :: trans
       character          :: direct
       character          :: storev
       integer       :: m
       integer       :: n
       integer       :: k
       magma_devptr_t:: dv
       integer       :: ldv
       magma_devptr_t:: dt
       integer       :: ldt
       magma_devptr_t:: dc
       integer       :: ldc
       magma_devptr_t:: dowrk
       integer       :: ldwork
     end subroutine magmaf_clarfb_gpu

     subroutine magmaf_cposv_gpu(  uplo, n, nrhs, dA, ldda, dB, lddb, info)
       character          :: uplo
       integer       :: n
       integer       :: nrhs
       magma_devptr_t:: dA
       integer       :: ldda
       magma_devptr_t:: dB
       integer       :: lddb
       integer       :: info
     end subroutine magmaf_cposv_gpu

     subroutine magmaf_cpotrf_gpu( uplo, n, dA, ldda, info)
       character          :: uplo
       integer       :: n
       magma_devptr_t:: dA
       integer       :: ldda
       integer       :: info
     end subroutine magmaf_cpotrf_gpu

     subroutine magmaf_cpotrs_gpu( uplo,  n, nrhs, dA, ldda, dB, lddb, info)
       character          :: uplo
       integer       :: n
       integer       :: nrhs
       magma_devptr_t:: dA
       integer       :: ldda
       magma_devptr_t:: dB
       integer       :: lddb
       integer       :: info
     end subroutine magmaf_cpotrs_gpu

     subroutine magmaf_cssssm_gpu( storev, m1, n1, m2, n2, k, ib, dA1, ldda1, dA2, ldda2, dL1, lddl1, dL2, lddl2, IPIV, info)
       character          :: storev
       integer       :: m1
       integer       :: n1
       integer       :: m2
       integer       :: n2
       integer       :: k
       integer       :: ib
       magma_devptr_t:: dA1
       integer       :: ldda1
       magma_devptr_t:: dA2
       integer       :: ldda2
       magma_devptr_t:: dL1
       integer       :: lddl1
       magma_devptr_t:: dL2
       integer       :: lddl2
       integer       :: IPIV(*)
       integer       :: info
     end subroutine magmaf_cssssm_gpu

     subroutine magmaf_cungqr_gpu( m, n, k, da, ldda, tau, dwork, nb, info)
       integer       :: m
       integer       :: n
       integer       :: k
       magma_devptr_t:: da
       integer       :: ldda
       complex    :: tau(*)
       magma_devptr_t:: dwork
       integer       :: nb
       integer       :: info
     end subroutine magmaf_cungqr_gpu

     subroutine magmaf_cunmqr_gpu( side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, td, nb, info)
       character          :: side
       character          :: trans
       integer       :: m
       integer       :: n
       integer       :: k
       magma_devptr_t:: a
       integer       :: lda
       complex    :: tau(*)
       magma_devptr_t:: c
       integer       :: ldc
       magma_devptr_t:: work
       integer       :: lwork
       magma_devptr_t:: td
       integer       :: nb
       integer       :: info
     end subroutine magmaf_cunmqr_gpu

  end interface

contains
  
  subroutine magmaf_coff1d( ptrNew, ptrOld, inc, i)
    magma_devptr_t :: ptrNew
    magma_devptr_t :: ptrOld
    integer        :: inc, i
    
    ptrNew = ptrOld + (i-1) * inc * sizeof_complex
    
  end subroutine magmaf_coff1d
  
  subroutine magmaf_coff2d( ptrNew, ptrOld, lda, i, j)
    magma_devptr_t :: ptrNew
    magma_devptr_t :: ptrOld
    integer        :: lda, i, j
    
    ptrNew = ptrOld + ((j-1) * lda + (i-1)) * sizeof_complex
    
  end subroutine magmaf_coff2d
  
end module magma_cfortran
