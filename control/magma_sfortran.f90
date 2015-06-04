!
!   -- MAGMA (version 1.2.1) --
!      Univ. of Tennessee, Knoxville
!      Univ. of California, Berkeley
!      Univ. of Colorado, Denver
!      June 2012
!
!   @generated s Thu Jun 28 12:30:06 2012
!

#define PRECISION_s

module magma_sfortran

  use magma_param, only: sizeof_real

  implicit none

  !---- Fortran interfaces to MAGMA subroutines ----
  interface

     subroutine magmaf_sgetptr( m, n, A, lda, d, e,tauq, taup, work, lwork, info)
       integer       :: m
       integer       :: n
       real    :: A(*)
       integer       :: lda
       real:: d(*)
       real:: e(*)
       real    :: tauq(*)
       real    :: taup(*)
       real    :: work(*)
       integer       :: lwork
       integer       :: info
     end subroutine magmaf_sgetptr


     subroutine magmaf_sgebrd( m, n, A, lda, d, e,tauq, taup, work, lwork, info)
       integer       :: m
       integer       :: n
       real    :: A(*)
       integer       :: lda
       real:: d(*)
       real:: e(*)
       real    :: tauq(*)
       real    :: taup(*)
       real    :: work(*)
       integer       :: lwork
       integer       :: info
     end subroutine magmaf_sgebrd

     subroutine magmaf_sgehrd2(n, ilo, ihi,A, lda, tau, work, lwork, info)
       integer       :: n
       integer       :: ilo
       integer       :: ihi
       real    :: A(*)
       integer       :: lda
       real    :: tau(*)
       real    :: work(*)
       integer       :: lwork
       integer       :: info
     end subroutine magmaf_sgehrd2

     subroutine magmaf_sgehrd(n, ilo, ihi,A, lda, tau, work, lwork, d_T, info)
       integer       :: n
       integer       :: ilo
       integer       :: ihi
       real    :: A(*)
       integer       :: lda
       real    :: tau(*)
       real    :: work(*)
       integer       :: lwork
       real    :: d_T(*)
       integer       :: info
     end subroutine magmaf_sgehrd

     subroutine magmaf_sgelqf( m, n, A,    lda,   tau, work, lwork, info)
       integer       :: m
       integer       :: n
       real    :: A(*)
       integer       :: lda
       real    :: tau(*)
       real    :: work(*)
       integer       :: lwork
       integer       :: info
     end subroutine magmaf_sgelqf

     subroutine magmaf_sgeqlf( m, n, A,    lda,   tau, work, lwork, info)
       integer       :: m
       integer       :: n
       real    :: A(*)
       integer       :: lda
       real    :: tau(*)
       real    :: work(*)
       integer       :: lwork
       integer       :: info
     end subroutine magmaf_sgeqlf

     subroutine magmaf_sgeqrf( m, n, A, lda, tau, work, lwork, info)
       integer       :: m
       integer       :: n
       real    :: A(*)
       integer       :: lda
       real    :: tau(*)
       real    :: work(*)
       integer       :: lwork
       integer       :: info
     end subroutine magmaf_sgeqrf

     subroutine magmaf_sgesv(  n, nrhs, A, lda, ipiv, B, ldb, info)
       integer       :: n
       integer       :: nrhs
       real    :: A
       integer       :: lda
       integer       :: ipiv(*)
       real    :: B
       integer       :: ldb
       integer       :: info
     end subroutine magmaf_sgesv

     subroutine magmaf_sgetrf( m, n, A, lda, ipiv, info)
       integer       :: m
       integer       :: n
       real    :: A(*)
       integer       :: lda
       integer       :: ipiv(*)
       integer       :: info
     end subroutine magmaf_sgetrf

     subroutine magmaf_ssygst( itype, uplo, n, A, lda, B, ldb, info)
       integer       :: itype
       character     :: uplo
       integer       :: n
       real    :: A(*)
       integer       :: lda
       real    :: B(*)
       integer       :: ldb
       integer       :: info
     end subroutine magmaf_ssygst

     subroutine magmaf_sposv(  uplo, n, nrhs, dA, ldda, dB, lddb, info)
       character     :: uplo
       integer       :: n
       integer       :: nrhs
       magma_devptr_t:: dA
       integer       :: ldda
       magma_devptr_t:: dB
       integer       :: lddb
       integer       :: info
     end subroutine magmaf_sposv
     
     subroutine magmaf_spotrf( uplo, n, A, lda, info)
       character          :: uplo
       integer       :: n
       real    :: A(*)
       integer       :: lda
       integer       :: info
     end subroutine magmaf_spotrf

     subroutine magmaf_ssytrd( uplo, n, A, lda, d, e, tau, work, lwork, info)
       character          :: uplo
       integer       :: n
       real    :: A(*)
       integer       :: lda
       real:: d(*)
       real:: e(*)
       real    :: tau(*)
       real    :: work(*)
       integer       :: lwork
       integer       :: info
     end subroutine magmaf_ssytrd

     subroutine magmaf_sormqr( side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info)
       character          :: side
       character          :: trans
       integer       :: m
       integer       :: n
       integer       :: k
       real    :: a(*)
       integer       :: lda
       real    :: tau(*)
       real    :: c(*)
       integer       :: ldc
       real    :: work(*)
       integer       :: lwork
       integer       :: info
     end subroutine magmaf_sormqr

     subroutine magmaf_sormtr( side, uplo, trans, m, n, a, lda,tau,c,    ldc,work, lwork,info)
       character          :: side
       character          :: uplo
       character          :: trans
       integer       :: m
       integer       :: n
       real    :: a(*)
       integer       :: lda
       real    :: tau(*)
       real    :: c(*)
       integer       :: ldc
       real    :: work(*)
       integer       :: lwork
       integer       :: info
     end subroutine magmaf_sormtr
#if defined(PRECISION_z) || defined(PRECISION_c)

     subroutine magmaf_sgeev( jobvl, jobvr, n, a, lda, w, vl, ldvl, vr, ldvr, work, lwork, rwork, info)
       character          :: jobvl
       character          :: jobvr
       integer       :: n
       real    :: a(*)
       integer       :: lda
       real    :: w(*)
       real    :: vl(*)
       integer       :: ldvl
       real    :: vr(*)
       integer       :: ldvr
       real    :: work(*)
       integer       :: lwork
       real:: rwork(*)
       integer       :: info
     end subroutine magmaf_sgeev

     subroutine magmaf_sgesvd( jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, rwork, info)
       character          :: jobu
       character          :: jobvt
       integer       :: m
       integer       :: n
       real    :: a(*)
       integer       :: lda
       real:: s(*)
       real    :: u(*)
       integer       :: ldu
       real    :: vt(*)
       integer       :: ldvt
       real    :: work(*)
       integer       :: lwork
       real:: rwork(*)
       integer       :: info
     end subroutine magmaf_sgesvd

     subroutine magmaf_ssyevd( jobz, uplo, n, a, lda, w, work, lwork, rwork, lrwork, iwork, liwork, info)
       character     :: jobz
       character     :: uplo
       integer       :: n
       real    :: a(*)
       integer       :: lda
       real:: w(*)
       real    :: work(*)
       integer       :: lwork
       real:: rwork(*)
       integer       :: lrwork
       integer       :: iwork(*)
       integer       :: liwork
       integer       :: info
     end subroutine magmaf_ssyevd

     subroutine magmaf_ssygvd( itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, rwork, lrwork, iwork, liwork, info)
       integer       :: itype
       character     :: jobz
       character     :: uplo
       integer       :: n
       real    :: a(*)
       integer       :: lda
       real    :: b(*)
       integer       :: ldb
       real:: w(*)
       real    :: work(*)
       integer       :: lwork
       real:: rwork(*)
       integer       :: lrwork
       integer       :: iwork(*)
       integer       :: liwork
       integer       :: info
     end subroutine magmaf_ssygvd

#else
     subroutine magmaf_sgeev( jobvl, jobvr, n, a, lda, wr, wi, vl, ldvl, vr, ldvr, work, lwork, info)
       character          :: jobvl
       character          :: jobvr
       integer       :: n
       real    :: a(*)
       integer       :: lda
       real    :: wr(*)
       real    :: wi(*)
       real    :: vl(*)
       integer       :: ldvl
       real    :: vr(*)
       integer       :: ldvr
       real    :: work(*)
       integer       :: lwork
       integer       :: info
     end subroutine magmaf_sgeev

     subroutine magmaf_sgesvd( jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, info)
       character          :: jobu
       character          :: jobvt
       integer       :: m
       integer       :: n
       real    :: a(*)
       integer       :: lda
       real:: s(*)
       real    :: u(*)
       integer       :: ldu
       real    :: vt(*)
       integer       :: ldvt
       real    :: work(*)
       integer       :: lwork
       integer       :: info
     end subroutine magmaf_sgesvd

     subroutine magmaf_ssyevd( jobz, uplo, n, a, lda, w, work, lwork, iwork, liwork, info)
       character          :: jobz
       character          :: uplo
       integer       :: n
       real    :: a(*)
       integer       :: lda
       real:: w(*)
       real    :: work(*)
       integer       :: lwork
       integer       :: iwork(*)
       integer       :: liwork
       integer       :: info
     end subroutine magmaf_ssyevd

     subroutine magmaf_ssygvd( itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, iwork, liwork, info)
       integer       :: itype
       character     :: jobz
       character     :: uplo
       integer       :: n
       real    :: a(*)
       integer       :: lda
       real    :: b(*)
       integer       :: ldb
       real:: w(*)
       real    :: work(*)
       integer       :: lwork
       integer       :: iwork(*)
       integer       :: liwork
       integer       :: info
     end subroutine magmaf_ssygvd
#endif

     subroutine magmaf_sgels_gpu(  trans, m, n, nrhs, dA, ldda, dB, lddb, hwork, lwork, info)
       character          :: trans
       integer       :: m
       integer       :: n
       integer       :: nrhs
       magma_devptr_t:: dA
       integer       :: ldda
       magma_devptr_t:: dB
       integer       :: lddb
       real    :: hwork(*)
       integer       :: lwork
       integer       :: info
     end subroutine magmaf_sgels_gpu

     subroutine magmaf_sgeqrf_gpu( m, n, dA, ldda, tau, dT, info)
       integer       :: m
       integer       :: n
       magma_devptr_t:: dA
       integer       :: ldda
       real    :: tau(*)
       magma_devptr_t:: dT
       integer       :: info
     end subroutine magmaf_sgeqrf_gpu

     subroutine magmaf_sgeqrf2_gpu(m, n, dA, ldda, tau, info)
       integer       :: m
       integer       :: n
       magma_devptr_t:: dA
       integer       :: ldda
       real    :: tau(*)
       integer       :: info
     end subroutine magmaf_sgeqrf2_gpu

     subroutine magmaf_sgeqrf3_gpu(m, n, dA, ldda, tau, dT, info)
       integer       :: m
       integer       :: n
       magma_devptr_t:: dA
       integer       :: ldda
       real    :: tau(*)
       magma_devptr_t:: dT
       integer       :: info
     end subroutine magmaf_sgeqrf3_gpu

     subroutine magmaf_sgeqrs_gpu( m, n, nrhs, dA, ldda, tau, dT, dB, lddb, hwork, lhwork, info)
       integer       :: m
       integer       :: n
       integer       :: nrhs
       magma_devptr_t:: dA
       integer       :: ldda
       real    :: tau
       magma_devptr_t:: dT
       magma_devptr_t:: dB
       integer       :: lddb
       real    :: hwork(*)
       integer       :: lhwork
       integer       :: info
     end subroutine magmaf_sgeqrs_gpu

     subroutine magmaf_sgeqrs3_gpu( m, n, nrhs, dA, ldda, tau, dT, dB, lddb, hwork, lhwork, info)
       integer       :: m
       integer       :: n
       integer       :: nrhs
       magma_devptr_t:: dA
       integer       :: ldda
       real    :: tau
       magma_devptr_t:: dT
       magma_devptr_t:: dB
       integer       :: lddb
       real    :: hwork(*)
       integer       :: lhwork
       integer       :: info
     end subroutine magmaf_sgeqrs3_gpu

     subroutine magmaf_sgessm_gpu( storev, m, n, k, ib, ipiv, dL1, lddl1, dL,  lddl, dA,  ldda, info)
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
     end subroutine magmaf_sgessm_gpu

     subroutine magmaf_sgesv_gpu(  n, nrhs, dA, ldda, ipiv, dB, lddb, info)
       integer       :: n
       integer       :: nrhs
       magma_devptr_t:: dA
       integer       :: ldda
       integer       :: ipiv(*)
       magma_devptr_t:: dB
       integer       :: lddb
       integer       :: info
     end subroutine magmaf_sgesv_gpu

     subroutine magmaf_sgetrf_gpu( m, n, dA, ldda, ipiv, info)
       integer       :: m
       integer       :: n
       magma_devptr_t:: dA
       integer       :: ldda
       integer       :: ipiv(*)
       integer       :: info
     end subroutine magmaf_sgetrf_gpu

     subroutine magmaf_sgetrs_gpu( trans, n, nrhs, dA, ldda, ipiv, dB, lddb, info)
       character          :: trans
       integer       :: n
       integer       :: nrhs
       magma_devptr_t:: dA
       integer       :: ldda
       integer       :: ipiv(*)
       magma_devptr_t:: dB
       integer       :: lddb
       integer       :: info
     end subroutine magmaf_sgetrs_gpu

     subroutine magmaf_ssygst_gpu( itype, uplo, n, dA, ldda, dB, lddb, info)
       integer       :: itype
       character     :: uplo
       integer       :: n
       magma_devptr_t:: dA(*)
       integer       :: ldda
       magma_devptr_t:: dB(*)
       integer       :: lddb
       integer       :: info
     end subroutine magmaf_ssygst_gpu

     subroutine magmaf_slabrd_gpu( m, n, nb, a, lda, da, ldda, d, e, tauq, taup, x, ldx, dx, lddx, y, ldy, dy, lddy)
       integer       :: m
       integer       :: n
       integer       :: nb
       real    :: a(*)
       integer       :: lda
       magma_devptr_t:: da
       integer       :: ldda
       real:: d(*)
       real:: e(*)
       real    :: tauq(*)
       real    :: taup(*)
       real    :: x(*)
       integer       :: ldx
       magma_devptr_t:: dx
       integer       :: lddx
       real    :: y(*)
       integer       :: ldy
       magma_devptr_t:: dy
       integer       :: lddy
     end subroutine magmaf_slabrd_gpu

     subroutine magmaf_slarfb_gpu( side, trans, direct, storev, m, n, k, dv, ldv, dt, ldt, dc, ldc, dowrk, ldwork)
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
     end subroutine magmaf_slarfb_gpu

     subroutine magmaf_sposv_gpu(  uplo, n, nrhs, dA, ldda, dB, lddb, info)
       character          :: uplo
       integer       :: n
       integer       :: nrhs
       magma_devptr_t:: dA
       integer       :: ldda
       magma_devptr_t:: dB
       integer       :: lddb
       integer       :: info
     end subroutine magmaf_sposv_gpu

     subroutine magmaf_spotrf_gpu( uplo, n, dA, ldda, info)
       character          :: uplo
       integer       :: n
       magma_devptr_t:: dA
       integer       :: ldda
       integer       :: info
     end subroutine magmaf_spotrf_gpu

     subroutine magmaf_spotrs_gpu( uplo,  n, nrhs, dA, ldda, dB, lddb, info)
       character          :: uplo
       integer       :: n
       integer       :: nrhs
       magma_devptr_t:: dA
       integer       :: ldda
       magma_devptr_t:: dB
       integer       :: lddb
       integer       :: info
     end subroutine magmaf_spotrs_gpu

     subroutine magmaf_sssssm_gpu( storev, m1, n1, m2, n2, k, ib, dA1, ldda1, dA2, ldda2, dL1, lddl1, dL2, lddl2, IPIV, info)
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
     end subroutine magmaf_sssssm_gpu

     subroutine magmaf_sorgqr_gpu( m, n, k, da, ldda, tau, dwork, nb, info)
       integer       :: m
       integer       :: n
       integer       :: k
       magma_devptr_t:: da
       integer       :: ldda
       real    :: tau(*)
       magma_devptr_t:: dwork
       integer       :: nb
       integer       :: info
     end subroutine magmaf_sorgqr_gpu

     subroutine magmaf_sormqr_gpu( side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, td, nb, info)
       character          :: side
       character          :: trans
       integer       :: m
       integer       :: n
       integer       :: k
       magma_devptr_t:: a
       integer       :: lda
       real    :: tau(*)
       magma_devptr_t:: c
       integer       :: ldc
       magma_devptr_t:: work
       integer       :: lwork
       magma_devptr_t:: td
       integer       :: nb
       integer       :: info
     end subroutine magmaf_sormqr_gpu

  end interface

contains
  
  subroutine magmaf_soff1d( ptrNew, ptrOld, inc, i)
    magma_devptr_t :: ptrNew
    magma_devptr_t :: ptrOld
    integer        :: inc, i
    
    ptrNew = ptrOld + (i-1) * inc * sizeof_real
    
  end subroutine magmaf_soff1d
  
  subroutine magmaf_soff2d( ptrNew, ptrOld, lda, i, j)
    magma_devptr_t :: ptrNew
    magma_devptr_t :: ptrOld
    integer        :: lda, i, j
    
    ptrNew = ptrOld + ((j-1) * lda + (i-1)) * sizeof_real
    
  end subroutine magmaf_soff2d
  
end module magma_sfortran
