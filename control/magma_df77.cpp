/*
    -- MAGMA (version 1.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2011

       @generated d Sun Nov 13 20:48:00 2011

*/

#include <ctype.h>
#include <stdio.h>
#include <string.h>
#include <stddef.h>
#include <stdlib.h>
#if defined(__GNUC__)
#include <stdint.h>
#endif /* __GNUC__ */
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cublas.h>
#include "magma.h"

/* 
 * typedef comming from fortran.h file provided in $CUDADIR/src directory
 * it will probably change with future release of cublas when they will use 64bits address
 */
typedef size_t devptr_t;

#define PRECISION_d

#ifdef PGI_FORTRAN
#define DEVPTR(__ptr) ((double*)(__ptr))
#else
#define DEVPTR(__ptr) ((double*)(uintptr_t)(*(__ptr)))
#endif


#ifndef MAGMA_FORTRAN_NAME
#if defined(ADD_)
#define MAGMA_FORTRAN_NAME(lcname, UCNAME)  magmaf_##lcname##_
#elif defined(NOCHANGE)
#define MAGMA_FORTRAN_NAME(lcname, UCNAME)  magmaf_##lcname
#elif defined(UPCASE)
#define MAGMA_FORTRAN_NAME(lcname, UCNAME)  MAGMAF_##UCNAME
#endif
#endif

#ifndef MAGMA_GPU_FORTRAN_NAME
#if defined(ADD_)
#define MAGMA_GPU_FORTRAN_NAME(lcname, UCNAME)  magmaf_##lcname##_gpu_
#elif defined(NOCHANGE)
#define MAGMA_GPU_FORTRAN_NAME(lcname, UCNAME)  magmaf_##lcname##_gpu
#elif defined(UPCASE)
#define MAGMA_GPU_FORTRAN_NAME(lcname, UCNAME)  MAGMAF_##UCNAME##_GPU
#endif
#endif

/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA function definitions / Data on CPU
*/
#define MAGMAF_DGEBRD  MAGMA_FORTRAN_NAME(dgebrd,  DGEBRD ) 
#define MAGMAF_DGEHRD2 MAGMA_FORTRAN_NAME(dgehrd2, DGEHRD2)
#define MAGMAF_DGEHRD  MAGMA_FORTRAN_NAME(dgehrd,  DGEHRD )
#define MAGMAF_DGELQF  MAGMA_FORTRAN_NAME(dgelqf,  DGELQF )
#define MAGMAF_ZGEQLF  MAGMA_FORTRAN_NAME(dgeqlf,  ZGEQLF )
#define MAGMAF_DGEQRF  MAGMA_FORTRAN_NAME(dgeqrf,  DGEQRF )
#define MAGMAF_DGESV   MAGMA_FORTRAN_NAME(dgesv,   DGESV  )
#define MAGMAF_DGETRF  MAGMA_FORTRAN_NAME(dgetrf,  DGETRF )
#define MAGMAF_DLATRD  MAGMA_FORTRAN_NAME(dlatrd,  DLATRD )
#define MAGMAF_ZLAHR2  MAGMA_FORTRAN_NAME(dlahr2,  ZLAHR2 )
#define MAGMAF_ZLAHRU  MAGMA_FORTRAN_NAME(dlahru,  ZLAHRU )
#define MAGMAF_DPOSV   MAGMA_FORTRAN_NAME(dposv,   DPOSV  )
#define MAGMAF_DPOTRF  MAGMA_FORTRAN_NAME(dpotrf,  DPOTRF )
#define MAGMAF_DSYTRD  MAGMA_FORTRAN_NAME(dsytrd,  DSYTRD )
#define MAGMAF_DORGQR  MAGMA_FORTRAN_NAME(dorgqr,  DORGQR )
#define MAGMAF_DORMQR  MAGMA_FORTRAN_NAME(dormqr,  DORMQR )
#define MAGMAF_DORMTR  MAGMA_FORTRAN_NAME(dormtr,  DORMTR )
#define MAGMAF_DORGHR  MAGMA_FORTRAN_NAME(dorghr,  DORGHR )
#define MAGMAF_DGEEV   MAGMA_FORTRAN_NAME(dgeev,   DGEEV  )
#define MAGMAF_DGESVD  MAGMA_FORTRAN_NAME(dgesvd,  DGESVD )
#define MAGMAF_DSYEVD  MAGMA_FORTRAN_NAME(dsyevd,  DSYEVD )
#define MAGMAF_DSYGVD  MAGMA_FORTRAN_NAME(dsygvd,  DSYGVD )

/* //////////////////////////////////////////////////////////////////////////// 
 -- MAGMA function definitions / Data on GPU
*/
#define MAGMAF_DGELS_GPU   MAGMA_GPU_FORTRAN_NAME(dgels,   DGELS  )
#define MAGMAF_DGEQRF_GPU  MAGMA_GPU_FORTRAN_NAME(dgeqrf,  DGEQRF ) 
#define MAGMAF_DGEQRF2_GPU MAGMA_GPU_FORTRAN_NAME(dgeqrf2, DGEQRF2)
#define MAGMAF_DGEQRF3_GPU MAGMA_GPU_FORTRAN_NAME(dgeqrf3, DGEQRF3)
#define MAGMAF_DGEQRS_GPU  MAGMA_GPU_FORTRAN_NAME(dgeqrs,  DGEQRS ) 
#define MAGMAF_DGEQRS3_GPU MAGMA_GPU_FORTRAN_NAME(dgeqrs3, DGEQRS3) 
#define MAGMAF_ZGESSM_GPU  MAGMA_GPU_FORTRAN_NAME(dgessm,  ZGESSM ) 
#define MAGMAF_DGESV_GPU   MAGMA_GPU_FORTRAN_NAME(dgesv,   DGESV  )  
#define MAGMAF_ZGETRL_GPU  MAGMA_GPU_FORTRAN_NAME(dgetrl,  ZGETRL ) 
#define MAGMAF_DGETRF_GPU  MAGMA_GPU_FORTRAN_NAME(dgetrf,  DGETRF ) 
#define MAGMAF_DGETRS_GPU  MAGMA_GPU_FORTRAN_NAME(dgetrs,  DGETRS ) 
#define MAGMAF_DLABRD_GPU  MAGMA_GPU_FORTRAN_NAME(dlabrd,  DLABRD ) 
#define MAGMAF_DLARFB_GPU  MAGMA_GPU_FORTRAN_NAME(dlarfb,  DLARFB ) 
#define MAGMAF_DPOSV_GPU   MAGMA_GPU_FORTRAN_NAME(dposv,   DPOSV  )  
#define MAGMAF_DPOTRF_GPU  MAGMA_GPU_FORTRAN_NAME(dpotrf,  DPOTRF ) 
#define MAGMAF_DPOTRS_GPU  MAGMA_GPU_FORTRAN_NAME(dpotrs,  DPOTRS ) 
#define MAGMAF_ZSSSSM_GPU  MAGMA_GPU_FORTRAN_NAME(dssssm,  ZSSSSM ) 
#define MAGMAF_ZTSTRF_GPU  MAGMA_GPU_FORTRAN_NAME(dtstrf,  ZTSTRF ) 
#define MAGMAF_DORGQR_GPU  MAGMA_GPU_FORTRAN_NAME(dorgqr,  DORGQR ) 
#define MAGMAF_DORMQR_GPU  MAGMA_GPU_FORTRAN_NAME(dormqr,  DORMQR ) 

#ifdef __cplusplus
extern "C" {
#endif

/***************************************************************************//**
 *  FORTRAN API - math functions (simple interface)
 **/

/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA function definitions / Data on CPU
*/
void MAGMAF_DGEBRD( magma_int_t *m, magma_int_t *n, double *A, 
                    magma_int_t *lda, double *d, double *e,
                    double *tauq, double *taup, 
                    double *work, magma_int_t *lwork, magma_int_t *info)
{
    magma_dgebrd( *m, *n, A, 
                  *lda, d, e,
                  tauq, taup, 
                  work, *lwork, info);
}
    
void MAGMAF_DGEHRD2(magma_int_t *n, magma_int_t *ilo, magma_int_t *ihi,
                    double *A, magma_int_t *lda, double *tau, 
                    double *work, magma_int_t *lwork, magma_int_t *info)
{
    magma_dgehrd2(*n, *ilo, *ihi,
                  A, *lda, tau, 
                  work, lwork, info);
}
    
void MAGMAF_DGEHRD( magma_int_t *n, magma_int_t *ilo, magma_int_t *ihi,
                    double *A, magma_int_t *lda, double *tau,
                    double *work, magma_int_t *lwork,
                    double *d_T, magma_int_t *info)
{
  magma_dgehrd( *n, *ilo, *ihi,
                A, *lda, tau,
                work, *lwork,
                d_T, info);
}

void MAGMAF_DGELQF( magma_int_t *m, magma_int_t *n, 
                    double *A,    magma_int_t *lda,   double *tau, 
                    double *work, magma_int_t *lwork, magma_int_t *info)
{
    magma_dgelqf( *m, *n, 
                  A,    *lda,   tau, 
                  work, *lwork, info);
}

void MAGMAF_ZGEQLF( magma_int_t *m, magma_int_t *n, 
                    double *A,    magma_int_t *lda,   double *tau, 
                    double *work, magma_int_t *lwork, magma_int_t *info)
{
    magma_dgeqlf( *m, *n, 
                  A,    *lda,   tau, 
                  work, *lwork, info);
}

void MAGMAF_DGEQRF( magma_int_t *m, magma_int_t *n, double *A, 
                    magma_int_t *lda, double *tau, double *work, 
                    magma_int_t *lwork, magma_int_t *info)
{
    magma_dgeqrf( *m, *n, A, 
                  *lda, tau, work, 
                  *lwork, info);
}

void MAGMAF_DGESV ( magma_int_t *n, magma_int_t *nrhs,
                    double *A, magma_int_t *lda, magma_int_t *ipiv,
                    double *B, magma_int_t *ldb, magma_int_t *info)
{
    magma_dgesv(  *n, *nrhs,
                  A, *lda, ipiv,
                  B, *ldb,
                  info);
}
    
void MAGMAF_DGETRF( magma_int_t *m, magma_int_t *n, double *A, 
                    magma_int_t *lda, magma_int_t *ipiv, 
                    magma_int_t *info)
{
    magma_dgetrf( *m, *n, A, 
                  *lda, ipiv, 
                  info);
}

// void MAGMAF_DLATRD( char *uplo, magma_int_t *n, magma_int_t *nb, double *a, 
//                     magma_int_t *lda, double *e, double *tau, 
//                     double *w, magma_int_t *ldw,
//                     double *da, magma_int_t *ldda, 
//                     double *dw, magma_int_t *lddw)
// {
//     magma_dlatrd( uplo[0], *n, *nb, a, 
//                   *lda, e, tau, 
//                   w, *ldw,
//                   da, *ldda, 
//                   dw, *lddw);
// }

  /* This has nothing to do here, it should be a GPU function */
// void MAGMAF_ZLAHR2( magma_int_t *m, magma_int_t *n, magma_int_t *nb, 
//                     double *da, double *dv, double *a, 
//                     magma_int_t *lda, double *tau, double *t, 
//                     magma_int_t *ldt, double *y, magma_int_t *ldy)
// {
//     magma_dlahr2( *m, *n, *nb, 
//                   da, dv, a, 
//                   *lda, tau, t, 
//                   *ldt, y, *ldy);
// }

// void MAGMAF_ZLAHRU( magma_int_t *m, magma_int_t *n, magma_int_t *nb, 
//                     double *a, magma_int_t *lda, 
//                     double *da, double *y, 
//                     double *v, double *t, 
//                     double *dwork)
// {
//     magma_dlahru( *m, *n, *nb, 
//                   a, *lda, 
//                   da, y, 
//                   v, t, 
//                   dwork);
// }

void MAGMAF_DPOSV(  char *uplo, magma_int_t *n, magma_int_t *nrhs,
                    double *A, magma_int_t *lda,
                    double *B, magma_int_t *ldb, magma_int_t *info)
{
    magma_dposv(  uplo[0], *n, *nrhs,
                  A, *lda,
                  B, *ldb, info);
}

void MAGMAF_DPOTRF( char *uplo, magma_int_t *n, double *A, 
                    magma_int_t *lda, magma_int_t *info)
{
    magma_dpotrf( uplo[0], *n, A, 
                  *lda, info);
}

void MAGMAF_DSYTRD( char *uplo, magma_int_t *n, double *A, 
                    magma_int_t *lda, double *d, double *e, 
                    double *tau, double *work, magma_int_t *lwork, 
                    magma_int_t *info)
{
    magma_dsytrd( uplo[0], *n, A, 
                  *lda, d, e, 
                  tau, work, *lwork, 
                  info);
}

// void MAGMAF_DORGQR( magma_int_t *m, magma_int_t *n, magma_int_t *k,
//                     double *a, magma_int_t *lda,
//                     double *tau, double *dwork,
//                     magma_int_t *nb, magma_int_t *info )
// {
//     magma_dorgqr( *m, *n, *k,
//                   a, *lda,
//                   tau, dwork,
//                   *nb, info );
// }

void MAGMAF_DORMQR( char *side, char *trans, 
                    magma_int_t *m, magma_int_t *n, magma_int_t *k, 
                    double *a, magma_int_t *lda, double *tau, 
                    double *c, magma_int_t *ldc, 
                    double *work, magma_int_t *lwork, magma_int_t *info)
{
    magma_dormqr( side[0], trans[0], 
                  *m, *n, *k, 
                  a, *lda, tau, 
                  c, *ldc, 
                  work, *lwork, info);
}

void MAGMAF_DORMTR( char *side, char *uplo, char *trans,
                    magma_int_t *m, magma_int_t *n,
                    double *a,    magma_int_t *lda,
                    double *tau,
                    double *c,    magma_int_t *ldc,
                    double *work, magma_int_t *lwork,
                    magma_int_t *info)
{
    magma_dormtr( side[0], uplo[0], trans[0],
                  *m, *n,
                  a,    *lda,
                  tau,
                  c,    *ldc,
                  work, *lwork,
                  info);
}

// void MAGMAF_DORGHR( magma_int_t *n, magma_int_t *ilo, magma_int_t *ihi,
//                     double *a, magma_int_t *lda,
//                     double *tau,
//                     double *dT, magma_int_t *nb,
//                     magma_int_t *info)
// {
//     magma_dorghr( *n, *ilo, *ihi,
//                   a, *lda,
//                   tau,
//                   dT, *nb,
//                   info);
// }

#if defined(PRECISION_z) || defined(PRECISION_c)
void MAGMAF_DGEEV( char *jobvl, char *jobvr, magma_int_t *n,
                   double *a, magma_int_t *lda,
                   double *w,
                   double *vl, magma_int_t *ldvl,
                   double *vr, magma_int_t *ldvr,
                   double *work, magma_int_t *lwork,
                   double *rwork, magma_int_t *info)
{
    magma_dgeev( jobvl[0], jobvr[0], *n,
                 a, *lda,
                 w,
                 vl, *ldvl,
                 vr, *ldvr,
                 work, *lwork,
                 rwork, info);
}

void MAGMAF_DGESVD( char *jobu, char *jobvt, magma_int_t *m, magma_int_t *n,
                    double *a,    magma_int_t *lda, double *s, 
                    double *u,    magma_int_t *ldu, 
                    double *vt,   magma_int_t *ldvt,
                    double *work, magma_int_t *lwork,
                    double *rwork, magma_int_t *info )
{
    magma_dgesvd( jobu[0], jobvt[0], *m, *n,
                  a,    *lda, s, 
                  u,    *ldu, 
                  vt,   *ldvt,
                  work, *lwork,
                  rwork, info );
}
    
void MAGMAF_DSYEVD( char *jobz, char *uplo, magma_int_t *n,
                    double *a,     magma_int_t *lda, double *w,
                    double *work,  magma_int_t *lwork,
                    double          *rwork, magma_int_t *lrwork,
                    magma_int_t *iwork, magma_int_t *liwork, magma_int_t *info)
{
    magma_dsyevd( jobz[0], uplo[0], *n,
                  a, *lda, w,
                  work, *lwork,
                  rwork, *lrwork,
                  iwork, *liwork, info);
}

void MAGMAF_DSYGVD(magma_int_t *itype, char *jobz, char *uplo, magma_int_t *n,
                   double *a, magma_int_t *lda, 
                   double *b, magma_int_t *ldb,
                   double *w, double *work, magma_int_t *lwork,
                   double *rwork, magma_int_t *lrwork,
                   magma_int_t *iwork, magma_int_t *liwork, magma_int_t *info)
{
  magma_dsygvd( *itype, jobz[0], uplo[0], *n,
                a, *lda, b, *ldb,
                w, work, *lwork,
                rwork, *lrwork,
                iwork, *liwork, info);
}
    
#else
void MAGMAF_DGEEV( char *jobvl, char *jobvr, magma_int_t *n,
                   double *a,    magma_int_t *lda,
                   double *wr, double *wi,
                   double *vl,   magma_int_t *ldvl,
                   double *vr,   magma_int_t *ldvr,
                   double *work, magma_int_t *lwork,
                   magma_int_t *info)
{
    magma_dgeev( jobvl[0], jobvr[0], *n,
                 a,    *lda,
                 wr, wi,
                 vl,   *ldvl,
                 vr,   *ldvr,
                 work, *lwork,
                 info);
}

void MAGMAF_DGESVD( char *jobu, char *jobvt, magma_int_t *m, magma_int_t *n,
                    double *a,    magma_int_t *lda, double *s,
                    double *u,    magma_int_t *ldu, 
                    double *vt,   magma_int_t *ldvt,
                    double *work, magma_int_t *lwork,
                    magma_int_t *info )
{
    magma_dgesvd( jobu[0], jobvt[0], *m, *n,
                  a,    *lda, s, 
                  u,    *ldu, 
                  vt,   *ldvt,
                  work, *lwork,
                  info );
}

void MAGMAF_DSYEVD( char *jobz, char *uplo, magma_int_t *n,
                    double *a, magma_int_t *lda, double *w,
                    double *work, magma_int_t *lwork,
                    magma_int_t *iwork, magma_int_t *liwork, magma_int_t *info)
{
    magma_dsyevd( jobz[0], uplo[0], *n,
                  a, *lda, w,
                  work, *lwork,
                  iwork, *liwork, info);
}

void MAGMAF_DSYGVD(magma_int_t *itype, char *jobz, char *uplo, magma_int_t *n,
                   double *a, magma_int_t *lda,
                   double *b, magma_int_t *ldb,
                   double *w, double *work, magma_int_t *lwork,
                   magma_int_t *iwork, magma_int_t *liwork, magma_int_t *info)
{
  magma_dsygvd( *itype, jobz[0], uplo[0], *n,
                a, *lda, b, *ldb,
                w, work, *lwork,
                iwork, *liwork, info);
}


#endif

/* //////////////////////////////////////////////////////////////////////////// 
 -- MAGMA function definitions / Data on GPU
*/
void MAGMAF_DGELS_GPU(  char *trans, magma_int_t *m, magma_int_t *n, magma_int_t *nrhs,
                        devptr_t *dA,    magma_int_t *ldda, 
                        devptr_t *dB,    magma_int_t *lddb, 
                        double *hwork, magma_int_t *lwork, 
                        magma_int_t *info)
{
    magma_dgels_gpu(  trans[0], *m, *n, *nrhs, 
                      DEVPTR(dA),    *ldda,  
                      DEVPTR(dB),    *lddb,  
                      hwork, *lwork,  info);
}

void MAGMAF_DGEQRF_GPU( magma_int_t *m, magma_int_t *n, 
                        devptr_t *dA,  magma_int_t *ldda, 
                        double *tau, devptr_t *dT, 
                        magma_int_t *info)
{
    magma_dgeqrf_gpu( *m, *n,  
                      DEVPTR(dA),  *ldda,  
                      tau, 
                      DEVPTR(dT),  info);
}

void MAGMAF_DGEQRF2_GPU(magma_int_t *m, magma_int_t *n, 
                        devptr_t *dA,  magma_int_t *ldda, 
                        double *tau, magma_int_t *info)
{
    magma_dgeqrf2_gpu(*m, *n,  
                      DEVPTR(dA),  *ldda,  
                      tau, info); 
}

void MAGMAF_DGEQRF3_GPU(magma_int_t *m, magma_int_t *n, 
                        devptr_t *dA,  magma_int_t *ldda, 
                        double *tau, devptr_t *dT,
                        magma_int_t *info)
{
    magma_dgeqrf3_gpu(*m, *n,  
                      DEVPTR(dA),  *ldda,  
                      tau, DEVPTR(dT), info); 
}

void MAGMAF_DGEQRS_GPU( magma_int_t *m, magma_int_t *n, magma_int_t *nrhs, 
                        devptr_t *dA,     magma_int_t *ldda, 
                        double *tau,   devptr_t *dT,
                        devptr_t *dB,    magma_int_t *lddb,
                        double *hwork, magma_int_t *lhwork, 
                        magma_int_t *info)
{
    magma_dgeqrs_gpu( *m, *n, *nrhs,  
                      DEVPTR(dA),     *ldda,  
                      tau,
                      DEVPTR(dT), 
                      DEVPTR(dB),    *lddb, 
                      hwork, *lhwork,  info);
}

void MAGMAF_DGEQRS3_GPU(magma_int_t *m, magma_int_t *n, magma_int_t *nrhs, 
                        devptr_t *dA,     magma_int_t *ldda, 
                        double *tau,   devptr_t *dT,
                        devptr_t *dB,    magma_int_t *lddb,
                        double *hwork, magma_int_t *lhwork, 
                        magma_int_t *info)
{
    magma_dgeqrs3_gpu(*m, *n, *nrhs,  
                      DEVPTR(dA),     *ldda,  
                      tau,
                      DEVPTR(dT), 
                      DEVPTR(dB),    *lddb, 
                      hwork, *lhwork,  info);
}

void MAGMAF_ZGESSM_GPU( char *storev, magma_int_t *m, magma_int_t *n, magma_int_t *k, magma_int_t *ib, 
                        magma_int_t *ipiv, 
                        devptr_t *dL1, magma_int_t *lddl1, 
                        devptr_t *dL,  magma_int_t *lddl, 
                        devptr_t *dA,  magma_int_t *ldda, 
                        magma_int_t *info)
{
    magma_dgessm_gpu( storev[0], *m, *n, *k, *ib, ipiv,  
                      DEVPTR(dL1), *lddl1,  
                      DEVPTR(dL),  *lddl,  
                      DEVPTR(dA),  *ldda,  info);
}

void MAGMAF_DGESV_GPU(  magma_int_t *n, magma_int_t *nrhs, 
                        devptr_t *dA, magma_int_t *ldda, magma_int_t *ipiv, 
                        devptr_t *dB, magma_int_t *lddb, magma_int_t *info)
{
    magma_dgesv_gpu(  *n, *nrhs,  
                      DEVPTR(dA), *ldda, ipiv,  
                      DEVPTR(dB), *lddb, info);
}

void MAGMAF_DGETRF_GPU( magma_int_t *m, magma_int_t *n, 
                        devptr_t *dA, magma_int_t *ldda, 
                        magma_int_t *ipiv, magma_int_t *info)
{
    magma_dgetrf_gpu( *m, *n,  
                      DEVPTR(dA), *ldda, ipiv, info);
}

void MAGMAF_DGETRS_GPU( char *trans, magma_int_t *n, magma_int_t *nrhs, 
                        devptr_t *dA, magma_int_t *ldda, magma_int_t *ipiv, 
                        devptr_t *dB, magma_int_t *lddb, magma_int_t *info)
{
    magma_dgetrs_gpu( trans[0], *n, *nrhs,  
                      DEVPTR(dA), *ldda, ipiv,  
                      DEVPTR(dB), *lddb, info);
}

void MAGMAF_DLABRD_GPU( magma_int_t *m, magma_int_t *n, magma_int_t *nb, 
                        double *a, magma_int_t *lda, devptr_t *da, magma_int_t *ldda,
                        double *d, double *e, double *tauq, double *taup,  
                        double *x, magma_int_t *ldx, devptr_t *dx, magma_int_t *lddx, 
                        double *y, magma_int_t *ldy, devptr_t *dy, magma_int_t *lddy)
{
    magma_dlabrd_gpu( *m, *n, *nb,  
                      a, *lda, DEVPTR(da), *ldda, 
                      d, e, tauq, taup,   
                      x, *ldx, DEVPTR(dx), *lddx,  
                      y, *ldy, DEVPTR(dy), *lddy);
}

void MAGMAF_DLARFB_GPU( char *side, char *trans, char *direct, char *storev, 
                        magma_int_t *m, magma_int_t *n, magma_int_t *k,
                        devptr_t *dv, magma_int_t *ldv, devptr_t *dt,    magma_int_t *ldt, 
                        devptr_t *dc, magma_int_t *ldc, devptr_t *dowrk, magma_int_t *ldwork )
{
    magma_dlarfb_gpu( side[0], trans[0], direct[0], storev[0],  *m, *n, *k, 
                      DEVPTR(dv), *ldv, DEVPTR(dt),    *ldt,  
                      DEVPTR(dc), *ldc, DEVPTR(dowrk), *ldwork);
}

void MAGMAF_DPOSV_GPU(  char *uplo, magma_int_t *n, magma_int_t *nrhs, 
                        devptr_t *dA, magma_int_t *ldda, 
                        devptr_t *dB, magma_int_t *lddb, magma_int_t *info)
{
    magma_dposv_gpu(  uplo[0], *n, *nrhs,  
                      DEVPTR(dA), *ldda,  
                      DEVPTR(dB), *lddb, info);
}

void MAGMAF_DPOTRF_GPU( char *uplo,  magma_int_t *n, 
                        devptr_t *dA, magma_int_t *ldda, magma_int_t *info)
{
    magma_dpotrf_gpu( uplo[0],  *n,  
                      DEVPTR(dA), *ldda, info); }

void MAGMAF_DPOTRS_GPU( char *uplo,  magma_int_t *n, magma_int_t *nrhs, 
                        devptr_t *dA, magma_int_t *ldda, 
                        devptr_t *dB, magma_int_t *lddb, magma_int_t *info)
{
    magma_dpotrs_gpu( uplo[0],  *n, *nrhs,  
                      DEVPTR(dA), *ldda,  
                      DEVPTR(dB), *lddb, info);
}

void MAGMAF_ZSSSSM_GPU( char *storev, magma_int_t *m1, magma_int_t *n1, 
                        magma_int_t *m2, magma_int_t *n2, magma_int_t *k, magma_int_t *ib, 
                        devptr_t *dA1, magma_int_t *ldda1, 
                        devptr_t *dA2, magma_int_t *ldda2, 
                        devptr_t *dL1, magma_int_t *lddl1, 
                        devptr_t *dL2, magma_int_t *lddl2,
                        magma_int_t *IPIV, magma_int_t *info)
{
    magma_dssssm_gpu( storev[0], *m1, *n1,  *m2, *n2, *k, *ib,  
                      DEVPTR(dA1), *ldda1,  
                      DEVPTR(dA2), *ldda2,  
                      DEVPTR(dL1), *lddl1,  
                      DEVPTR(dL2), *lddl2,
                      IPIV, info);
}

void MAGMAF_DORGQR_GPU( magma_int_t *m, magma_int_t *n, magma_int_t *k, 
                        devptr_t *da, magma_int_t *ldda, 
                        double *tau, devptr_t *dwork, 
                        magma_int_t *nb, magma_int_t *info )
{
    magma_dorgqr_gpu( *m, *n, *k,  
                      DEVPTR(da), *ldda, tau, 
                      DEVPTR(dwork), *nb, info );
}

void MAGMAF_DORMQR_GPU( char *side, char *trans, 
                        magma_int_t *m, magma_int_t *n, magma_int_t *k,
                        devptr_t *a,    magma_int_t *lda, double *tau, 
                        devptr_t *c,    magma_int_t *ldc,
                        devptr_t *work, magma_int_t *lwork, 
                        devptr_t *td,   magma_int_t *nb, magma_int_t *info)
{
    magma_dormqr_gpu( side[0], trans[0], *m, *n, *k, 
                      DEVPTR(a),    *lda, tau,  
                      DEVPTR(c),    *ldc, 
                      DEVPTR(work), *lwork,  
                      DEVPTR(td),   *nb, info);
}

#ifdef __cplusplus
}
#endif
