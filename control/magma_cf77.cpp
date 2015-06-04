/*
    -- MAGMA (version 1.2.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       May 2012

       @generated c Tue May 15 18:17:09 2012

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

#define PRECISION_c

#ifdef PGI_FORTRAN
#define DEVPTR(__ptr) ((cuFloatComplex*)(__ptr))
#else
#define DEVPTR(__ptr) ((cuFloatComplex*)(uintptr_t)(*(__ptr)))
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
#define MAGMAF_CGEBRD  MAGMA_FORTRAN_NAME(cgebrd,  CGEBRD ) 
#define MAGMAF_CGEHRD2 MAGMA_FORTRAN_NAME(cgehrd2, CGEHRD2)
#define MAGMAF_CGEHRD  MAGMA_FORTRAN_NAME(cgehrd,  CGEHRD )
#define MAGMAF_CGELQF  MAGMA_FORTRAN_NAME(cgelqf,  CGELQF )
#define MAGMAF_ZGEQLF  MAGMA_FORTRAN_NAME(cgeqlf,  ZGEQLF )
#define MAGMAF_CGEQRF  MAGMA_FORTRAN_NAME(cgeqrf,  CGEQRF )
#define MAGMAF_CGESV   MAGMA_FORTRAN_NAME(cgesv,   CGESV  )
#define MAGMAF_CGETRF  MAGMA_FORTRAN_NAME(cgetrf,  CGETRF )
#define MAGMAF_CLATRD  MAGMA_FORTRAN_NAME(clatrd,  CLATRD )
#define MAGMAF_ZLAHR2  MAGMA_FORTRAN_NAME(clahr2,  ZLAHR2 )
#define MAGMAF_ZLAHRU  MAGMA_FORTRAN_NAME(clahru,  ZLAHRU )
#define MAGMAF_CPOSV   MAGMA_FORTRAN_NAME(cposv,   CPOSV  )
#define MAGMAF_CPOTRF  MAGMA_FORTRAN_NAME(cpotrf,  CPOTRF )
#define MAGMAF_CHETRD  MAGMA_FORTRAN_NAME(chetrd,  CHETRD )
#define MAGMAF_CUNGQR  MAGMA_FORTRAN_NAME(cungqr,  CUNGQR )
#define MAGMAF_CUNMQR  MAGMA_FORTRAN_NAME(cunmqr,  CUNMQR )
#define MAGMAF_CUNMTR  MAGMA_FORTRAN_NAME(cunmtr,  CUNMTR )
#define MAGMAF_CUNGHR  MAGMA_FORTRAN_NAME(cunghr,  CUNGHR )
#define MAGMAF_CGEEV   MAGMA_FORTRAN_NAME(cgeev,   CGEEV  )
#define MAGMAF_CGESVD  MAGMA_FORTRAN_NAME(cgesvd,  CGESVD )
#define MAGMAF_CHEEVD  MAGMA_FORTRAN_NAME(cheevd,  CHEEVD )
#define MAGMAF_CHEGVD  MAGMA_FORTRAN_NAME(chegvd,  CHEGVD )

/* //////////////////////////////////////////////////////////////////////////// 
 -- MAGMA function definitions / Data on GPU
*/
#define MAGMAF_CGELS_GPU   MAGMA_GPU_FORTRAN_NAME(cgels,   CGELS  )
#define MAGMAF_CGEQRF_GPU  MAGMA_GPU_FORTRAN_NAME(cgeqrf,  CGEQRF ) 
#define MAGMAF_CGEQRF2_GPU MAGMA_GPU_FORTRAN_NAME(cgeqrf2, CGEQRF2)
#define MAGMAF_CGEQRF3_GPU MAGMA_GPU_FORTRAN_NAME(cgeqrf3, CGEQRF3)
#define MAGMAF_CGEQRS_GPU  MAGMA_GPU_FORTRAN_NAME(cgeqrs,  CGEQRS ) 
#define MAGMAF_CGEQRS3_GPU MAGMA_GPU_FORTRAN_NAME(cgeqrs3, CGEQRS3) 
#define MAGMAF_ZGESSM_GPU  MAGMA_GPU_FORTRAN_NAME(cgessm,  ZGESSM ) 
#define MAGMAF_CGESV_GPU   MAGMA_GPU_FORTRAN_NAME(cgesv,   CGESV  )  
#define MAGMAF_ZGETRL_GPU  MAGMA_GPU_FORTRAN_NAME(cgetrl,  ZGETRL ) 
#define MAGMAF_CGETRF_GPU  MAGMA_GPU_FORTRAN_NAME(cgetrf,  CGETRF ) 
#define MAGMAF_CGETRS_GPU  MAGMA_GPU_FORTRAN_NAME(cgetrs,  CGETRS ) 
#define MAGMAF_CLABRD_GPU  MAGMA_GPU_FORTRAN_NAME(clabrd,  CLABRD ) 
#define MAGMAF_CLARFB_GPU  MAGMA_GPU_FORTRAN_NAME(clarfb,  CLARFB ) 
#define MAGMAF_CPOSV_GPU   MAGMA_GPU_FORTRAN_NAME(cposv,   CPOSV  )  
#define MAGMAF_CPOTRF_GPU  MAGMA_GPU_FORTRAN_NAME(cpotrf,  CPOTRF ) 
#define MAGMAF_CPOTRS_GPU  MAGMA_GPU_FORTRAN_NAME(cpotrs,  CPOTRS ) 
#define MAGMAF_ZSSSSM_GPU  MAGMA_GPU_FORTRAN_NAME(cssssm,  ZSSSSM ) 
#define MAGMAF_ZTSTRF_GPU  MAGMA_GPU_FORTRAN_NAME(ctstrf,  ZTSTRF ) 
#define MAGMAF_CUNGQR_GPU  MAGMA_GPU_FORTRAN_NAME(cungqr,  CUNGQR ) 
#define MAGMAF_CUNMQR_GPU  MAGMA_GPU_FORTRAN_NAME(cunmqr,  CUNMQR ) 

#ifdef __cplusplus
extern "C" {
#endif

/***************************************************************************//**
 *  FORTRAN API - math functions (simple interface)
 **/

/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA function definitions / Data on CPU
*/
void MAGMAF_CGEBRD( magma_int_t *m, magma_int_t *n, cuFloatComplex *A, 
                    magma_int_t *lda, float *d, float *e,
                    cuFloatComplex *tauq, cuFloatComplex *taup, 
                    cuFloatComplex *work, magma_int_t *lwork, magma_int_t *info)
{
    magma_cgebrd( *m, *n, A, 
                  *lda, d, e,
                  tauq, taup, 
                  work, *lwork, info);
}
    
void MAGMAF_CGEHRD2(magma_int_t *n, magma_int_t *ilo, magma_int_t *ihi,
                    cuFloatComplex *A, magma_int_t *lda, cuFloatComplex *tau, 
                    cuFloatComplex *work, magma_int_t *lwork, magma_int_t *info)
{
    magma_cgehrd2(*n, *ilo, *ihi,
                  A, *lda, tau, 
                  work, lwork, info);
}
    
void MAGMAF_CGEHRD( magma_int_t *n, magma_int_t *ilo, magma_int_t *ihi,
                    cuFloatComplex *A, magma_int_t *lda, cuFloatComplex *tau,
                    cuFloatComplex *work, magma_int_t *lwork,
                    cuFloatComplex *d_T, magma_int_t *info)
{
  magma_cgehrd( *n, *ilo, *ihi,
                A, *lda, tau,
                work, *lwork,
                d_T, info);
}

void MAGMAF_CGELQF( magma_int_t *m, magma_int_t *n, 
                    cuFloatComplex *A,    magma_int_t *lda,   cuFloatComplex *tau, 
                    cuFloatComplex *work, magma_int_t *lwork, magma_int_t *info)
{
    magma_cgelqf( *m, *n, 
                  A,    *lda,   tau, 
                  work, *lwork, info);
}

void MAGMAF_ZGEQLF( magma_int_t *m, magma_int_t *n, 
                    cuFloatComplex *A,    magma_int_t *lda,   cuFloatComplex *tau, 
                    cuFloatComplex *work, magma_int_t *lwork, magma_int_t *info)
{
    magma_cgeqlf( *m, *n, 
                  A,    *lda,   tau, 
                  work, *lwork, info);
}

void MAGMAF_CGEQRF( magma_int_t *m, magma_int_t *n, cuFloatComplex *A, 
                    magma_int_t *lda, cuFloatComplex *tau, cuFloatComplex *work, 
                    magma_int_t *lwork, magma_int_t *info)
{
    magma_cgeqrf( *m, *n, A, 
                  *lda, tau, work, 
                  *lwork, info);
}

void MAGMAF_CGESV ( magma_int_t *n, magma_int_t *nrhs,
                    cuFloatComplex *A, magma_int_t *lda, magma_int_t *ipiv,
                    cuFloatComplex *B, magma_int_t *ldb, magma_int_t *info)
{
    magma_cgesv(  *n, *nrhs,
                  A, *lda, ipiv,
                  B, *ldb,
                  info);
}
    
void MAGMAF_CGETRF( magma_int_t *m, magma_int_t *n, cuFloatComplex *A, 
                    magma_int_t *lda, magma_int_t *ipiv, 
                    magma_int_t *info)
{
    magma_cgetrf( *m, *n, A, 
                  *lda, ipiv, 
                  info);
}

// void MAGMAF_CLATRD( char *uplo, magma_int_t *n, magma_int_t *nb, cuFloatComplex *a, 
//                     magma_int_t *lda, float *e, cuFloatComplex *tau, 
//                     cuFloatComplex *w, magma_int_t *ldw,
//                     cuFloatComplex *da, magma_int_t *ldda, 
//                     cuFloatComplex *dw, magma_int_t *lddw)
// {
//     magma_clatrd( uplo[0], *n, *nb, a, 
//                   *lda, e, tau, 
//                   w, *ldw,
//                   da, *ldda, 
//                   dw, *lddw);
// }

  /* This has nothing to do here, it should be a GPU function */
// void MAGMAF_ZLAHR2( magma_int_t *m, magma_int_t *n, magma_int_t *nb, 
//                     cuFloatComplex *da, cuFloatComplex *dv, cuFloatComplex *a, 
//                     magma_int_t *lda, cuFloatComplex *tau, cuFloatComplex *t, 
//                     magma_int_t *ldt, cuFloatComplex *y, magma_int_t *ldy)
// {
//     magma_clahr2( *m, *n, *nb, 
//                   da, dv, a, 
//                   *lda, tau, t, 
//                   *ldt, y, *ldy);
// }

// void MAGMAF_ZLAHRU( magma_int_t *m, magma_int_t *n, magma_int_t *nb, 
//                     cuFloatComplex *a, magma_int_t *lda, 
//                     cuFloatComplex *da, cuFloatComplex *y, 
//                     cuFloatComplex *v, cuFloatComplex *t, 
//                     cuFloatComplex *dwork)
// {
//     magma_clahru( *m, *n, *nb, 
//                   a, *lda, 
//                   da, y, 
//                   v, t, 
//                   dwork);
// }

void MAGMAF_CPOSV(  char *uplo, magma_int_t *n, magma_int_t *nrhs,
                    cuFloatComplex *A, magma_int_t *lda,
                    cuFloatComplex *B, magma_int_t *ldb, magma_int_t *info)
{
    magma_cposv(  uplo[0], *n, *nrhs,
                  A, *lda,
                  B, *ldb, info);
}

void MAGMAF_CPOTRF( char *uplo, magma_int_t *n, cuFloatComplex *A, 
                    magma_int_t *lda, magma_int_t *info)
{
    magma_cpotrf( uplo[0], *n, A, 
                  *lda, info);
}

void MAGMAF_CHETRD( char *uplo, magma_int_t *n, cuFloatComplex *A, 
                    magma_int_t *lda, float *d, float *e, 
                    cuFloatComplex *tau, cuFloatComplex *work, magma_int_t *lwork, 
                    magma_int_t *info)
{
    magma_chetrd( uplo[0], *n, A, 
                  *lda, d, e, 
                  tau, work, *lwork, 
                  info);
}

// void MAGMAF_CUNGQR( magma_int_t *m, magma_int_t *n, magma_int_t *k,
//                     cuFloatComplex *a, magma_int_t *lda,
//                     cuFloatComplex *tau, cuFloatComplex *dwork,
//                     magma_int_t *nb, magma_int_t *info )
// {
//     magma_cungqr( *m, *n, *k,
//                   a, *lda,
//                   tau, dwork,
//                   *nb, info );
// }

void MAGMAF_CUNMQR( char *side, char *trans, 
                    magma_int_t *m, magma_int_t *n, magma_int_t *k, 
                    cuFloatComplex *a, magma_int_t *lda, cuFloatComplex *tau, 
                    cuFloatComplex *c, magma_int_t *ldc, 
                    cuFloatComplex *work, magma_int_t *lwork, magma_int_t *info)
{
    magma_cunmqr( side[0], trans[0], 
                  *m, *n, *k, 
                  a, *lda, tau, 
                  c, *ldc, 
                  work, *lwork, info);
}

void MAGMAF_CUNMTR( char *side, char *uplo, char *trans,
                    magma_int_t *m, magma_int_t *n,
                    cuFloatComplex *a,    magma_int_t *lda,
                    cuFloatComplex *tau,
                    cuFloatComplex *c,    magma_int_t *ldc,
                    cuFloatComplex *work, magma_int_t *lwork,
                    magma_int_t *info)
{
    magma_cunmtr( side[0], uplo[0], trans[0],
                  *m, *n,
                  a,    *lda,
                  tau,
                  c,    *ldc,
                  work, *lwork,
                  info);
}

// void MAGMAF_CUNGHR( magma_int_t *n, magma_int_t *ilo, magma_int_t *ihi,
//                     cuFloatComplex *a, magma_int_t *lda,
//                     cuFloatComplex *tau,
//                     cuFloatComplex *dT, magma_int_t *nb,
//                     magma_int_t *info)
// {
//     magma_cunghr( *n, *ilo, *ihi,
//                   a, *lda,
//                   tau,
//                   dT, *nb,
//                   info);
// }

#if defined(PRECISION_z) || defined(PRECISION_c)
void MAGMAF_CGEEV( char *jobvl, char *jobvr, magma_int_t *n,
                   cuFloatComplex *a, magma_int_t *lda,
                   cuFloatComplex *w,
                   cuFloatComplex *vl, magma_int_t *ldvl,
                   cuFloatComplex *vr, magma_int_t *ldvr,
                   cuFloatComplex *work, magma_int_t *lwork,
                   float *rwork, magma_int_t *info)
{
    magma_cgeev( jobvl[0], jobvr[0], *n,
                 a, *lda,
                 w,
                 vl, *ldvl,
                 vr, *ldvr,
                 work, *lwork,
                 rwork, info);
}

void MAGMAF_CGESVD( char *jobu, char *jobvt, magma_int_t *m, magma_int_t *n,
                    cuFloatComplex *a,    magma_int_t *lda, float *s, 
                    cuFloatComplex *u,    magma_int_t *ldu, 
                    cuFloatComplex *vt,   magma_int_t *ldvt,
                    cuFloatComplex *work, magma_int_t *lwork,
                    float *rwork, magma_int_t *info )
{
    magma_cgesvd( jobu[0], jobvt[0], *m, *n,
                  a,    *lda, s, 
                  u,    *ldu, 
                  vt,   *ldvt,
                  work, *lwork,
                  rwork, info );
}
    
void MAGMAF_CHEEVD( char *jobz, char *uplo, magma_int_t *n,
                    cuFloatComplex *a,     magma_int_t *lda, float *w,
                    cuFloatComplex *work,  magma_int_t *lwork,
                    float          *rwork, magma_int_t *lrwork,
                    magma_int_t *iwork, magma_int_t *liwork, magma_int_t *info)
{
    magma_cheevd( jobz[0], uplo[0], *n,
                  a, *lda, w,
                  work, *lwork,
                  rwork, *lrwork,
                  iwork, *liwork, info);
}

void MAGMAF_CHEGVD(magma_int_t *itype, char *jobz, char *uplo, magma_int_t *n,
                   cuFloatComplex *a, magma_int_t *lda, 
                   cuFloatComplex *b, magma_int_t *ldb,
                   float *w, cuFloatComplex *work, magma_int_t *lwork,
                   float *rwork, magma_int_t *lrwork,
                   magma_int_t *iwork, magma_int_t *liwork, magma_int_t *info)
{
  magma_chegvd( *itype, jobz[0], uplo[0], *n,
                a, *lda, b, *ldb,
                w, work, *lwork,
                rwork, *lrwork,
                iwork, *liwork, info);
}
    
#else
void MAGMAF_CGEEV( char *jobvl, char *jobvr, magma_int_t *n,
                   cuFloatComplex *a,    magma_int_t *lda,
                   cuFloatComplex *wr, cuFloatComplex *wi,
                   cuFloatComplex *vl,   magma_int_t *ldvl,
                   cuFloatComplex *vr,   magma_int_t *ldvr,
                   cuFloatComplex *work, magma_int_t *lwork,
                   magma_int_t *info)
{
    magma_cgeev( jobvl[0], jobvr[0], *n,
                 a,    *lda,
                 wr, wi,
                 vl,   *ldvl,
                 vr,   *ldvr,
                 work, *lwork,
                 info);
}

void MAGMAF_CGESVD( char *jobu, char *jobvt, magma_int_t *m, magma_int_t *n,
                    cuFloatComplex *a,    magma_int_t *lda, float *s,
                    cuFloatComplex *u,    magma_int_t *ldu, 
                    cuFloatComplex *vt,   magma_int_t *ldvt,
                    cuFloatComplex *work, magma_int_t *lwork,
                    magma_int_t *info )
{
    magma_cgesvd( jobu[0], jobvt[0], *m, *n,
                  a,    *lda, s, 
                  u,    *ldu, 
                  vt,   *ldvt,
                  work, *lwork,
                  info );
}

void MAGMAF_CHEEVD( char *jobz, char *uplo, magma_int_t *n,
                    cuFloatComplex *a, magma_int_t *lda, float *w,
                    cuFloatComplex *work, magma_int_t *lwork,
                    magma_int_t *iwork, magma_int_t *liwork, magma_int_t *info)
{
    magma_cheevd( jobz[0], uplo[0], *n,
                  a, *lda, w,
                  work, *lwork,
                  iwork, *liwork, info);
}

void MAGMAF_CHEGVD(magma_int_t *itype, char *jobz, char *uplo, magma_int_t *n,
                   cuFloatComplex *a, magma_int_t *lda,
                   cuFloatComplex *b, magma_int_t *ldb,
                   float *w, cuFloatComplex *work, magma_int_t *lwork,
                   magma_int_t *iwork, magma_int_t *liwork, magma_int_t *info)
{
  magma_chegvd( *itype, jobz[0], uplo[0], *n,
                a, *lda, b, *ldb,
                w, work, *lwork,
                iwork, *liwork, info);
}


#endif

/* //////////////////////////////////////////////////////////////////////////// 
 -- MAGMA function definitions / Data on GPU
*/
void MAGMAF_CGELS_GPU(  char *trans, magma_int_t *m, magma_int_t *n, magma_int_t *nrhs,
                        devptr_t *dA,    magma_int_t *ldda, 
                        devptr_t *dB,    magma_int_t *lddb, 
                        cuFloatComplex *hwork, magma_int_t *lwork, 
                        magma_int_t *info)
{
    magma_cgels_gpu(  trans[0], *m, *n, *nrhs, 
                      DEVPTR(dA),    *ldda,  
                      DEVPTR(dB),    *lddb,  
                      hwork, *lwork,  info);
}

void MAGMAF_CGEQRF_GPU( magma_int_t *m, magma_int_t *n, 
                        devptr_t *dA,  magma_int_t *ldda, 
                        cuFloatComplex *tau, devptr_t *dT, 
                        magma_int_t *info)
{
    magma_cgeqrf_gpu( *m, *n,  
                      DEVPTR(dA),  *ldda,  
                      tau, 
                      DEVPTR(dT),  info);
}

void MAGMAF_CGEQRF2_GPU(magma_int_t *m, magma_int_t *n, 
                        devptr_t *dA,  magma_int_t *ldda, 
                        cuFloatComplex *tau, magma_int_t *info)
{
    magma_cgeqrf2_gpu(*m, *n,  
                      DEVPTR(dA),  *ldda,  
                      tau, info); 
}

void MAGMAF_CGEQRF3_GPU(magma_int_t *m, magma_int_t *n, 
                        devptr_t *dA,  magma_int_t *ldda, 
                        cuFloatComplex *tau, devptr_t *dT,
                        magma_int_t *info)
{
    magma_cgeqrf3_gpu(*m, *n,  
                      DEVPTR(dA),  *ldda,  
                      tau, DEVPTR(dT), info); 
}

void MAGMAF_CGEQRS_GPU( magma_int_t *m, magma_int_t *n, magma_int_t *nrhs, 
                        devptr_t *dA,     magma_int_t *ldda, 
                        cuFloatComplex *tau,   devptr_t *dT,
                        devptr_t *dB,    magma_int_t *lddb,
                        cuFloatComplex *hwork, magma_int_t *lhwork, 
                        magma_int_t *info)
{
    magma_cgeqrs_gpu( *m, *n, *nrhs,  
                      DEVPTR(dA),     *ldda,  
                      tau,
                      DEVPTR(dT), 
                      DEVPTR(dB),    *lddb, 
                      hwork, *lhwork,  info);
}

void MAGMAF_CGEQRS3_GPU(magma_int_t *m, magma_int_t *n, magma_int_t *nrhs, 
                        devptr_t *dA,     magma_int_t *ldda, 
                        cuFloatComplex *tau,   devptr_t *dT,
                        devptr_t *dB,    magma_int_t *lddb,
                        cuFloatComplex *hwork, magma_int_t *lhwork, 
                        magma_int_t *info)
{
    magma_cgeqrs3_gpu(*m, *n, *nrhs,  
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
    magma_cgessm_gpu( storev[0], *m, *n, *k, *ib, ipiv,  
                      DEVPTR(dL1), *lddl1,  
                      DEVPTR(dL),  *lddl,  
                      DEVPTR(dA),  *ldda,  info);
}

void MAGMAF_CGESV_GPU(  magma_int_t *n, magma_int_t *nrhs, 
                        devptr_t *dA, magma_int_t *ldda, magma_int_t *ipiv, 
                        devptr_t *dB, magma_int_t *lddb, magma_int_t *info)
{
    magma_cgesv_gpu(  *n, *nrhs,  
                      DEVPTR(dA), *ldda, ipiv,  
                      DEVPTR(dB), *lddb, info);
}

void MAGMAF_CGETRF_GPU( magma_int_t *m, magma_int_t *n, 
                        devptr_t *dA, magma_int_t *ldda, 
                        magma_int_t *ipiv, magma_int_t *info)
{
    magma_cgetrf_gpu( *m, *n,  
                      DEVPTR(dA), *ldda, ipiv, info);
}

void MAGMAF_CGETRS_GPU( char *trans, magma_int_t *n, magma_int_t *nrhs, 
                        devptr_t *dA, magma_int_t *ldda, magma_int_t *ipiv, 
                        devptr_t *dB, magma_int_t *lddb, magma_int_t *info)
{
    magma_cgetrs_gpu( trans[0], *n, *nrhs,  
                      DEVPTR(dA), *ldda, ipiv,  
                      DEVPTR(dB), *lddb, info);
}

void MAGMAF_CLABRD_GPU( magma_int_t *m, magma_int_t *n, magma_int_t *nb, 
                        cuFloatComplex *a, magma_int_t *lda, devptr_t *da, magma_int_t *ldda,
                        float *d, float *e, cuFloatComplex *tauq, cuFloatComplex *taup,  
                        cuFloatComplex *x, magma_int_t *ldx, devptr_t *dx, magma_int_t *lddx, 
                        cuFloatComplex *y, magma_int_t *ldy, devptr_t *dy, magma_int_t *lddy)
{
    magma_clabrd_gpu( *m, *n, *nb,  
                      a, *lda, DEVPTR(da), *ldda, 
                      d, e, tauq, taup,   
                      x, *ldx, DEVPTR(dx), *lddx,  
                      y, *ldy, DEVPTR(dy), *lddy);
}

void MAGMAF_CLARFB_GPU( char *side, char *trans, char *direct, char *storev, 
                        magma_int_t *m, magma_int_t *n, magma_int_t *k,
                        devptr_t *dv, magma_int_t *ldv, devptr_t *dt,    magma_int_t *ldt, 
                        devptr_t *dc, magma_int_t *ldc, devptr_t *dowrk, magma_int_t *ldwork )
{
    magma_clarfb_gpu( side[0], trans[0], direct[0], storev[0],  *m, *n, *k, 
                      DEVPTR(dv), *ldv, DEVPTR(dt),    *ldt,  
                      DEVPTR(dc), *ldc, DEVPTR(dowrk), *ldwork);
}

void MAGMAF_CPOSV_GPU(  char *uplo, magma_int_t *n, magma_int_t *nrhs, 
                        devptr_t *dA, magma_int_t *ldda, 
                        devptr_t *dB, magma_int_t *lddb, magma_int_t *info)
{
    magma_cposv_gpu(  uplo[0], *n, *nrhs,  
                      DEVPTR(dA), *ldda,  
                      DEVPTR(dB), *lddb, info);
}

void MAGMAF_CPOTRF_GPU( char *uplo,  magma_int_t *n, 
                        devptr_t *dA, magma_int_t *ldda, magma_int_t *info)
{
    magma_cpotrf_gpu( uplo[0],  *n,  
                      DEVPTR(dA), *ldda, info); }

void MAGMAF_CPOTRS_GPU( char *uplo,  magma_int_t *n, magma_int_t *nrhs, 
                        devptr_t *dA, magma_int_t *ldda, 
                        devptr_t *dB, magma_int_t *lddb, magma_int_t *info)
{
    magma_cpotrs_gpu( uplo[0],  *n, *nrhs,  
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
    magma_cssssm_gpu( storev[0], *m1, *n1,  *m2, *n2, *k, *ib,  
                      DEVPTR(dA1), *ldda1,  
                      DEVPTR(dA2), *ldda2,  
                      DEVPTR(dL1), *lddl1,  
                      DEVPTR(dL2), *lddl2,
                      IPIV, info);
}

void MAGMAF_CUNGQR_GPU( magma_int_t *m, magma_int_t *n, magma_int_t *k, 
                        devptr_t *da, magma_int_t *ldda, 
                        cuFloatComplex *tau, devptr_t *dwork, 
                        magma_int_t *nb, magma_int_t *info )
{
    magma_cungqr_gpu( *m, *n, *k,  
                      DEVPTR(da), *ldda, tau, 
                      DEVPTR(dwork), *nb, info );
}

void MAGMAF_CUNMQR_GPU( char *side, char *trans, 
                        magma_int_t *m, magma_int_t *n, magma_int_t *k,
                        devptr_t *a,    magma_int_t *lda, cuFloatComplex *tau, 
                        devptr_t *c,    magma_int_t *ldc,
                        devptr_t *work, magma_int_t *lwork, 
                        devptr_t *td,   magma_int_t *nb, magma_int_t *info)
{
    magma_cunmqr_gpu( side[0], trans[0], *m, *n, *k, 
                      DEVPTR(a),    *lda, tau,  
                      DEVPTR(c),    *ldc, 
                      DEVPTR(work), *lwork,  
                      DEVPTR(td),   *nb, info);
}

#ifdef __cplusplus
}
#endif
