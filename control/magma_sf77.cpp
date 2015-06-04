/*
    -- MAGMA (version 1.2.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       June 2012

       @generated s Thu Jun 28 12:30:05 2012

*/

#include "common_magma.h"

/* 
 * typedef comming from fortran.h file provided in $CUDADIR/src directory
 * it will probably change with future release of cublas when they will use 64bits address
 */
typedef size_t devptr_t;

#define PRECISION_s

#ifdef PGI_FORTRAN
#define DEVPTR(__ptr) ((float*)(__ptr))
#else
#define DEVPTR(__ptr) ((float*)(uintptr_t)(*(__ptr)))
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
#define MAGMAF_SGEBRD  MAGMA_FORTRAN_NAME(sgebrd,  SGEBRD ) 
#define MAGMAF_SGEHRD2 MAGMA_FORTRAN_NAME(sgehrd2, SGEHRD2)
#define MAGMAF_SGEHRD  MAGMA_FORTRAN_NAME(sgehrd,  SGEHRD )
#define MAGMAF_SGELQF  MAGMA_FORTRAN_NAME(sgelqf,  SGELQF )
#define MAGMAF_SGEQLF  MAGMA_FORTRAN_NAME(sgeqlf,  SGEQLF )
#define MAGMAF_SGEQRF  MAGMA_FORTRAN_NAME(sgeqrf,  SGEQRF )
#define MAGMAF_SGESV   MAGMA_FORTRAN_NAME(sgesv,   SGESV  )
#define MAGMAF_SGETRF  MAGMA_FORTRAN_NAME(sgetrf,  SGETRF )
#define MAGMAF_SSYGST  MAGMA_FORTRAN_NAME(ssygst,  SSYGST )
#define MAGMAF_SLATRD  MAGMA_FORTRAN_NAME(slatrd,  SLATRD )
#define MAGMAF_SLAHR2  MAGMA_FORTRAN_NAME(slahr2,  SLAHR2 )
#define MAGMAF_SLAHRU  MAGMA_FORTRAN_NAME(slahru,  SLAHRU )
#define MAGMAF_SPOSV   MAGMA_FORTRAN_NAME(sposv,   SPOSV  )
#define MAGMAF_SPOTRF  MAGMA_FORTRAN_NAME(spotrf,  SPOTRF )
#define MAGMAF_SSYTRD  MAGMA_FORTRAN_NAME(ssytrd,  SSYTRD )
#define MAGMAF_SORGQR  MAGMA_FORTRAN_NAME(sorgqr,  SORGQR )
#define MAGMAF_SORMQR  MAGMA_FORTRAN_NAME(sormqr,  SORMQR )
#define MAGMAF_SORMTR  MAGMA_FORTRAN_NAME(sormtr,  SORMTR )
#define MAGMAF_SORGHR  MAGMA_FORTRAN_NAME(sorghr,  SORGHR )
#define MAGMAF_SGEEV   MAGMA_FORTRAN_NAME(sgeev,   SGEEV  )
#define MAGMAF_SGESVD  MAGMA_FORTRAN_NAME(sgesvd,  SGESVD )
#define MAGMAF_SSYEVD  MAGMA_FORTRAN_NAME(ssyevd,  SSYEVD )
#define MAGMAF_SSYGVD  MAGMA_FORTRAN_NAME(ssygvd,  SSYGVD )

/* //////////////////////////////////////////////////////////////////////////// 
 -- MAGMA function definitions / Data on GPU
*/
#define MAGMAF_SGELS_GPU   MAGMA_GPU_FORTRAN_NAME(sgels,   SGELS  )
#define MAGMAF_SGEQRF_GPU  MAGMA_GPU_FORTRAN_NAME(sgeqrf,  SGEQRF ) 
#define MAGMAF_SGEQRF2_GPU MAGMA_GPU_FORTRAN_NAME(sgeqrf2, SGEQRF2)
#define MAGMAF_SGEQRF3_GPU MAGMA_GPU_FORTRAN_NAME(sgeqrf3, SGEQRF3)
#define MAGMAF_SGEQRS_GPU  MAGMA_GPU_FORTRAN_NAME(sgeqrs,  SGEQRS ) 
#define MAGMAF_SGEQRS3_GPU MAGMA_GPU_FORTRAN_NAME(sgeqrs3, SGEQRS3) 
#define MAGMAF_SGESSM_GPU  MAGMA_GPU_FORTRAN_NAME(sgessm,  SGESSM ) 
#define MAGMAF_SGESV_GPU   MAGMA_GPU_FORTRAN_NAME(sgesv,   SGESV  )  
#define MAGMAF_SGETRL_GPU  MAGMA_GPU_FORTRAN_NAME(sgetrl,  SGETRL ) 
#define MAGMAF_SGETRF_GPU  MAGMA_GPU_FORTRAN_NAME(sgetrf,  SGETRF ) 
#define MAGMAF_SGETRS_GPU  MAGMA_GPU_FORTRAN_NAME(sgetrs,  SGETRS ) 
#define MAGMAF_SSYGST_GPU  MAGMA_GPU_FORTRAN_NAME(ssygst,  SSYGST )
#define MAGMAF_SLABRD_GPU  MAGMA_GPU_FORTRAN_NAME(slabrd,  SLABRD ) 
#define MAGMAF_SLARFB_GPU  MAGMA_GPU_FORTRAN_NAME(slarfb,  SLARFB ) 
#define MAGMAF_SPOSV_GPU   MAGMA_GPU_FORTRAN_NAME(sposv,   SPOSV  )  
#define MAGMAF_SPOTRF_GPU  MAGMA_GPU_FORTRAN_NAME(spotrf,  SPOTRF ) 
#define MAGMAF_SPOTRS_GPU  MAGMA_GPU_FORTRAN_NAME(spotrs,  SPOTRS ) 
#define MAGMAF_SSSSSM_GPU  MAGMA_GPU_FORTRAN_NAME(sssssm,  SSSSSM ) 
#define MAGMAF_STSTRF_GPU  MAGMA_GPU_FORTRAN_NAME(ststrf,  STSTRF ) 
#define MAGMAF_SORGQR_GPU  MAGMA_GPU_FORTRAN_NAME(sorgqr,  SORGQR ) 
#define MAGMAF_SORMQR_GPU  MAGMA_GPU_FORTRAN_NAME(sormqr,  SORMQR ) 

#ifdef __cplusplus
extern "C" {
#endif

/***************************************************************************//**
 *  FORTRAN API - math functions (simple interface)
 **/

/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA function definitions / Data on CPU
*/
void MAGMAF_SGEBRD( magma_int_t *m, magma_int_t *n, float *A, 
                    magma_int_t *lda, float *d, float *e,
                    float *tauq, float *taup, 
                    float *work, magma_int_t *lwork, magma_int_t *info)
{
    magma_sgebrd( *m, *n, A, 
                  *lda, d, e,
                  tauq, taup, 
                  work, *lwork, info);
}
    
void MAGMAF_SGEHRD2(magma_int_t *n, magma_int_t *ilo, magma_int_t *ihi,
                    float *A, magma_int_t *lda, float *tau, 
                    float *work, magma_int_t *lwork, magma_int_t *info)
{
    magma_sgehrd2(*n, *ilo, *ihi,
                  A, *lda, tau, 
                  work, lwork, info);
}
    
void MAGMAF_SGEHRD( magma_int_t *n, magma_int_t *ilo, magma_int_t *ihi,
                    float *A, magma_int_t *lda, float *tau,
                    float *work, magma_int_t *lwork,
                    float *d_T, magma_int_t *info)
{
  magma_sgehrd( *n, *ilo, *ihi,
                A, *lda, tau,
                work, *lwork,
                d_T, info);
}

void MAGMAF_SGELQF( magma_int_t *m, magma_int_t *n, 
                    float *A,    magma_int_t *lda,   float *tau, 
                    float *work, magma_int_t *lwork, magma_int_t *info)
{
    magma_sgelqf( *m, *n, 
                  A,    *lda,   tau, 
                  work, *lwork, info);
}

void MAGMAF_SGEQLF( magma_int_t *m, magma_int_t *n, 
                    float *A,    magma_int_t *lda,   float *tau, 
                    float *work, magma_int_t *lwork, magma_int_t *info)
{
    magma_sgeqlf( *m, *n, 
                  A,    *lda,   tau, 
                  work, *lwork, info);
}

void MAGMAF_SGEQRF( magma_int_t *m, magma_int_t *n, float *A, 
                    magma_int_t *lda, float *tau, float *work, 
                    magma_int_t *lwork, magma_int_t *info)
{
    magma_sgeqrf( *m, *n, A, 
                  *lda, tau, work, 
                  *lwork, info);
}

void MAGMAF_SGESV ( magma_int_t *n, magma_int_t *nrhs,
                    float *A, magma_int_t *lda, magma_int_t *ipiv,
                    float *B, magma_int_t *ldb, magma_int_t *info)
{
    magma_sgesv(  *n, *nrhs,
                  A, *lda, ipiv,
                  B, *ldb,
                  info);
}
    
void MAGMAF_SGETRF( magma_int_t *m, magma_int_t *n, float *A, 
                    magma_int_t *lda, magma_int_t *ipiv, 
                    magma_int_t *info)
{
    magma_sgetrf( *m, *n, A, 
                  *lda, ipiv, 
                  info);
}

void MAGMAF_SSYGST( magma_int_t *itype, char *uplo, magma_int_t *n,
                    float *A, magma_int_t *lda,
                    float *B, magma_int_t *ldb, magma_int_t *info )
{
    magma_ssygst( *itype, *uplo, *n, A, *lda, B, *ldb, info );
}

// void MAGMAF_SLATRD( char *uplo, magma_int_t *n, magma_int_t *nb, float *a, 
//                     magma_int_t *lda, float *e, float *tau, 
//                     float *w, magma_int_t *ldw,
//                     float *da, magma_int_t *ldda, 
//                     float *dw, magma_int_t *lddw)
// {
//     magma_slatrd( uplo[0], *n, *nb, a, 
//                   *lda, e, tau, 
//                   w, *ldw,
//                   da, *ldda, 
//                   dw, *lddw);
// }

  /* This has nothing to do here, it should be a GPU function */
// void MAGMAF_SLAHR2( magma_int_t *m, magma_int_t *n, magma_int_t *nb, 
//                     float *da, float *dv, float *a, 
//                     magma_int_t *lda, float *tau, float *t, 
//                     magma_int_t *ldt, float *y, magma_int_t *ldy)
// {
//     magma_slahr2( *m, *n, *nb, 
//                   da, dv, a, 
//                   *lda, tau, t, 
//                   *ldt, y, *ldy);
// }

// void MAGMAF_SLAHRU( magma_int_t *m, magma_int_t *n, magma_int_t *nb, 
//                     float *a, magma_int_t *lda, 
//                     float *da, float *y, 
//                     float *v, float *t, 
//                     float *dwork)
// {
//     magma_slahru( *m, *n, *nb, 
//                   a, *lda, 
//                   da, y, 
//                   v, t, 
//                   dwork);
// }

void MAGMAF_SPOSV(  char *uplo, magma_int_t *n, magma_int_t *nrhs,
                    float *A, magma_int_t *lda,
                    float *B, magma_int_t *ldb, magma_int_t *info)
{
    magma_sposv(  uplo[0], *n, *nrhs,
                  A, *lda,
                  B, *ldb, info);
}

void MAGMAF_SPOTRF( char *uplo, magma_int_t *n, float *A, 
                    magma_int_t *lda, magma_int_t *info)
{
    magma_spotrf( uplo[0], *n, A, 
                  *lda, info);
}

void MAGMAF_SSYTRD( char *uplo, magma_int_t *n, float *A, 
                    magma_int_t *lda, float *d, float *e, 
                    float *tau, float *work, magma_int_t *lwork, 
                    magma_int_t *info)
{
    magma_ssytrd( uplo[0], *n, A, 
                  *lda, d, e, 
                  tau, work, *lwork, 
                  info);
}

// void MAGMAF_SORGQR( magma_int_t *m, magma_int_t *n, magma_int_t *k,
//                     float *a, magma_int_t *lda,
//                     float *tau, float *dwork,
//                     magma_int_t *nb, magma_int_t *info )
// {
//     magma_sorgqr( *m, *n, *k,
//                   a, *lda,
//                   tau, dwork,
//                   *nb, info );
// }

void MAGMAF_SORMQR( char *side, char *trans, 
                    magma_int_t *m, magma_int_t *n, magma_int_t *k, 
                    float *a, magma_int_t *lda, float *tau, 
                    float *c, magma_int_t *ldc, 
                    float *work, magma_int_t *lwork, magma_int_t *info)
{
    magma_sormqr( side[0], trans[0], 
                  *m, *n, *k, 
                  a, *lda, tau, 
                  c, *ldc, 
                  work, *lwork, info);
}

void MAGMAF_SORMTR( char *side, char *uplo, char *trans,
                    magma_int_t *m, magma_int_t *n,
                    float *a,    magma_int_t *lda,
                    float *tau,
                    float *c,    magma_int_t *ldc,
                    float *work, magma_int_t *lwork,
                    magma_int_t *info)
{
    magma_sormtr( side[0], uplo[0], trans[0],
                  *m, *n,
                  a,    *lda,
                  tau,
                  c,    *ldc,
                  work, *lwork,
                  info);
}

// void MAGMAF_SORGHR( magma_int_t *n, magma_int_t *ilo, magma_int_t *ihi,
//                     float *a, magma_int_t *lda,
//                     float *tau,
//                     float *dT, magma_int_t *nb,
//                     magma_int_t *info)
// {
//     magma_sorghr( *n, *ilo, *ihi,
//                   a, *lda,
//                   tau,
//                   dT, *nb,
//                   info);
// }

#if defined(PRECISION_z) || defined(PRECISION_c)
void MAGMAF_SGEEV( char *jobvl, char *jobvr, magma_int_t *n,
                   float *a, magma_int_t *lda,
                   float *w,
                   float *vl, magma_int_t *ldvl,
                   float *vr, magma_int_t *ldvr,
                   float *work, magma_int_t *lwork,
                   float *rwork, magma_int_t *info)
{
    magma_sgeev( jobvl[0], jobvr[0], *n,
                 a, *lda,
                 w,
                 vl, *ldvl,
                 vr, *ldvr,
                 work, *lwork,
                 rwork, info);
}

void MAGMAF_SGESVD( char *jobu, char *jobvt, magma_int_t *m, magma_int_t *n,
                    float *a,    magma_int_t *lda, float *s, 
                    float *u,    magma_int_t *ldu, 
                    float *vt,   magma_int_t *ldvt,
                    float *work, magma_int_t *lwork,
                    float *rwork, magma_int_t *info )
{
    magma_sgesvd( jobu[0], jobvt[0], *m, *n,
                  a,    *lda, s, 
                  u,    *ldu, 
                  vt,   *ldvt,
                  work, *lwork,
                  rwork, info );
}
    
void MAGMAF_SSYEVD( char *jobz, char *uplo, magma_int_t *n,
                    float *a,     magma_int_t *lda, float *w,
                    float *work,  magma_int_t *lwork,
                    float          *rwork, magma_int_t *lrwork,
                    magma_int_t *iwork, magma_int_t *liwork, magma_int_t *info)
{
    magma_ssyevd( jobz[0], uplo[0], *n,
                  a, *lda, w,
                  work, *lwork,
                  rwork, *lrwork,
                  iwork, *liwork, info);
}

void MAGMAF_SSYGVD(magma_int_t *itype, char *jobz, char *uplo, magma_int_t *n,
                   float *a, magma_int_t *lda, 
                   float *b, magma_int_t *ldb,
                   float *w, float *work, magma_int_t *lwork,
                   float *rwork, magma_int_t *lrwork,
                   magma_int_t *iwork, magma_int_t *liwork, magma_int_t *info)
{
  magma_ssygvd( *itype, jobz[0], uplo[0], *n,
                a, *lda, b, *ldb,
                w, work, *lwork,
                rwork, *lrwork,
                iwork, *liwork, info);
}
    
#else
void MAGMAF_SGEEV( char *jobvl, char *jobvr, magma_int_t *n,
                   float *a,    magma_int_t *lda,
                   float *wr, float *wi,
                   float *vl,   magma_int_t *ldvl,
                   float *vr,   magma_int_t *ldvr,
                   float *work, magma_int_t *lwork,
                   magma_int_t *info)
{
    magma_sgeev( jobvl[0], jobvr[0], *n,
                 a,    *lda,
                 wr, wi,
                 vl,   *ldvl,
                 vr,   *ldvr,
                 work, *lwork,
                 info);
}

void MAGMAF_SGESVD( char *jobu, char *jobvt, magma_int_t *m, magma_int_t *n,
                    float *a,    magma_int_t *lda, float *s,
                    float *u,    magma_int_t *ldu, 
                    float *vt,   magma_int_t *ldvt,
                    float *work, magma_int_t *lwork,
                    magma_int_t *info )
{
    magma_sgesvd( jobu[0], jobvt[0], *m, *n,
                  a,    *lda, s, 
                  u,    *ldu, 
                  vt,   *ldvt,
                  work, *lwork,
                  info );
}

void MAGMAF_SSYEVD( char *jobz, char *uplo, magma_int_t *n,
                    float *a, magma_int_t *lda, float *w,
                    float *work, magma_int_t *lwork,
                    magma_int_t *iwork, magma_int_t *liwork, magma_int_t *info)
{
    magma_ssyevd( jobz[0], uplo[0], *n,
                  a, *lda, w,
                  work, *lwork,
                  iwork, *liwork, info);
}

void MAGMAF_SSYGVD(magma_int_t *itype, char *jobz, char *uplo, magma_int_t *n,
                   float *a, magma_int_t *lda,
                   float *b, magma_int_t *ldb,
                   float *w, float *work, magma_int_t *lwork,
                   magma_int_t *iwork, magma_int_t *liwork, magma_int_t *info)
{
  magma_ssygvd( *itype, jobz[0], uplo[0], *n,
                a, *lda, b, *ldb,
                w, work, *lwork,
                iwork, *liwork, info);
}


#endif

/* //////////////////////////////////////////////////////////////////////////// 
 -- MAGMA function definitions / Data on GPU
*/
void MAGMAF_SGELS_GPU(  char *trans, magma_int_t *m, magma_int_t *n, magma_int_t *nrhs,
                        devptr_t *dA,    magma_int_t *ldda, 
                        devptr_t *dB,    magma_int_t *lddb, 
                        float *hwork, magma_int_t *lwork, 
                        magma_int_t *info)
{
    magma_sgels_gpu(  trans[0], *m, *n, *nrhs, 
                      DEVPTR(dA),    *ldda,  
                      DEVPTR(dB),    *lddb,  
                      hwork, *lwork,  info);
}

void MAGMAF_SGEQRF_GPU( magma_int_t *m, magma_int_t *n, 
                        devptr_t *dA,  magma_int_t *ldda, 
                        float *tau, devptr_t *dT, 
                        magma_int_t *info)
{
    magma_sgeqrf_gpu( *m, *n,  
                      DEVPTR(dA),  *ldda,  
                      tau, 
                      DEVPTR(dT),  info);
}

void MAGMAF_SGEQRF2_GPU(magma_int_t *m, magma_int_t *n, 
                        devptr_t *dA,  magma_int_t *ldda, 
                        float *tau, magma_int_t *info)
{
    magma_sgeqrf2_gpu(*m, *n,  
                      DEVPTR(dA),  *ldda,  
                      tau, info); 
}

void MAGMAF_SGEQRF3_GPU(magma_int_t *m, magma_int_t *n, 
                        devptr_t *dA,  magma_int_t *ldda, 
                        float *tau, devptr_t *dT,
                        magma_int_t *info)
{
    magma_sgeqrf3_gpu(*m, *n,  
                      DEVPTR(dA),  *ldda,  
                      tau, DEVPTR(dT), info); 
}

void MAGMAF_SGEQRS_GPU( magma_int_t *m, magma_int_t *n, magma_int_t *nrhs, 
                        devptr_t *dA,     magma_int_t *ldda, 
                        float *tau,   devptr_t *dT,
                        devptr_t *dB,    magma_int_t *lddb,
                        float *hwork, magma_int_t *lhwork, 
                        magma_int_t *info)
{
    magma_sgeqrs_gpu( *m, *n, *nrhs,  
                      DEVPTR(dA),     *ldda,  
                      tau,
                      DEVPTR(dT), 
                      DEVPTR(dB),    *lddb, 
                      hwork, *lhwork,  info);
}

void MAGMAF_SGEQRS3_GPU(magma_int_t *m, magma_int_t *n, magma_int_t *nrhs, 
                        devptr_t *dA,     magma_int_t *ldda, 
                        float *tau,   devptr_t *dT,
                        devptr_t *dB,    magma_int_t *lddb,
                        float *hwork, magma_int_t *lhwork, 
                        magma_int_t *info)
{
    magma_sgeqrs3_gpu(*m, *n, *nrhs,  
                      DEVPTR(dA),     *ldda,  
                      tau,
                      DEVPTR(dT), 
                      DEVPTR(dB),    *lddb, 
                      hwork, *lhwork,  info);
}

void MAGMAF_SGESSM_GPU( char *storev, magma_int_t *m, magma_int_t *n, magma_int_t *k, magma_int_t *ib, 
                        magma_int_t *ipiv, 
                        devptr_t *dL1, magma_int_t *lddl1, 
                        devptr_t *dL,  magma_int_t *lddl, 
                        devptr_t *dA,  magma_int_t *ldda, 
                        magma_int_t *info)
{
    magma_sgessm_gpu( storev[0], *m, *n, *k, *ib, ipiv,  
                      DEVPTR(dL1), *lddl1,  
                      DEVPTR(dL),  *lddl,  
                      DEVPTR(dA),  *ldda,  info);
}

void MAGMAF_SGESV_GPU(  magma_int_t *n, magma_int_t *nrhs, 
                        devptr_t *dA, magma_int_t *ldda, magma_int_t *ipiv, 
                        devptr_t *dB, magma_int_t *lddb, magma_int_t *info)
{
    magma_sgesv_gpu(  *n, *nrhs,  
                      DEVPTR(dA), *ldda, ipiv,  
                      DEVPTR(dB), *lddb, info);
}

void MAGMAF_SGETRF_GPU( magma_int_t *m, magma_int_t *n, 
                        devptr_t *dA, magma_int_t *ldda, 
                        magma_int_t *ipiv, magma_int_t *info)
{
    magma_sgetrf_gpu( *m, *n,  
                      DEVPTR(dA), *ldda, ipiv, info);
}

void MAGMAF_SGETRS_GPU( char *trans, magma_int_t *n, magma_int_t *nrhs, 
                        devptr_t *dA, magma_int_t *ldda, magma_int_t *ipiv, 
                        devptr_t *dB, magma_int_t *lddb, magma_int_t *info)
{
    magma_sgetrs_gpu( trans[0], *n, *nrhs,  
                      DEVPTR(dA), *ldda, ipiv,  
                      DEVPTR(dB), *lddb, info);
}

void MAGMAF_SSYGST_GPU( magma_int_t *itype, char *uplo, magma_int_t *n,
                        devptr_t *dA, magma_int_t *ldda,
                        devptr_t *dB, magma_int_t *lddb, magma_int_t *info )
{
    magma_ssygst_gpu( *itype, *uplo, *n,
                      DEVPTR(dA), *ldda,
                      DEVPTR(dB), *lddb, info );
}

void MAGMAF_SLABRD_GPU( magma_int_t *m, magma_int_t *n, magma_int_t *nb, 
                        float *a, magma_int_t *lda, devptr_t *da, magma_int_t *ldda,
                        float *d, float *e, float *tauq, float *taup,  
                        float *x, magma_int_t *ldx, devptr_t *dx, magma_int_t *lddx, 
                        float *y, magma_int_t *ldy, devptr_t *dy, magma_int_t *lddy)
{
    magma_slabrd_gpu( *m, *n, *nb,  
                      a, *lda, DEVPTR(da), *ldda, 
                      d, e, tauq, taup,   
                      x, *ldx, DEVPTR(dx), *lddx,  
                      y, *ldy, DEVPTR(dy), *lddy);
}

void MAGMAF_SLARFB_GPU( char *side, char *trans, char *direct, char *storev, 
                        magma_int_t *m, magma_int_t *n, magma_int_t *k,
                        devptr_t *dv, magma_int_t *ldv, devptr_t *dt,    magma_int_t *ldt, 
                        devptr_t *dc, magma_int_t *ldc, devptr_t *dowrk, magma_int_t *ldwork )
{
    magma_slarfb_gpu( side[0], trans[0], direct[0], storev[0],  *m, *n, *k, 
                      DEVPTR(dv), *ldv, DEVPTR(dt),    *ldt,  
                      DEVPTR(dc), *ldc, DEVPTR(dowrk), *ldwork);
}

void MAGMAF_SPOSV_GPU(  char *uplo, magma_int_t *n, magma_int_t *nrhs, 
                        devptr_t *dA, magma_int_t *ldda, 
                        devptr_t *dB, magma_int_t *lddb, magma_int_t *info)
{
    magma_sposv_gpu(  uplo[0], *n, *nrhs,  
                      DEVPTR(dA), *ldda,  
                      DEVPTR(dB), *lddb, info);
}

void MAGMAF_SPOTRF_GPU( char *uplo,  magma_int_t *n, 
                        devptr_t *dA, magma_int_t *ldda, magma_int_t *info)
{
    magma_spotrf_gpu( uplo[0],  *n,  
                      DEVPTR(dA), *ldda, info); }

void MAGMAF_SPOTRS_GPU( char *uplo,  magma_int_t *n, magma_int_t *nrhs, 
                        devptr_t *dA, magma_int_t *ldda, 
                        devptr_t *dB, magma_int_t *lddb, magma_int_t *info)
{
    magma_spotrs_gpu( uplo[0],  *n, *nrhs,  
                      DEVPTR(dA), *ldda,  
                      DEVPTR(dB), *lddb, info);
}

void MAGMAF_SSSSSM_GPU( char *storev, magma_int_t *m1, magma_int_t *n1, 
                        magma_int_t *m2, magma_int_t *n2, magma_int_t *k, magma_int_t *ib, 
                        devptr_t *dA1, magma_int_t *ldda1, 
                        devptr_t *dA2, magma_int_t *ldda2, 
                        devptr_t *dL1, magma_int_t *lddl1, 
                        devptr_t *dL2, magma_int_t *lddl2,
                        magma_int_t *IPIV, magma_int_t *info)
{
    magma_sssssm_gpu( storev[0], *m1, *n1,  *m2, *n2, *k, *ib,  
                      DEVPTR(dA1), *ldda1,  
                      DEVPTR(dA2), *ldda2,  
                      DEVPTR(dL1), *lddl1,  
                      DEVPTR(dL2), *lddl2,
                      IPIV, info);
}

void MAGMAF_SORGQR_GPU( magma_int_t *m, magma_int_t *n, magma_int_t *k, 
                        devptr_t *da, magma_int_t *ldda, 
                        float *tau, devptr_t *dwork, 
                        magma_int_t *nb, magma_int_t *info )
{
    magma_sorgqr_gpu( *m, *n, *k,  
                      DEVPTR(da), *ldda, tau, 
                      DEVPTR(dwork), *nb, info );
}

void MAGMAF_SORMQR_GPU( char *side, char *trans, 
                        magma_int_t *m, magma_int_t *n, magma_int_t *k,
                        devptr_t *a,    magma_int_t *lda, float *tau, 
                        devptr_t *c,    magma_int_t *ldc,
                        devptr_t *work, magma_int_t *lwork, 
                        devptr_t *td,   magma_int_t *nb, magma_int_t *info)
{
    magma_sormqr_gpu( side[0], trans[0], *m, *n, *k, 
                      DEVPTR(a),    *lda, tau,  
                      DEVPTR(c),    *ldc, 
                      DEVPTR(work), *lwork,  
                      DEVPTR(td),   *nb, info);
}

#ifdef __cplusplus
}
#endif
