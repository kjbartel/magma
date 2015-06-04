/*
 *   -- MAGMA (version 1.3.0) --
 *      Univ. of Tennessee, Knoxville
 *      Univ. of California, Berkeley
 *      Univ. of Colorado, Denver
 *      November 2012
 *
 * @generated s Wed Nov 14 22:52:26 2012
 */

#ifndef _MAGMA_S_H_
#define _MAGMA_S_H_
#define PRECISION_s

#ifdef __cplusplus
extern "C" {
#endif

/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA Auxiliary functions to get the NB used
*/
int magma_get_spotrf_nb(int m);
int magma_get_sgetrf_nb(int m);
int magma_get_sgetri_nb(int m);
int magma_get_sgeqp3_nb(int m);
int magma_get_sgeqrf_nb(int m);
int magma_get_sgeqlf_nb(int m);
int magma_get_sgehrd_nb(int m);
int magma_get_ssytrd_nb(int m);
int magma_get_sgelqf_nb(int m);
int magma_get_sgebrd_nb(int m);
int magma_get_ssygst_nb(int m);
int magma_get_sgesvd_nb(int m);

/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA function definitions / Data on CPU
*/
magma_int_t magma_sgebrd( magma_int_t m, magma_int_t n, float *A, 
                          magma_int_t lda, float *d, float *e,
                          float *tauq,  float *taup, 
                          float *work, magma_int_t lwork, magma_int_t *info);
magma_int_t magma_sgehrd2(magma_int_t n, magma_int_t ilo, magma_int_t ihi,
                          float *A, magma_int_t lda, float *tau, 
                          float *work, magma_int_t *lwork, magma_int_t *info);
magma_int_t magma_sgehrd( magma_int_t n, magma_int_t ilo, magma_int_t ihi,
                          float *A, magma_int_t lda, float *tau,
                          float *work, magma_int_t lwork,
                          float *d_T, magma_int_t *info);
magma_int_t magma_sgelqf( magma_int_t m, magma_int_t n, 
                          float *A,    magma_int_t lda,   float *tau, 
                          float *work, magma_int_t lwork, magma_int_t *info);
magma_int_t magma_sgeqlf( magma_int_t m, magma_int_t n, 
                          float *A,    magma_int_t lda,   float *tau, 
                          float *work, magma_int_t lwork, magma_int_t *info);
magma_int_t magma_sgeqrf( magma_int_t m, magma_int_t n, float *A, 
                          magma_int_t lda, float *tau, float *work, 
                          magma_int_t lwork, magma_int_t *info);
magma_int_t magma_sgeqr2( magma_int_t *m, magma_int_t *n, float *a,
                          magma_int_t *lda, float *tau, float *work, 
                          magma_int_t *info);
magma_int_t magma_sgeqrf4(magma_int_t num_gpus, magma_int_t m, magma_int_t n,
                          float *a,    magma_int_t lda, float *tau,
                          float *work, magma_int_t lwork, magma_int_t *info );
magma_int_t magma_sgeqrf_ooc( magma_int_t m, magma_int_t n, float *A,
                          magma_int_t lda, float *tau, float *work,
                          magma_int_t lwork, magma_int_t *info);
magma_int_t magma_sgesv ( magma_int_t n, magma_int_t nrhs, 
                          float *A, magma_int_t lda, magma_int_t *ipiv, 
                          float *B, magma_int_t ldb, magma_int_t *info);
magma_int_t magma_sgetrf( magma_int_t m, magma_int_t n, float *A, 
                          magma_int_t lda, magma_int_t *ipiv, 
                          magma_int_t *info);
magma_int_t magma_sgetrf2(magma_int_t m, magma_int_t n, float *a, 
                          magma_int_t lda, magma_int_t *ipiv, magma_int_t *info);

magma_int_t magma_slaqps( magma_int_t m, magma_int_t n, magma_int_t offset, 
                          magma_int_t nb, magma_int_t *kb, 
                          float *A,  magma_int_t lda,
                          float *dA, magma_int_t ldda,
                          magma_int_t *jpvt, float *tau, float *vn1, float *vn2, 
                          float *auxv, 
                          float *F,  magma_int_t ldf,
                          float *dF, magma_int_t lddf );
magma_int_t magma_slarf(  char *side, magma_int_t *m, magma_int_t *n,
                          float *v, magma_int_t *incv, float *tau,
                          float *c__, magma_int_t *ldc, float *work);
magma_int_t magma_slarfg( magma_int_t *n, float *alpha, float *x,
                          magma_int_t *incx, float *tau);
magma_int_t magma_slatrd( char uplo, magma_int_t n, magma_int_t nb, float *a, 
                          magma_int_t lda, float *e, float *tau, 
                          float *w, magma_int_t ldw,
                          float *da, magma_int_t ldda, 
                          float *dw, magma_int_t lddw);
magma_int_t magma_slatrd2(char uplo, magma_int_t n, magma_int_t nb,
                          float *a,  magma_int_t lda,
                          float *e, float *tau,
                          float *w,  magma_int_t ldw,
                          float *da, magma_int_t ldda,
                          float *dw, magma_int_t lddw,
                          float *dwork, magma_int_t ldwork);
magma_int_t magma_slahr2( magma_int_t m, magma_int_t n, magma_int_t nb, 
                          float *da, float *dv, float *a, 
                          magma_int_t lda, float *tau, float *t, 
                          magma_int_t ldt, float *y, magma_int_t ldy);
magma_int_t magma_slahru( magma_int_t n, magma_int_t ihi, magma_int_t k, magma_int_t nb, 
                          float *a, magma_int_t lda, 
                          float *da, float *y, 
                          float *v, float *t, 
                          float *dwork);
magma_int_t magma_sposv ( char uplo, magma_int_t n, magma_int_t nrhs, 
                          float *A, magma_int_t lda, 
                          float *B, magma_int_t ldb, magma_int_t *info);
magma_int_t magma_spotrf( char uplo, magma_int_t n, float *A, 
                          magma_int_t lda, magma_int_t *info);
magma_int_t magma_spotri( char uplo, magma_int_t n, float *A,
                          magma_int_t lda, magma_int_t *info);
magma_int_t magma_slauum( char uplo, magma_int_t n, float *A,
                          magma_int_t lda, magma_int_t *info);
magma_int_t magma_strtri( char uplo, char diag, magma_int_t n, float *A, 
                          magma_int_t lda, magma_int_t *info);
magma_int_t magma_ssytrd( char uplo, magma_int_t n, float *A, 
                          magma_int_t lda, float *d, float *e, 
                          float *tau, float *work, magma_int_t lwork, 
                          magma_int_t *info);
magma_int_t magma_sorgqr( magma_int_t m, magma_int_t n, magma_int_t k,
                          float *a, magma_int_t lda,
                          float *tau, float *dwork,
                          magma_int_t nb, magma_int_t *info );
magma_int_t magma_sormql( const char side, const char trans,
                          magma_int_t m, magma_int_t n, magma_int_t k,
                          float *a, magma_int_t lda,
                          float *tau,
                          float *c, magma_int_t ldc,
                          float *work, magma_int_t lwork,
                          magma_int_t *info);
magma_int_t magma_sormqr( char side, char trans, 
                          magma_int_t m, magma_int_t n, magma_int_t k, 
                          float *a, magma_int_t lda, float *tau, 
                          float *c, magma_int_t ldc, 
                          float *work, magma_int_t lwork, magma_int_t *info);
magma_int_t magma_sormtr( char side, char uplo, char trans,
                          magma_int_t m, magma_int_t n,
                          float *a,    magma_int_t lda,
                          float *tau,
                          float *c,    magma_int_t ldc,
                          float *work, magma_int_t lwork,
                          magma_int_t *info);
magma_int_t magma_sorghr( magma_int_t n, magma_int_t ilo, magma_int_t ihi,
                          float *a, magma_int_t lda,
                          float *tau,
                          float *dT, magma_int_t nb,
                          magma_int_t *info);
magma_int_t magma_ssyev( char jobz, char uplo, magma_int_t n,
                         float *a, magma_int_t lda, float *w,
                         float *work, magma_int_t lwork,
                         float *rwork, magma_int_t *info);
magma_int_t magma_ssyevx(char jobz, char range, char uplo, magma_int_t n,
                         float *a, magma_int_t lda, float vl, float vu,
                         magma_int_t il, magma_int_t iu, float abstol, magma_int_t *m,
                         float *w, float *z, magma_int_t ldz, 
                         float *work, magma_int_t lwork,
                         float *rwork, magma_int_t *iwork, 
                         magma_int_t *ifail, magma_int_t *info);
#if defined(PRECISION_z) || defined(PRECISION_c)
magma_int_t  magma_sgeev( char jobvl, char jobvr, magma_int_t n,
                          float *a, magma_int_t lda,
                          float *w,
                          float *vl, magma_int_t ldvl,
                          float *vr, magma_int_t ldvr,
                          float *work, magma_int_t lwork,
                          float *rwork, magma_int_t *info);
magma_int_t magma_sgeqp3( magma_int_t m, magma_int_t n,
                          float *a, magma_int_t lda, 
                          magma_int_t *jpvt, float *tau,
                          float *work, magma_int_t lwork, 
                          float *rwork, magma_int_t *info);
magma_int_t magma_sgesvd( char jobu, char jobvt, magma_int_t m, magma_int_t n,
                          float *a,    magma_int_t lda, float *s, 
                          float *u,    magma_int_t ldu, 
                          float *vt,   magma_int_t ldvt,
                          float *work, magma_int_t lwork,
                          float *rwork, magma_int_t *info );
magma_int_t magma_ssyevd( char jobz, char uplo, magma_int_t n,
                          float *a, magma_int_t lda, float *w,
                          float *work, magma_int_t lwork,
                          float *rwork, magma_int_t lrwork,
                          magma_int_t *iwork, magma_int_t liwork, magma_int_t *info);
magma_int_t magma_ssyevr( char jobz, char range, char uplo, magma_int_t n,
                          float *a, magma_int_t lda, float vl, float vu,
                          magma_int_t il, magma_int_t iu, float abstol, magma_int_t *m,
                          float *w, float *z, magma_int_t ldz, 
                          magma_int_t *isuppz,
                          float *work, magma_int_t lwork,
                          float *rwork, magma_int_t lrwork, magma_int_t *iwork,
                          magma_int_t liwork, magma_int_t *info);
magma_int_t magma_ssygvd( magma_int_t itype, char jobz, char uplo, magma_int_t n,
                          float *a, magma_int_t lda,
                          float *b, magma_int_t ldb,
                          float *w, float *work, magma_int_t lwork,
                          float *rwork, magma_int_t lrwork, magma_int_t *iwork,
                          magma_int_t liwork, magma_int_t *info);
magma_int_t magma_ssygvdx(magma_int_t itype, char jobz, char range, char uplo, 
                          magma_int_t n, float *a, magma_int_t lda,
                          float *b, magma_int_t ldb,
                          float vl, float vu, magma_int_t il, magma_int_t iu,
                          magma_int_t *m, float *w, float *work, 
                          magma_int_t lwork, float *rwork,
                          magma_int_t lrwork, magma_int_t *iwork,
                          magma_int_t liwork, magma_int_t *info);
magma_int_t magma_ssygvx( magma_int_t itype, char jobz, char range, char uplo, 
                          magma_int_t n, float *a, magma_int_t lda, 
                          float *b, magma_int_t ldb,
                          float vl, float vu, magma_int_t il, magma_int_t iu,
                          float abstol, magma_int_t *m, float *w, 
                          float *z, magma_int_t ldz,
                          float *work, magma_int_t lwork, float *rwork,
                          magma_int_t *iwork, magma_int_t *ifail, magma_int_t *info);
magma_int_t magma_ssygvr( magma_int_t itype, char jobz, char range, char uplo, 
                          magma_int_t n, float *a, magma_int_t lda,
                          float *b, magma_int_t ldb,
                          float vl, float vu, magma_int_t il, magma_int_t iu,
                          float abstol, magma_int_t *m, float *w, 
                          float *z, magma_int_t ldz,
                          magma_int_t *isuppz, float *work, magma_int_t lwork,
                          float *rwork, magma_int_t lrwork, magma_int_t *iwork,
                          magma_int_t liwork, magma_int_t *info);
magma_int_t magma_sstedx( char range, magma_int_t n, float vl, float vu,
                          magma_int_t il, magma_int_t iu, float *D, float *E,
                          float *Z, magma_int_t ldz,
                          float *rwork, magma_int_t ldrwork, magma_int_t *iwork,
                          magma_int_t liwork, float* dwork, magma_int_t *info);
#else
magma_int_t  magma_sgeev( char jobvl, char jobvr, magma_int_t n,
                          float *a,    magma_int_t lda,
                          float *wr, float *wi,
                          float *vl,   magma_int_t ldvl,
                          float *vr,   magma_int_t ldvr,
                          float *work, magma_int_t lwork,
                          magma_int_t *info);
magma_int_t magma_sgeqp3( magma_int_t m, magma_int_t n,
                          float *a, magma_int_t lda, 
                          magma_int_t *jpvt, float *tau,
                          float *work, magma_int_t lwork,
                          magma_int_t *info);
magma_int_t magma_sgesvd( char jobu, char jobvt, magma_int_t m, magma_int_t n,
                          float *a,    magma_int_t lda, float *s, 
                          float *u,    magma_int_t ldu, 
                          float *vt,   magma_int_t ldvt,
                          float *work, magma_int_t lwork,
                          magma_int_t *info );
magma_int_t magma_ssyevd( char jobz, char uplo, magma_int_t n,
                          float *a, magma_int_t lda, float *w,
                          float *work, magma_int_t lwork,
                          magma_int_t *iwork, magma_int_t liwork, magma_int_t *info);
magma_int_t magma_ssygvd( magma_int_t itype, char jobz, char uplo, magma_int_t n,
                          float *a, magma_int_t lda,
                          float *b, magma_int_t ldb,
                          float *w, float *work, magma_int_t lwork,
                          magma_int_t *iwork, magma_int_t liwork, magma_int_t *info);
magma_int_t magma_sstedx( char range, magma_int_t n, float vl, float vu,
                          magma_int_t il, magma_int_t iu, float* d, float* e,
                          float* z, magma_int_t ldz,
                          float* work, magma_int_t lwork,
                          magma_int_t* iwork, magma_int_t liwork,
                          float* dwork, magma_int_t* info);
magma_int_t magma_slaex0( magma_int_t n, float* d, float* e, float* q, magma_int_t ldq,
                          float* work, magma_int_t* iwork, float* dwork,
                          char range, float vl, float vu,
                          magma_int_t il, magma_int_t iu, magma_int_t* info);
magma_int_t magma_slaex1( magma_int_t n, float* d, float* q, magma_int_t ldq,
                          magma_int_t* indxq, float rho, magma_int_t cutpnt,
                          float* work, magma_int_t* iwork, float* dwork,
                          char range, float vl, float vu,
                          magma_int_t il, magma_int_t iu, magma_int_t* info);
magma_int_t magma_slaex3( magma_int_t k, magma_int_t n, magma_int_t n1, float* d,
                          float* q, magma_int_t ldq, float rho,
                          float* dlamda, float* q2, magma_int_t* indx,
                          magma_int_t* ctot, float* w, float* s, magma_int_t* indxq,
                          float* dwork,
                          char range, float vl, float vu, magma_int_t il, magma_int_t iu,
                          magma_int_t* info );
#endif

magma_int_t magma_ssygst( magma_int_t itype, char uplo, magma_int_t n,
                          float *a, magma_int_t lda,
                          float *b, magma_int_t ldb, magma_int_t *info);

/* //////////////////////////////////////////////////////////////////////////// 
 -- MAGMA function definitions / Data on GPU
*/
magma_int_t magma_sgels_gpu(  char trans, magma_int_t m, magma_int_t n, magma_int_t nrhs,
                              float *dA,    magma_int_t ldda, 
                              float *dB,    magma_int_t lddb, 
                              float *hwork, magma_int_t lwork, 
                              magma_int_t *info);
magma_int_t magma_sgels3_gpu( char trans, magma_int_t m, magma_int_t n, magma_int_t nrhs,
                              float *dA,    magma_int_t ldda,
                              float *dB,    magma_int_t lddb,
                              float *hwork, magma_int_t lwork,
                              magma_int_t *info);
magma_int_t magma_sgelqf_gpu( magma_int_t m, magma_int_t n,
                              float *dA,    magma_int_t ldda,   float *tau,
                              float *work, magma_int_t lwork, magma_int_t *info);
magma_int_t magma_sgeqrf_gpu( magma_int_t m, magma_int_t n, 
                              float *dA,  magma_int_t ldda, 
                              float *tau, float *dT, 
                              magma_int_t *info);
magma_int_t magma_sgeqrf2_gpu(magma_int_t m, magma_int_t n, 
                              float *dA,  magma_int_t ldda, 
                              float *tau, magma_int_t *info);
magma_int_t magma_sgeqrf2_mgpu(magma_int_t num_gpus, magma_int_t m, magma_int_t n,
                               float **dlA, magma_int_t ldda,
                               float *tau, magma_int_t *info );
magma_int_t magma_sgeqrf3_gpu(magma_int_t m, magma_int_t n, 
                              float *dA,  magma_int_t ldda, 
                              float *tau, float *dT, 
                              magma_int_t *info);
magma_int_t magma_sgeqrs_gpu( magma_int_t m, magma_int_t n, magma_int_t nrhs, 
                              float *dA,     magma_int_t ldda, 
                              float *tau,   float *dT,
                              float *dB,    magma_int_t lddb,
                              float *hwork, magma_int_t lhwork, 
                              magma_int_t *info);
magma_int_t magma_sgeqrs3_gpu( magma_int_t m, magma_int_t n, magma_int_t nrhs, 
                              float *dA,     magma_int_t ldda, 
                              float *tau,   float *dT,
                              float *dB,    magma_int_t lddb,
                              float *hwork, magma_int_t lhwork, 
                              magma_int_t *info);
magma_int_t magma_sgessm_gpu( char storev, magma_int_t m, magma_int_t n, magma_int_t k, magma_int_t ib, 
                              magma_int_t *ipiv, 
                              float *dL1, magma_int_t lddl1, 
                              float *dL,  magma_int_t lddl, 
                              float *dA,  magma_int_t ldda, 
                              magma_int_t *info);
magma_int_t magma_sgesv_gpu(  magma_int_t n, magma_int_t nrhs, 
                              float *dA, magma_int_t ldda, magma_int_t *ipiv, 
                              float *dB, magma_int_t lddb, magma_int_t *info);
magma_int_t magma_sgetrf_incpiv_gpu( char storev, magma_int_t m, magma_int_t n, magma_int_t ib,
                              float *hA, magma_int_t ldha, float *dA, magma_int_t ldda,
                              float *hL, magma_int_t ldhl, float *dL, magma_int_t lddl,
                              magma_int_t *ipiv, 
                              float *dwork, magma_int_t lddwork,
                              magma_int_t *info);
magma_int_t magma_sgetrf_gpu( magma_int_t m, magma_int_t n, 
                              float *dA, magma_int_t ldda, 
                              magma_int_t *ipiv, magma_int_t *info);
magma_int_t 
magma_sgetrf_nopiv_gpu      ( magma_int_t m, magma_int_t n,
                              float *dA, magma_int_t ldda,
                              magma_int_t *info);
magma_int_t magma_sgetri_gpu( magma_int_t n, 
                              float *dA, magma_int_t ldda, magma_int_t *ipiv, 
                              float *dwork, magma_int_t lwork, magma_int_t *info);
magma_int_t magma_sgetrs_gpu( char trans, magma_int_t n, magma_int_t nrhs, 
                              float *dA, magma_int_t ldda, magma_int_t *ipiv, 
                              float *dB, magma_int_t lddb, magma_int_t *info);
magma_int_t magma_slabrd_gpu( magma_int_t m, magma_int_t n, magma_int_t nb, 
                              float *a, magma_int_t lda, float *da, magma_int_t ldda,
                              float *d, float *e, float *tauq, float *taup,  
                              float *x, magma_int_t ldx, float *dx, magma_int_t lddx, 
                              float *y, magma_int_t ldy, float *dy, magma_int_t lddy);
magma_int_t magma_slarfb_gpu( char side, char trans, char direct, char storev, 
                              magma_int_t m, magma_int_t n, magma_int_t k,
                              const float *dv, magma_int_t ldv,
                              const float *dt, magma_int_t ldt, 
                              float *dc,       magma_int_t ldc,
                              float *dwork,    magma_int_t ldwork );
magma_int_t magma_sposv_gpu(  char uplo, magma_int_t n, magma_int_t nrhs, 
                              float *dA, magma_int_t ldda, 
                              float *dB, magma_int_t lddb, magma_int_t *info);
magma_int_t magma_spotrf_gpu( char uplo,  magma_int_t n, 
                              float *dA, magma_int_t ldda, magma_int_t *info);
magma_int_t magma_spotri_gpu( char uplo,  magma_int_t n,
                              float *dA, magma_int_t ldda, magma_int_t *info);
magma_int_t magma_slauum_gpu( char uplo,  magma_int_t n,
                              float *dA, magma_int_t ldda, magma_int_t *info);
magma_int_t magma_strtri_gpu( char uplo,  char diag, magma_int_t n,
                              float *dA, magma_int_t ldda, magma_int_t *info);
magma_int_t magma_ssytrd_gpu( char uplo, magma_int_t n,
                              float *da, magma_int_t ldda,
                              float *d, float *e, float *tau,
                              float *wa,  magma_int_t ldwa,
                              float *work, magma_int_t lwork,
                              magma_int_t *info);
magma_int_t magma_ssytrd2_gpu(char uplo, magma_int_t n,
                              float *da, magma_int_t ldda,
                              float *d, float *e, float *tau,
                              float *wa,  magma_int_t ldwa,
                              float *work, magma_int_t lwork,
                              float *dwork, magma_int_t ldwork,
                              magma_int_t *info);
magma_int_t magma_spotrs_gpu( char uplo,  magma_int_t n, magma_int_t nrhs, 
                              float *dA, magma_int_t ldda, 
                              float *dB, magma_int_t lddb, magma_int_t *info);
magma_int_t magma_sssssm_gpu( char storev, magma_int_t m1, magma_int_t n1, 
                              magma_int_t m2, magma_int_t n2, magma_int_t k, magma_int_t ib, 
                              float *dA1, magma_int_t ldda1, 
                              float *dA2, magma_int_t ldda2, 
                              float *dL1, magma_int_t lddl1, 
                              float *dL2, magma_int_t lddl2,
                              magma_int_t *IPIV, magma_int_t *info);
magma_int_t magma_ststrf_gpu( char storev, magma_int_t m, magma_int_t n, magma_int_t ib, magma_int_t nb,
                              float *hU, magma_int_t ldhu, float *dU, magma_int_t lddu, 
                              float *hA, magma_int_t ldha, float *dA, magma_int_t ldda, 
                              float *hL, magma_int_t ldhl, float *dL, magma_int_t lddl,
                              magma_int_t *ipiv, 
                              float *hwork, magma_int_t ldhwork, 
                              float *dwork, magma_int_t lddwork,
                              magma_int_t *info);
magma_int_t magma_sorgqr_gpu( magma_int_t m, magma_int_t n, magma_int_t k, 
                              float *da, magma_int_t ldda, 
                              float *tau, float *dwork, 
                              magma_int_t nb, magma_int_t *info );
magma_int_t magma_sormql2_gpu(const char side, const char trans,
                              magma_int_t m, magma_int_t n, magma_int_t k,
                              float *da, magma_int_t ldda,
                              float *tau,
                              float *dc, magma_int_t lddc,
                              float *wa, magma_int_t ldwa,
                              magma_int_t *info);
magma_int_t magma_sormqr_gpu( char side, char trans, 
                              magma_int_t m, magma_int_t n, magma_int_t k,
                              float *a,    magma_int_t lda, float *tau, 
                              float *c,    magma_int_t ldc,
                              float *work, magma_int_t lwork, 
                              float *td,   magma_int_t nb, magma_int_t *info);
magma_int_t magma_sormqr2_gpu(const char side, const char trans,
                              magma_int_t m, magma_int_t n, magma_int_t k,
                              float *da,   magma_int_t ldda,
                              float *tau,
                              float *dc,    magma_int_t lddc,
                              float *wa,    magma_int_t ldwa,
                              magma_int_t *info);
magma_int_t magma_sormtr_gpu( char side, char uplo, char trans,
                              magma_int_t m, magma_int_t n,
                              float *da,    magma_int_t ldda,
                              float *tau,
                              float *dc,    magma_int_t lddc,
                              float *wa,    magma_int_t ldwa,
                              magma_int_t *info);

#if defined(PRECISION_z) || defined(PRECISION_c)
magma_int_t magma_ssyevd_gpu( char jobz, char uplo,
                              magma_int_t n,
                              float *da, magma_int_t ldda,
                              float *w,
                              float *wa,  magma_int_t ldwa,
                              float *work, magma_int_t lwork,
                              float *rwork, magma_int_t lrwork,
                              magma_int_t *iwork, magma_int_t liwork,
                              magma_int_t *info);
magma_int_t magma_ssyevdx_gpu(char jobz, char range, char uplo,
                              magma_int_t n, float *da, 
                              magma_int_t ldda, float vl, float vu, 
                              magma_int_t il, magma_int_t iu,
                              magma_int_t *m, float *w,
                              float *wa,  magma_int_t ldwa,
                              float *work, magma_int_t lwork,
                              float *rwork, magma_int_t lrwork,
                              magma_int_t *iwork, magma_int_t liwork,
                              magma_int_t *info);
magma_int_t magma_ssyevr_gpu( char jobz, char range, char uplo, magma_int_t n,
                              float *da, magma_int_t ldda, float vl, float vu,
                              magma_int_t il, magma_int_t iu, float abstol, magma_int_t *m,
                              float *w, float *dz, magma_int_t lddz,
                              magma_int_t *isuppz,
                              float *wa, magma_int_t ldwa,
                              float *wz, magma_int_t ldwz,
                              float *work, magma_int_t lwork,
                              float *rwork, magma_int_t lrwork, magma_int_t *iwork,
                              magma_int_t liwork, magma_int_t *info);
#else
magma_int_t magma_ssyevd_gpu( char jobz, char uplo,
                              magma_int_t n,
                              float *da, magma_int_t ldda,
                              float *w,
                              float *wa,  magma_int_t ldwa,
                              float *work, magma_int_t lwork,
                              magma_int_t *iwork, magma_int_t liwork,
                              magma_int_t *info);
#endif

magma_int_t magma_ssyevx_gpu( char jobz, char range, char uplo, magma_int_t n,
                              float *da, magma_int_t ldda, float vl, 
                              float vu, magma_int_t il, magma_int_t iu, 
                              float abstol, magma_int_t *m,
                              float *w, float *dz, magma_int_t lddz,
                              float *wa, magma_int_t ldwa,
                              float *wz, magma_int_t ldwz,
                              float *work, magma_int_t lwork,
                              float *rwork, magma_int_t *iwork, 
                              magma_int_t *ifail, magma_int_t *info);
magma_int_t magma_ssygst_gpu(magma_int_t itype, char uplo, magma_int_t n,
                             float *da, magma_int_t ldda,
                             float *db, magma_int_t lddb, magma_int_t *info);


/* //////////////////////////////////////////////////////////////////////////// 
 -- MAGMA utility function definitions
*/

void magma_sprint    ( magma_int_t m, magma_int_t n, const float  *A, magma_int_t lda  );
void magma_sprint_gpu( magma_int_t m, magma_int_t n, const float *dA, magma_int_t ldda );

void spanel_to_q( char uplo, magma_int_t ib, float *A, magma_int_t lda, float *work );
void sq_to_panel( char uplo, magma_int_t ib, float *A, magma_int_t lda, float *work );

#ifdef __cplusplus
}
#endif

#undef PRECISION_s
#endif /* _MAGMA_S_H_ */
