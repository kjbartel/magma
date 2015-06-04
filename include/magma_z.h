/*
 *   -- MAGMA (version 1.3.0) --
 *      Univ. of Tennessee, Knoxville
 *      Univ. of California, Berkeley
 *      Univ. of Colorado, Denver
 *      November 2012
 *
 * @precisions normal z -> s d c
 */

#ifndef _MAGMA_Z_H_
#define _MAGMA_Z_H_
#define PRECISION_z

#ifdef __cplusplus
extern "C" {
#endif

/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA Auxiliary functions to get the NB used
*/
int magma_get_zpotrf_nb(int m);
int magma_get_zgetrf_nb(int m);
int magma_get_zgetri_nb(int m);
int magma_get_zgeqp3_nb(int m);
int magma_get_zgeqrf_nb(int m);
int magma_get_zgeqlf_nb(int m);
int magma_get_zgehrd_nb(int m);
int magma_get_zhetrd_nb(int m);
int magma_get_zgelqf_nb(int m);
int magma_get_zgebrd_nb(int m);
int magma_get_zhegst_nb(int m);
int magma_get_zgesvd_nb(int m);

/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA function definitions / Data on CPU
*/
magma_int_t magma_zgebrd( magma_int_t m, magma_int_t n, cuDoubleComplex *A, 
                          magma_int_t lda, double *d, double *e,
                          cuDoubleComplex *tauq,  cuDoubleComplex *taup, 
                          cuDoubleComplex *work, magma_int_t lwork, magma_int_t *info);
magma_int_t magma_zgehrd2(magma_int_t n, magma_int_t ilo, magma_int_t ihi,
                          cuDoubleComplex *A, magma_int_t lda, cuDoubleComplex *tau, 
                          cuDoubleComplex *work, magma_int_t *lwork, magma_int_t *info);
magma_int_t magma_zgehrd( magma_int_t n, magma_int_t ilo, magma_int_t ihi,
                          cuDoubleComplex *A, magma_int_t lda, cuDoubleComplex *tau,
                          cuDoubleComplex *work, magma_int_t lwork,
                          cuDoubleComplex *d_T, magma_int_t *info);
magma_int_t magma_zgelqf( magma_int_t m, magma_int_t n, 
                          cuDoubleComplex *A,    magma_int_t lda,   cuDoubleComplex *tau, 
                          cuDoubleComplex *work, magma_int_t lwork, magma_int_t *info);
magma_int_t magma_zgeqlf( magma_int_t m, magma_int_t n, 
                          cuDoubleComplex *A,    magma_int_t lda,   cuDoubleComplex *tau, 
                          cuDoubleComplex *work, magma_int_t lwork, magma_int_t *info);
magma_int_t magma_zgeqrf( magma_int_t m, magma_int_t n, cuDoubleComplex *A, 
                          magma_int_t lda, cuDoubleComplex *tau, cuDoubleComplex *work, 
                          magma_int_t lwork, magma_int_t *info);
magma_int_t magma_zgeqr2( magma_int_t *m, magma_int_t *n, cuDoubleComplex *a,
                          magma_int_t *lda, cuDoubleComplex *tau, cuDoubleComplex *work, 
                          magma_int_t *info);
magma_int_t magma_zgeqrf4(magma_int_t num_gpus, magma_int_t m, magma_int_t n,
                          cuDoubleComplex *a,    magma_int_t lda, cuDoubleComplex *tau,
                          cuDoubleComplex *work, magma_int_t lwork, magma_int_t *info );
magma_int_t magma_zgeqrf_ooc( magma_int_t m, magma_int_t n, cuDoubleComplex *A,
                          magma_int_t lda, cuDoubleComplex *tau, cuDoubleComplex *work,
                          magma_int_t lwork, magma_int_t *info);
magma_int_t magma_zgesv ( magma_int_t n, magma_int_t nrhs, 
                          cuDoubleComplex *A, magma_int_t lda, magma_int_t *ipiv, 
                          cuDoubleComplex *B, magma_int_t ldb, magma_int_t *info);
magma_int_t magma_zgetrf( magma_int_t m, magma_int_t n, cuDoubleComplex *A, 
                          magma_int_t lda, magma_int_t *ipiv, 
                          magma_int_t *info);
magma_int_t magma_zgetrf2(magma_int_t m, magma_int_t n, cuDoubleComplex *a, 
                          magma_int_t lda, magma_int_t *ipiv, magma_int_t *info);

magma_int_t magma_zlaqps( magma_int_t m, magma_int_t n, magma_int_t offset, 
                          magma_int_t nb, magma_int_t *kb, 
                          cuDoubleComplex *A,  magma_int_t lda,
                          cuDoubleComplex *dA, magma_int_t ldda,
                          magma_int_t *jpvt, cuDoubleComplex *tau, double *vn1, double *vn2, 
                          cuDoubleComplex *auxv, 
                          cuDoubleComplex *F,  magma_int_t ldf,
                          cuDoubleComplex *dF, magma_int_t lddf );
magma_int_t magma_zlarf(  char *side, magma_int_t *m, magma_int_t *n,
                          cuDoubleComplex *v, magma_int_t *incv, cuDoubleComplex *tau,
                          cuDoubleComplex *c__, magma_int_t *ldc, cuDoubleComplex *work);
magma_int_t magma_zlarfg( magma_int_t *n, cuDoubleComplex *alpha, cuDoubleComplex *x,
                          magma_int_t *incx, cuDoubleComplex *tau);
magma_int_t magma_zlatrd( char uplo, magma_int_t n, magma_int_t nb, cuDoubleComplex *a, 
                          magma_int_t lda, double *e, cuDoubleComplex *tau, 
                          cuDoubleComplex *w, magma_int_t ldw,
                          cuDoubleComplex *da, magma_int_t ldda, 
                          cuDoubleComplex *dw, magma_int_t lddw);
magma_int_t magma_zlatrd2(char uplo, magma_int_t n, magma_int_t nb,
                          cuDoubleComplex *a,  magma_int_t lda,
                          double *e, cuDoubleComplex *tau,
                          cuDoubleComplex *w,  magma_int_t ldw,
                          cuDoubleComplex *da, magma_int_t ldda,
                          cuDoubleComplex *dw, magma_int_t lddw,
                          cuDoubleComplex *dwork, magma_int_t ldwork);
magma_int_t magma_zlahr2( magma_int_t m, magma_int_t n, magma_int_t nb, 
                          cuDoubleComplex *da, cuDoubleComplex *dv, cuDoubleComplex *a, 
                          magma_int_t lda, cuDoubleComplex *tau, cuDoubleComplex *t, 
                          magma_int_t ldt, cuDoubleComplex *y, magma_int_t ldy);
magma_int_t magma_zlahru( magma_int_t n, magma_int_t ihi, magma_int_t k, magma_int_t nb, 
                          cuDoubleComplex *a, magma_int_t lda, 
                          cuDoubleComplex *da, cuDoubleComplex *y, 
                          cuDoubleComplex *v, cuDoubleComplex *t, 
                          cuDoubleComplex *dwork);
magma_int_t magma_zposv ( char uplo, magma_int_t n, magma_int_t nrhs, 
                          cuDoubleComplex *A, magma_int_t lda, 
                          cuDoubleComplex *B, magma_int_t ldb, magma_int_t *info);
magma_int_t magma_zpotrf( char uplo, magma_int_t n, cuDoubleComplex *A, 
                          magma_int_t lda, magma_int_t *info);
magma_int_t magma_zpotri( char uplo, magma_int_t n, cuDoubleComplex *A,
                          magma_int_t lda, magma_int_t *info);
magma_int_t magma_zlauum( char uplo, magma_int_t n, cuDoubleComplex *A,
                          magma_int_t lda, magma_int_t *info);
magma_int_t magma_ztrtri( char uplo, char diag, magma_int_t n, cuDoubleComplex *A, 
                          magma_int_t lda, magma_int_t *info);
magma_int_t magma_zhetrd( char uplo, magma_int_t n, cuDoubleComplex *A, 
                          magma_int_t lda, double *d, double *e, 
                          cuDoubleComplex *tau, cuDoubleComplex *work, magma_int_t lwork, 
                          magma_int_t *info);
magma_int_t magma_zungqr( magma_int_t m, magma_int_t n, magma_int_t k,
                          cuDoubleComplex *a, magma_int_t lda,
                          cuDoubleComplex *tau, cuDoubleComplex *dwork,
                          magma_int_t nb, magma_int_t *info );
magma_int_t magma_zunmql( const char side, const char trans,
                          magma_int_t m, magma_int_t n, magma_int_t k,
                          cuDoubleComplex *a, magma_int_t lda,
                          cuDoubleComplex *tau,
                          cuDoubleComplex *c, magma_int_t ldc,
                          cuDoubleComplex *work, magma_int_t lwork,
                          magma_int_t *info);
magma_int_t magma_zunmqr( char side, char trans, 
                          magma_int_t m, magma_int_t n, magma_int_t k, 
                          cuDoubleComplex *a, magma_int_t lda, cuDoubleComplex *tau, 
                          cuDoubleComplex *c, magma_int_t ldc, 
                          cuDoubleComplex *work, magma_int_t lwork, magma_int_t *info);
magma_int_t magma_zunmtr( char side, char uplo, char trans,
                          magma_int_t m, magma_int_t n,
                          cuDoubleComplex *a,    magma_int_t lda,
                          cuDoubleComplex *tau,
                          cuDoubleComplex *c,    magma_int_t ldc,
                          cuDoubleComplex *work, magma_int_t lwork,
                          magma_int_t *info);
magma_int_t magma_zunghr( magma_int_t n, magma_int_t ilo, magma_int_t ihi,
                          cuDoubleComplex *a, magma_int_t lda,
                          cuDoubleComplex *tau,
                          cuDoubleComplex *dT, magma_int_t nb,
                          magma_int_t *info);
magma_int_t magma_zheev( char jobz, char uplo, magma_int_t n,
                         cuDoubleComplex *a, magma_int_t lda, double *w,
                         cuDoubleComplex *work, magma_int_t lwork,
                         double *rwork, magma_int_t *info);
magma_int_t magma_zheevx(char jobz, char range, char uplo, magma_int_t n,
                         cuDoubleComplex *a, magma_int_t lda, double vl, double vu,
                         magma_int_t il, magma_int_t iu, double abstol, magma_int_t *m,
                         double *w, cuDoubleComplex *z, magma_int_t ldz, 
                         cuDoubleComplex *work, magma_int_t lwork,
                         double *rwork, magma_int_t *iwork, 
                         magma_int_t *ifail, magma_int_t *info);
#if defined(PRECISION_z) || defined(PRECISION_c)
magma_int_t  magma_zgeev( char jobvl, char jobvr, magma_int_t n,
                          cuDoubleComplex *a, magma_int_t lda,
                          cuDoubleComplex *w,
                          cuDoubleComplex *vl, magma_int_t ldvl,
                          cuDoubleComplex *vr, magma_int_t ldvr,
                          cuDoubleComplex *work, magma_int_t lwork,
                          double *rwork, magma_int_t *info);
magma_int_t magma_zgeqp3( magma_int_t m, magma_int_t n,
                          cuDoubleComplex *a, magma_int_t lda, 
                          magma_int_t *jpvt, cuDoubleComplex *tau,
                          cuDoubleComplex *work, magma_int_t lwork, 
                          double *rwork, magma_int_t *info);
magma_int_t magma_zgesvd( char jobu, char jobvt, magma_int_t m, magma_int_t n,
                          cuDoubleComplex *a,    magma_int_t lda, double *s, 
                          cuDoubleComplex *u,    magma_int_t ldu, 
                          cuDoubleComplex *vt,   magma_int_t ldvt,
                          cuDoubleComplex *work, magma_int_t lwork,
                          double *rwork, magma_int_t *info );
magma_int_t magma_zheevd( char jobz, char uplo, magma_int_t n,
                          cuDoubleComplex *a, magma_int_t lda, double *w,
                          cuDoubleComplex *work, magma_int_t lwork,
                          double *rwork, magma_int_t lrwork,
                          magma_int_t *iwork, magma_int_t liwork, magma_int_t *info);
magma_int_t magma_zheevr( char jobz, char range, char uplo, magma_int_t n,
                          cuDoubleComplex *a, magma_int_t lda, double vl, double vu,
                          magma_int_t il, magma_int_t iu, double abstol, magma_int_t *m,
                          double *w, cuDoubleComplex *z, magma_int_t ldz, 
                          magma_int_t *isuppz,
                          cuDoubleComplex *work, magma_int_t lwork,
                          double *rwork, magma_int_t lrwork, magma_int_t *iwork,
                          magma_int_t liwork, magma_int_t *info);
magma_int_t magma_zhegvd( magma_int_t itype, char jobz, char uplo, magma_int_t n,
                          cuDoubleComplex *a, magma_int_t lda,
                          cuDoubleComplex *b, magma_int_t ldb,
                          double *w, cuDoubleComplex *work, magma_int_t lwork,
                          double *rwork, magma_int_t lrwork, magma_int_t *iwork,
                          magma_int_t liwork, magma_int_t *info);
magma_int_t magma_zhegvdx(magma_int_t itype, char jobz, char range, char uplo, 
                          magma_int_t n, cuDoubleComplex *a, magma_int_t lda,
                          cuDoubleComplex *b, magma_int_t ldb,
                          double vl, double vu, magma_int_t il, magma_int_t iu,
                          magma_int_t *m, double *w, cuDoubleComplex *work, 
                          magma_int_t lwork, double *rwork,
                          magma_int_t lrwork, magma_int_t *iwork,
                          magma_int_t liwork, magma_int_t *info);
magma_int_t magma_zhegvx( magma_int_t itype, char jobz, char range, char uplo, 
                          magma_int_t n, cuDoubleComplex *a, magma_int_t lda, 
                          cuDoubleComplex *b, magma_int_t ldb,
                          double vl, double vu, magma_int_t il, magma_int_t iu,
                          double abstol, magma_int_t *m, double *w, 
                          cuDoubleComplex *z, magma_int_t ldz,
                          cuDoubleComplex *work, magma_int_t lwork, double *rwork,
                          magma_int_t *iwork, magma_int_t *ifail, magma_int_t *info);
magma_int_t magma_zhegvr( magma_int_t itype, char jobz, char range, char uplo, 
                          magma_int_t n, cuDoubleComplex *a, magma_int_t lda,
                          cuDoubleComplex *b, magma_int_t ldb,
                          double vl, double vu, magma_int_t il, magma_int_t iu,
                          double abstol, magma_int_t *m, double *w, 
                          cuDoubleComplex *z, magma_int_t ldz,
                          magma_int_t *isuppz, cuDoubleComplex *work, magma_int_t lwork,
                          double *rwork, magma_int_t lrwork, magma_int_t *iwork,
                          magma_int_t liwork, magma_int_t *info);
magma_int_t magma_zstedx( char range, magma_int_t n, double vl, double vu,
                          magma_int_t il, magma_int_t iu, double *D, double *E,
                          cuDoubleComplex *Z, magma_int_t ldz,
                          double *rwork, magma_int_t ldrwork, magma_int_t *iwork,
                          magma_int_t liwork, double* dwork, magma_int_t *info);
#else
magma_int_t  magma_zgeev( char jobvl, char jobvr, magma_int_t n,
                          cuDoubleComplex *a,    magma_int_t lda,
                          cuDoubleComplex *wr, cuDoubleComplex *wi,
                          cuDoubleComplex *vl,   magma_int_t ldvl,
                          cuDoubleComplex *vr,   magma_int_t ldvr,
                          cuDoubleComplex *work, magma_int_t lwork,
                          magma_int_t *info);
magma_int_t magma_zgeqp3( magma_int_t m, magma_int_t n,
                          cuDoubleComplex *a, magma_int_t lda, 
                          magma_int_t *jpvt, cuDoubleComplex *tau,
                          cuDoubleComplex *work, magma_int_t lwork,
                          magma_int_t *info);
magma_int_t magma_zgesvd( char jobu, char jobvt, magma_int_t m, magma_int_t n,
                          cuDoubleComplex *a,    magma_int_t lda, double *s, 
                          cuDoubleComplex *u,    magma_int_t ldu, 
                          cuDoubleComplex *vt,   magma_int_t ldvt,
                          cuDoubleComplex *work, magma_int_t lwork,
                          magma_int_t *info );
magma_int_t magma_zheevd( char jobz, char uplo, magma_int_t n,
                          cuDoubleComplex *a, magma_int_t lda, double *w,
                          cuDoubleComplex *work, magma_int_t lwork,
                          magma_int_t *iwork, magma_int_t liwork, magma_int_t *info);
magma_int_t magma_zhegvd( magma_int_t itype, char jobz, char uplo, magma_int_t n,
                          cuDoubleComplex *a, magma_int_t lda,
                          cuDoubleComplex *b, magma_int_t ldb,
                          double *w, cuDoubleComplex *work, magma_int_t lwork,
                          magma_int_t *iwork, magma_int_t liwork, magma_int_t *info);
magma_int_t magma_zstedx( char range, magma_int_t n, double vl, double vu,
                          magma_int_t il, magma_int_t iu, double* d, double* e,
                          double* z, magma_int_t ldz,
                          double* work, magma_int_t lwork,
                          magma_int_t* iwork, magma_int_t liwork,
                          double* dwork, magma_int_t* info);
magma_int_t magma_zlaex0( magma_int_t n, double* d, double* e, double* q, magma_int_t ldq,
                          double* work, magma_int_t* iwork, double* dwork,
                          char range, double vl, double vu,
                          magma_int_t il, magma_int_t iu, magma_int_t* info);
magma_int_t magma_zlaex1( magma_int_t n, double* d, double* q, magma_int_t ldq,
                          magma_int_t* indxq, double rho, magma_int_t cutpnt,
                          double* work, magma_int_t* iwork, double* dwork,
                          char range, double vl, double vu,
                          magma_int_t il, magma_int_t iu, magma_int_t* info);
magma_int_t magma_zlaex3( magma_int_t k, magma_int_t n, magma_int_t n1, double* d,
                          double* q, magma_int_t ldq, double rho,
                          double* dlamda, double* q2, magma_int_t* indx,
                          magma_int_t* ctot, double* w, double* s, magma_int_t* indxq,
                          double* dwork,
                          char range, double vl, double vu, magma_int_t il, magma_int_t iu,
                          magma_int_t* info );
#endif

magma_int_t magma_zhegst( magma_int_t itype, char uplo, magma_int_t n,
                          cuDoubleComplex *a, magma_int_t lda,
                          cuDoubleComplex *b, magma_int_t ldb, magma_int_t *info);

/* //////////////////////////////////////////////////////////////////////////// 
 -- MAGMA function definitions / Data on GPU
*/
magma_int_t magma_zgels_gpu(  char trans, magma_int_t m, magma_int_t n, magma_int_t nrhs,
                              cuDoubleComplex *dA,    magma_int_t ldda, 
                              cuDoubleComplex *dB,    magma_int_t lddb, 
                              cuDoubleComplex *hwork, magma_int_t lwork, 
                              magma_int_t *info);
magma_int_t magma_zgels3_gpu( char trans, magma_int_t m, magma_int_t n, magma_int_t nrhs,
                              cuDoubleComplex *dA,    magma_int_t ldda,
                              cuDoubleComplex *dB,    magma_int_t lddb,
                              cuDoubleComplex *hwork, magma_int_t lwork,
                              magma_int_t *info);
magma_int_t magma_zgelqf_gpu( magma_int_t m, magma_int_t n,
                              cuDoubleComplex *dA,    magma_int_t ldda,   cuDoubleComplex *tau,
                              cuDoubleComplex *work, magma_int_t lwork, magma_int_t *info);
magma_int_t magma_zgeqrf_gpu( magma_int_t m, magma_int_t n, 
                              cuDoubleComplex *dA,  magma_int_t ldda, 
                              cuDoubleComplex *tau, cuDoubleComplex *dT, 
                              magma_int_t *info);
magma_int_t magma_zgeqrf2_gpu(magma_int_t m, magma_int_t n, 
                              cuDoubleComplex *dA,  magma_int_t ldda, 
                              cuDoubleComplex *tau, magma_int_t *info);
magma_int_t magma_zgeqrf2_mgpu(magma_int_t num_gpus, magma_int_t m, magma_int_t n,
                               cuDoubleComplex **dlA, magma_int_t ldda,
                               cuDoubleComplex *tau, magma_int_t *info );
magma_int_t magma_zgeqrf3_gpu(magma_int_t m, magma_int_t n, 
                              cuDoubleComplex *dA,  magma_int_t ldda, 
                              cuDoubleComplex *tau, cuDoubleComplex *dT, 
                              magma_int_t *info);
magma_int_t magma_zgeqrs_gpu( magma_int_t m, magma_int_t n, magma_int_t nrhs, 
                              cuDoubleComplex *dA,     magma_int_t ldda, 
                              cuDoubleComplex *tau,   cuDoubleComplex *dT,
                              cuDoubleComplex *dB,    magma_int_t lddb,
                              cuDoubleComplex *hwork, magma_int_t lhwork, 
                              magma_int_t *info);
magma_int_t magma_zgeqrs3_gpu( magma_int_t m, magma_int_t n, magma_int_t nrhs, 
                              cuDoubleComplex *dA,     magma_int_t ldda, 
                              cuDoubleComplex *tau,   cuDoubleComplex *dT,
                              cuDoubleComplex *dB,    magma_int_t lddb,
                              cuDoubleComplex *hwork, magma_int_t lhwork, 
                              magma_int_t *info);
magma_int_t magma_zgessm_gpu( char storev, magma_int_t m, magma_int_t n, magma_int_t k, magma_int_t ib, 
                              magma_int_t *ipiv, 
                              cuDoubleComplex *dL1, magma_int_t lddl1, 
                              cuDoubleComplex *dL,  magma_int_t lddl, 
                              cuDoubleComplex *dA,  magma_int_t ldda, 
                              magma_int_t *info);
magma_int_t magma_zgesv_gpu(  magma_int_t n, magma_int_t nrhs, 
                              cuDoubleComplex *dA, magma_int_t ldda, magma_int_t *ipiv, 
                              cuDoubleComplex *dB, magma_int_t lddb, magma_int_t *info);
magma_int_t magma_zgetrf_incpiv_gpu( char storev, magma_int_t m, magma_int_t n, magma_int_t ib,
                              cuDoubleComplex *hA, magma_int_t ldha, cuDoubleComplex *dA, magma_int_t ldda,
                              cuDoubleComplex *hL, magma_int_t ldhl, cuDoubleComplex *dL, magma_int_t lddl,
                              magma_int_t *ipiv, 
                              cuDoubleComplex *dwork, magma_int_t lddwork,
                              magma_int_t *info);
magma_int_t magma_zgetrf_gpu( magma_int_t m, magma_int_t n, 
                              cuDoubleComplex *dA, magma_int_t ldda, 
                              magma_int_t *ipiv, magma_int_t *info);
magma_int_t 
magma_zgetrf_nopiv_gpu      ( magma_int_t m, magma_int_t n,
                              cuDoubleComplex *dA, magma_int_t ldda,
                              magma_int_t *info);
magma_int_t magma_zgetri_gpu( magma_int_t n, 
                              cuDoubleComplex *dA, magma_int_t ldda, magma_int_t *ipiv, 
                              cuDoubleComplex *dwork, magma_int_t lwork, magma_int_t *info);
magma_int_t magma_zgetrs_gpu( char trans, magma_int_t n, magma_int_t nrhs, 
                              cuDoubleComplex *dA, magma_int_t ldda, magma_int_t *ipiv, 
                              cuDoubleComplex *dB, magma_int_t lddb, magma_int_t *info);
magma_int_t magma_zlabrd_gpu( magma_int_t m, magma_int_t n, magma_int_t nb, 
                              cuDoubleComplex *a, magma_int_t lda, cuDoubleComplex *da, magma_int_t ldda,
                              double *d, double *e, cuDoubleComplex *tauq, cuDoubleComplex *taup,  
                              cuDoubleComplex *x, magma_int_t ldx, cuDoubleComplex *dx, magma_int_t lddx, 
                              cuDoubleComplex *y, magma_int_t ldy, cuDoubleComplex *dy, magma_int_t lddy);
magma_int_t magma_zlarfb_gpu( char side, char trans, char direct, char storev, 
                              magma_int_t m, magma_int_t n, magma_int_t k,
                              const cuDoubleComplex *dv, magma_int_t ldv,
                              const cuDoubleComplex *dt, magma_int_t ldt, 
                              cuDoubleComplex *dc,       magma_int_t ldc,
                              cuDoubleComplex *dwork,    magma_int_t ldwork );
magma_int_t magma_zposv_gpu(  char uplo, magma_int_t n, magma_int_t nrhs, 
                              cuDoubleComplex *dA, magma_int_t ldda, 
                              cuDoubleComplex *dB, magma_int_t lddb, magma_int_t *info);
magma_int_t magma_zpotrf_gpu( char uplo,  magma_int_t n, 
                              cuDoubleComplex *dA, magma_int_t ldda, magma_int_t *info);
magma_int_t magma_zpotri_gpu( char uplo,  magma_int_t n,
                              cuDoubleComplex *dA, magma_int_t ldda, magma_int_t *info);
magma_int_t magma_zlauum_gpu( char uplo,  magma_int_t n,
                              cuDoubleComplex *dA, magma_int_t ldda, magma_int_t *info);
magma_int_t magma_ztrtri_gpu( char uplo,  char diag, magma_int_t n,
                              cuDoubleComplex *dA, magma_int_t ldda, magma_int_t *info);
magma_int_t magma_zhetrd_gpu( char uplo, magma_int_t n,
                              cuDoubleComplex *da, magma_int_t ldda,
                              double *d, double *e, cuDoubleComplex *tau,
                              cuDoubleComplex *wa,  magma_int_t ldwa,
                              cuDoubleComplex *work, magma_int_t lwork,
                              magma_int_t *info);
magma_int_t magma_zhetrd2_gpu(char uplo, magma_int_t n,
                              cuDoubleComplex *da, magma_int_t ldda,
                              double *d, double *e, cuDoubleComplex *tau,
                              cuDoubleComplex *wa,  magma_int_t ldwa,
                              cuDoubleComplex *work, magma_int_t lwork,
                              cuDoubleComplex *dwork, magma_int_t ldwork,
                              magma_int_t *info);
magma_int_t magma_zpotrs_gpu( char uplo,  magma_int_t n, magma_int_t nrhs, 
                              cuDoubleComplex *dA, magma_int_t ldda, 
                              cuDoubleComplex *dB, magma_int_t lddb, magma_int_t *info);
magma_int_t magma_zssssm_gpu( char storev, magma_int_t m1, magma_int_t n1, 
                              magma_int_t m2, magma_int_t n2, magma_int_t k, magma_int_t ib, 
                              cuDoubleComplex *dA1, magma_int_t ldda1, 
                              cuDoubleComplex *dA2, magma_int_t ldda2, 
                              cuDoubleComplex *dL1, magma_int_t lddl1, 
                              cuDoubleComplex *dL2, magma_int_t lddl2,
                              magma_int_t *IPIV, magma_int_t *info);
magma_int_t magma_ztstrf_gpu( char storev, magma_int_t m, magma_int_t n, magma_int_t ib, magma_int_t nb,
                              cuDoubleComplex *hU, magma_int_t ldhu, cuDoubleComplex *dU, magma_int_t lddu, 
                              cuDoubleComplex *hA, magma_int_t ldha, cuDoubleComplex *dA, magma_int_t ldda, 
                              cuDoubleComplex *hL, magma_int_t ldhl, cuDoubleComplex *dL, magma_int_t lddl,
                              magma_int_t *ipiv, 
                              cuDoubleComplex *hwork, magma_int_t ldhwork, 
                              cuDoubleComplex *dwork, magma_int_t lddwork,
                              magma_int_t *info);
magma_int_t magma_zungqr_gpu( magma_int_t m, magma_int_t n, magma_int_t k, 
                              cuDoubleComplex *da, magma_int_t ldda, 
                              cuDoubleComplex *tau, cuDoubleComplex *dwork, 
                              magma_int_t nb, magma_int_t *info );
magma_int_t magma_zunmql2_gpu(const char side, const char trans,
                              magma_int_t m, magma_int_t n, magma_int_t k,
                              cuDoubleComplex *da, magma_int_t ldda,
                              cuDoubleComplex *tau,
                              cuDoubleComplex *dc, magma_int_t lddc,
                              cuDoubleComplex *wa, magma_int_t ldwa,
                              magma_int_t *info);
magma_int_t magma_zunmqr_gpu( char side, char trans, 
                              magma_int_t m, magma_int_t n, magma_int_t k,
                              cuDoubleComplex *a,    magma_int_t lda, cuDoubleComplex *tau, 
                              cuDoubleComplex *c,    magma_int_t ldc,
                              cuDoubleComplex *work, magma_int_t lwork, 
                              cuDoubleComplex *td,   magma_int_t nb, magma_int_t *info);
magma_int_t magma_zunmqr2_gpu(const char side, const char trans,
                              magma_int_t m, magma_int_t n, magma_int_t k,
                              cuDoubleComplex *da,   magma_int_t ldda,
                              cuDoubleComplex *tau,
                              cuDoubleComplex *dc,    magma_int_t lddc,
                              cuDoubleComplex *wa,    magma_int_t ldwa,
                              magma_int_t *info);
magma_int_t magma_zunmtr_gpu( char side, char uplo, char trans,
                              magma_int_t m, magma_int_t n,
                              cuDoubleComplex *da,    magma_int_t ldda,
                              cuDoubleComplex *tau,
                              cuDoubleComplex *dc,    magma_int_t lddc,
                              cuDoubleComplex *wa,    magma_int_t ldwa,
                              magma_int_t *info);

#if defined(PRECISION_z) || defined(PRECISION_c)
magma_int_t magma_zheevd_gpu( char jobz, char uplo,
                              magma_int_t n,
                              cuDoubleComplex *da, magma_int_t ldda,
                              double *w,
                              cuDoubleComplex *wa,  magma_int_t ldwa,
                              cuDoubleComplex *work, magma_int_t lwork,
                              double *rwork, magma_int_t lrwork,
                              magma_int_t *iwork, magma_int_t liwork,
                              magma_int_t *info);
magma_int_t magma_zheevdx_gpu(char jobz, char range, char uplo,
                              magma_int_t n, cuDoubleComplex *da, 
                              magma_int_t ldda, double vl, double vu, 
                              magma_int_t il, magma_int_t iu,
                              magma_int_t *m, double *w,
                              cuDoubleComplex *wa,  magma_int_t ldwa,
                              cuDoubleComplex *work, magma_int_t lwork,
                              double *rwork, magma_int_t lrwork,
                              magma_int_t *iwork, magma_int_t liwork,
                              magma_int_t *info);
magma_int_t magma_zheevr_gpu( char jobz, char range, char uplo, magma_int_t n,
                              cuDoubleComplex *da, magma_int_t ldda, double vl, double vu,
                              magma_int_t il, magma_int_t iu, double abstol, magma_int_t *m,
                              double *w, cuDoubleComplex *dz, magma_int_t lddz,
                              magma_int_t *isuppz,
                              cuDoubleComplex *wa, magma_int_t ldwa,
                              cuDoubleComplex *wz, magma_int_t ldwz,
                              cuDoubleComplex *work, magma_int_t lwork,
                              double *rwork, magma_int_t lrwork, magma_int_t *iwork,
                              magma_int_t liwork, magma_int_t *info);
#else
magma_int_t magma_zheevd_gpu( char jobz, char uplo,
                              magma_int_t n,
                              cuDoubleComplex *da, magma_int_t ldda,
                              cuDoubleComplex *w,
                              cuDoubleComplex *wa,  magma_int_t ldwa,
                              cuDoubleComplex *work, magma_int_t lwork,
                              magma_int_t *iwork, magma_int_t liwork,
                              magma_int_t *info);
#endif

magma_int_t magma_zheevx_gpu( char jobz, char range, char uplo, magma_int_t n,
                              cuDoubleComplex *da, magma_int_t ldda, double vl, 
                              double vu, magma_int_t il, magma_int_t iu, 
                              double abstol, magma_int_t *m,
                              double *w, cuDoubleComplex *dz, magma_int_t lddz,
                              cuDoubleComplex *wa, magma_int_t ldwa,
                              cuDoubleComplex *wz, magma_int_t ldwz,
                              cuDoubleComplex *work, magma_int_t lwork,
                              double *rwork, magma_int_t *iwork, 
                              magma_int_t *ifail, magma_int_t *info);
magma_int_t magma_zhegst_gpu(magma_int_t itype, char uplo, magma_int_t n,
                             cuDoubleComplex *da, magma_int_t ldda,
                             cuDoubleComplex *db, magma_int_t lddb, magma_int_t *info);


/* //////////////////////////////////////////////////////////////////////////// 
 -- MAGMA utility function definitions
*/

void magma_zprint    ( magma_int_t m, magma_int_t n, const cuDoubleComplex  *A, magma_int_t lda  );
void magma_zprint_gpu( magma_int_t m, magma_int_t n, const cuDoubleComplex *dA, magma_int_t ldda );

void zpanel_to_q( char uplo, magma_int_t ib, cuDoubleComplex *A, magma_int_t lda, cuDoubleComplex *work );
void zq_to_panel( char uplo, magma_int_t ib, cuDoubleComplex *A, magma_int_t lda, cuDoubleComplex *work );

#ifdef __cplusplus
}
#endif

#undef PRECISION_z
#endif /* _MAGMA_Z_H_ */
