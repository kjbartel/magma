/*
 *   -- MAGMA (version 1.3.0) --
 *      Univ. of Tennessee, Knoxville
 *      Univ. of California, Berkeley
 *      Univ. of Colorado, Denver
 *      November 2012
 *
 * @generated c Wed Nov 14 22:52:26 2012
 */

#ifndef _MAGMA_C_H_
#define _MAGMA_C_H_
#define PRECISION_c

#ifdef __cplusplus
extern "C" {
#endif

/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA Auxiliary functions to get the NB used
*/
int magma_get_cpotrf_nb(int m);
int magma_get_cgetrf_nb(int m);
int magma_get_cgetri_nb(int m);
int magma_get_cgeqp3_nb(int m);
int magma_get_cgeqrf_nb(int m);
int magma_get_cgeqlf_nb(int m);
int magma_get_cgehrd_nb(int m);
int magma_get_chetrd_nb(int m);
int magma_get_cgelqf_nb(int m);
int magma_get_cgebrd_nb(int m);
int magma_get_chegst_nb(int m);
int magma_get_cgesvd_nb(int m);

/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA function definitions / Data on CPU
*/
magma_int_t magma_cgebrd( magma_int_t m, magma_int_t n, cuFloatComplex *A, 
                          magma_int_t lda, float *d, float *e,
                          cuFloatComplex *tauq,  cuFloatComplex *taup, 
                          cuFloatComplex *work, magma_int_t lwork, magma_int_t *info);
magma_int_t magma_cgehrd2(magma_int_t n, magma_int_t ilo, magma_int_t ihi,
                          cuFloatComplex *A, magma_int_t lda, cuFloatComplex *tau, 
                          cuFloatComplex *work, magma_int_t *lwork, magma_int_t *info);
magma_int_t magma_cgehrd( magma_int_t n, magma_int_t ilo, magma_int_t ihi,
                          cuFloatComplex *A, magma_int_t lda, cuFloatComplex *tau,
                          cuFloatComplex *work, magma_int_t lwork,
                          cuFloatComplex *d_T, magma_int_t *info);
magma_int_t magma_cgelqf( magma_int_t m, magma_int_t n, 
                          cuFloatComplex *A,    magma_int_t lda,   cuFloatComplex *tau, 
                          cuFloatComplex *work, magma_int_t lwork, magma_int_t *info);
magma_int_t magma_cgeqlf( magma_int_t m, magma_int_t n, 
                          cuFloatComplex *A,    magma_int_t lda,   cuFloatComplex *tau, 
                          cuFloatComplex *work, magma_int_t lwork, magma_int_t *info);
magma_int_t magma_cgeqrf( magma_int_t m, magma_int_t n, cuFloatComplex *A, 
                          magma_int_t lda, cuFloatComplex *tau, cuFloatComplex *work, 
                          magma_int_t lwork, magma_int_t *info);
magma_int_t magma_cgeqr2( magma_int_t *m, magma_int_t *n, cuFloatComplex *a,
                          magma_int_t *lda, cuFloatComplex *tau, cuFloatComplex *work, 
                          magma_int_t *info);
magma_int_t magma_cgeqrf4(magma_int_t num_gpus, magma_int_t m, magma_int_t n,
                          cuFloatComplex *a,    magma_int_t lda, cuFloatComplex *tau,
                          cuFloatComplex *work, magma_int_t lwork, magma_int_t *info );
magma_int_t magma_cgeqrf_ooc( magma_int_t m, magma_int_t n, cuFloatComplex *A,
                          magma_int_t lda, cuFloatComplex *tau, cuFloatComplex *work,
                          magma_int_t lwork, magma_int_t *info);
magma_int_t magma_cgesv ( magma_int_t n, magma_int_t nrhs, 
                          cuFloatComplex *A, magma_int_t lda, magma_int_t *ipiv, 
                          cuFloatComplex *B, magma_int_t ldb, magma_int_t *info);
magma_int_t magma_cgetrf( magma_int_t m, magma_int_t n, cuFloatComplex *A, 
                          magma_int_t lda, magma_int_t *ipiv, 
                          magma_int_t *info);
magma_int_t magma_cgetrf2(magma_int_t m, magma_int_t n, cuFloatComplex *a, 
                          magma_int_t lda, magma_int_t *ipiv, magma_int_t *info);

magma_int_t magma_claqps( magma_int_t m, magma_int_t n, magma_int_t offset, 
                          magma_int_t nb, magma_int_t *kb, 
                          cuFloatComplex *A,  magma_int_t lda,
                          cuFloatComplex *dA, magma_int_t ldda,
                          magma_int_t *jpvt, cuFloatComplex *tau, float *vn1, float *vn2, 
                          cuFloatComplex *auxv, 
                          cuFloatComplex *F,  magma_int_t ldf,
                          cuFloatComplex *dF, magma_int_t lddf );
magma_int_t magma_clarf(  char *side, magma_int_t *m, magma_int_t *n,
                          cuFloatComplex *v, magma_int_t *incv, cuFloatComplex *tau,
                          cuFloatComplex *c__, magma_int_t *ldc, cuFloatComplex *work);
magma_int_t magma_clarfg( magma_int_t *n, cuFloatComplex *alpha, cuFloatComplex *x,
                          magma_int_t *incx, cuFloatComplex *tau);
magma_int_t magma_clatrd( char uplo, magma_int_t n, magma_int_t nb, cuFloatComplex *a, 
                          magma_int_t lda, float *e, cuFloatComplex *tau, 
                          cuFloatComplex *w, magma_int_t ldw,
                          cuFloatComplex *da, magma_int_t ldda, 
                          cuFloatComplex *dw, magma_int_t lddw);
magma_int_t magma_clatrd2(char uplo, magma_int_t n, magma_int_t nb,
                          cuFloatComplex *a,  magma_int_t lda,
                          float *e, cuFloatComplex *tau,
                          cuFloatComplex *w,  magma_int_t ldw,
                          cuFloatComplex *da, magma_int_t ldda,
                          cuFloatComplex *dw, magma_int_t lddw,
                          cuFloatComplex *dwork, magma_int_t ldwork);
magma_int_t magma_clahr2( magma_int_t m, magma_int_t n, magma_int_t nb, 
                          cuFloatComplex *da, cuFloatComplex *dv, cuFloatComplex *a, 
                          magma_int_t lda, cuFloatComplex *tau, cuFloatComplex *t, 
                          magma_int_t ldt, cuFloatComplex *y, magma_int_t ldy);
magma_int_t magma_clahru( magma_int_t n, magma_int_t ihi, magma_int_t k, magma_int_t nb, 
                          cuFloatComplex *a, magma_int_t lda, 
                          cuFloatComplex *da, cuFloatComplex *y, 
                          cuFloatComplex *v, cuFloatComplex *t, 
                          cuFloatComplex *dwork);
magma_int_t magma_cposv ( char uplo, magma_int_t n, magma_int_t nrhs, 
                          cuFloatComplex *A, magma_int_t lda, 
                          cuFloatComplex *B, magma_int_t ldb, magma_int_t *info);
magma_int_t magma_cpotrf( char uplo, magma_int_t n, cuFloatComplex *A, 
                          magma_int_t lda, magma_int_t *info);
magma_int_t magma_cpotri( char uplo, magma_int_t n, cuFloatComplex *A,
                          magma_int_t lda, magma_int_t *info);
magma_int_t magma_clauum( char uplo, magma_int_t n, cuFloatComplex *A,
                          magma_int_t lda, magma_int_t *info);
magma_int_t magma_ctrtri( char uplo, char diag, magma_int_t n, cuFloatComplex *A, 
                          magma_int_t lda, magma_int_t *info);
magma_int_t magma_chetrd( char uplo, magma_int_t n, cuFloatComplex *A, 
                          magma_int_t lda, float *d, float *e, 
                          cuFloatComplex *tau, cuFloatComplex *work, magma_int_t lwork, 
                          magma_int_t *info);
magma_int_t magma_cungqr( magma_int_t m, magma_int_t n, magma_int_t k,
                          cuFloatComplex *a, magma_int_t lda,
                          cuFloatComplex *tau, cuFloatComplex *dwork,
                          magma_int_t nb, magma_int_t *info );
magma_int_t magma_cunmql( const char side, const char trans,
                          magma_int_t m, magma_int_t n, magma_int_t k,
                          cuFloatComplex *a, magma_int_t lda,
                          cuFloatComplex *tau,
                          cuFloatComplex *c, magma_int_t ldc,
                          cuFloatComplex *work, magma_int_t lwork,
                          magma_int_t *info);
magma_int_t magma_cunmqr( char side, char trans, 
                          magma_int_t m, magma_int_t n, magma_int_t k, 
                          cuFloatComplex *a, magma_int_t lda, cuFloatComplex *tau, 
                          cuFloatComplex *c, magma_int_t ldc, 
                          cuFloatComplex *work, magma_int_t lwork, magma_int_t *info);
magma_int_t magma_cunmtr( char side, char uplo, char trans,
                          magma_int_t m, magma_int_t n,
                          cuFloatComplex *a,    magma_int_t lda,
                          cuFloatComplex *tau,
                          cuFloatComplex *c,    magma_int_t ldc,
                          cuFloatComplex *work, magma_int_t lwork,
                          magma_int_t *info);
magma_int_t magma_cunghr( magma_int_t n, magma_int_t ilo, magma_int_t ihi,
                          cuFloatComplex *a, magma_int_t lda,
                          cuFloatComplex *tau,
                          cuFloatComplex *dT, magma_int_t nb,
                          magma_int_t *info);
magma_int_t magma_cheev( char jobz, char uplo, magma_int_t n,
                         cuFloatComplex *a, magma_int_t lda, float *w,
                         cuFloatComplex *work, magma_int_t lwork,
                         float *rwork, magma_int_t *info);
magma_int_t magma_cheevx(char jobz, char range, char uplo, magma_int_t n,
                         cuFloatComplex *a, magma_int_t lda, float vl, float vu,
                         magma_int_t il, magma_int_t iu, float abstol, magma_int_t *m,
                         float *w, cuFloatComplex *z, magma_int_t ldz, 
                         cuFloatComplex *work, magma_int_t lwork,
                         float *rwork, magma_int_t *iwork, 
                         magma_int_t *ifail, magma_int_t *info);
#if defined(PRECISION_z) || defined(PRECISION_c)
magma_int_t  magma_cgeev( char jobvl, char jobvr, magma_int_t n,
                          cuFloatComplex *a, magma_int_t lda,
                          cuFloatComplex *w,
                          cuFloatComplex *vl, magma_int_t ldvl,
                          cuFloatComplex *vr, magma_int_t ldvr,
                          cuFloatComplex *work, magma_int_t lwork,
                          float *rwork, magma_int_t *info);
magma_int_t magma_cgeqp3( magma_int_t m, magma_int_t n,
                          cuFloatComplex *a, magma_int_t lda, 
                          magma_int_t *jpvt, cuFloatComplex *tau,
                          cuFloatComplex *work, magma_int_t lwork, 
                          float *rwork, magma_int_t *info);
magma_int_t magma_cgesvd( char jobu, char jobvt, magma_int_t m, magma_int_t n,
                          cuFloatComplex *a,    magma_int_t lda, float *s, 
                          cuFloatComplex *u,    magma_int_t ldu, 
                          cuFloatComplex *vt,   magma_int_t ldvt,
                          cuFloatComplex *work, magma_int_t lwork,
                          float *rwork, magma_int_t *info );
magma_int_t magma_cheevd( char jobz, char uplo, magma_int_t n,
                          cuFloatComplex *a, magma_int_t lda, float *w,
                          cuFloatComplex *work, magma_int_t lwork,
                          float *rwork, magma_int_t lrwork,
                          magma_int_t *iwork, magma_int_t liwork, magma_int_t *info);
magma_int_t magma_cheevr( char jobz, char range, char uplo, magma_int_t n,
                          cuFloatComplex *a, magma_int_t lda, float vl, float vu,
                          magma_int_t il, magma_int_t iu, float abstol, magma_int_t *m,
                          float *w, cuFloatComplex *z, magma_int_t ldz, 
                          magma_int_t *isuppz,
                          cuFloatComplex *work, magma_int_t lwork,
                          float *rwork, magma_int_t lrwork, magma_int_t *iwork,
                          magma_int_t liwork, magma_int_t *info);
magma_int_t magma_chegvd( magma_int_t itype, char jobz, char uplo, magma_int_t n,
                          cuFloatComplex *a, magma_int_t lda,
                          cuFloatComplex *b, magma_int_t ldb,
                          float *w, cuFloatComplex *work, magma_int_t lwork,
                          float *rwork, magma_int_t lrwork, magma_int_t *iwork,
                          magma_int_t liwork, magma_int_t *info);
magma_int_t magma_chegvdx(magma_int_t itype, char jobz, char range, char uplo, 
                          magma_int_t n, cuFloatComplex *a, magma_int_t lda,
                          cuFloatComplex *b, magma_int_t ldb,
                          float vl, float vu, magma_int_t il, magma_int_t iu,
                          magma_int_t *m, float *w, cuFloatComplex *work, 
                          magma_int_t lwork, float *rwork,
                          magma_int_t lrwork, magma_int_t *iwork,
                          magma_int_t liwork, magma_int_t *info);
magma_int_t magma_chegvx( magma_int_t itype, char jobz, char range, char uplo, 
                          magma_int_t n, cuFloatComplex *a, magma_int_t lda, 
                          cuFloatComplex *b, magma_int_t ldb,
                          float vl, float vu, magma_int_t il, magma_int_t iu,
                          float abstol, magma_int_t *m, float *w, 
                          cuFloatComplex *z, magma_int_t ldz,
                          cuFloatComplex *work, magma_int_t lwork, float *rwork,
                          magma_int_t *iwork, magma_int_t *ifail, magma_int_t *info);
magma_int_t magma_chegvr( magma_int_t itype, char jobz, char range, char uplo, 
                          magma_int_t n, cuFloatComplex *a, magma_int_t lda,
                          cuFloatComplex *b, magma_int_t ldb,
                          float vl, float vu, magma_int_t il, magma_int_t iu,
                          float abstol, magma_int_t *m, float *w, 
                          cuFloatComplex *z, magma_int_t ldz,
                          magma_int_t *isuppz, cuFloatComplex *work, magma_int_t lwork,
                          float *rwork, magma_int_t lrwork, magma_int_t *iwork,
                          magma_int_t liwork, magma_int_t *info);
magma_int_t magma_cstedx( char range, magma_int_t n, float vl, float vu,
                          magma_int_t il, magma_int_t iu, float *D, float *E,
                          cuFloatComplex *Z, magma_int_t ldz,
                          float *rwork, magma_int_t ldrwork, magma_int_t *iwork,
                          magma_int_t liwork, float* dwork, magma_int_t *info);
#else
magma_int_t  magma_cgeev( char jobvl, char jobvr, magma_int_t n,
                          cuFloatComplex *a,    magma_int_t lda,
                          cuFloatComplex *wr, cuFloatComplex *wi,
                          cuFloatComplex *vl,   magma_int_t ldvl,
                          cuFloatComplex *vr,   magma_int_t ldvr,
                          cuFloatComplex *work, magma_int_t lwork,
                          magma_int_t *info);
magma_int_t magma_cgeqp3( magma_int_t m, magma_int_t n,
                          cuFloatComplex *a, magma_int_t lda, 
                          magma_int_t *jpvt, cuFloatComplex *tau,
                          cuFloatComplex *work, magma_int_t lwork,
                          magma_int_t *info);
magma_int_t magma_cgesvd( char jobu, char jobvt, magma_int_t m, magma_int_t n,
                          cuFloatComplex *a,    magma_int_t lda, float *s, 
                          cuFloatComplex *u,    magma_int_t ldu, 
                          cuFloatComplex *vt,   magma_int_t ldvt,
                          cuFloatComplex *work, magma_int_t lwork,
                          magma_int_t *info );
magma_int_t magma_cheevd( char jobz, char uplo, magma_int_t n,
                          cuFloatComplex *a, magma_int_t lda, float *w,
                          cuFloatComplex *work, magma_int_t lwork,
                          magma_int_t *iwork, magma_int_t liwork, magma_int_t *info);
magma_int_t magma_chegvd( magma_int_t itype, char jobz, char uplo, magma_int_t n,
                          cuFloatComplex *a, magma_int_t lda,
                          cuFloatComplex *b, magma_int_t ldb,
                          float *w, cuFloatComplex *work, magma_int_t lwork,
                          magma_int_t *iwork, magma_int_t liwork, magma_int_t *info);
magma_int_t magma_cstedx( char range, magma_int_t n, float vl, float vu,
                          magma_int_t il, magma_int_t iu, float* d, float* e,
                          float* z, magma_int_t ldz,
                          float* work, magma_int_t lwork,
                          magma_int_t* iwork, magma_int_t liwork,
                          float* dwork, magma_int_t* info);
magma_int_t magma_claex0( magma_int_t n, float* d, float* e, float* q, magma_int_t ldq,
                          float* work, magma_int_t* iwork, float* dwork,
                          char range, float vl, float vu,
                          magma_int_t il, magma_int_t iu, magma_int_t* info);
magma_int_t magma_claex1( magma_int_t n, float* d, float* q, magma_int_t ldq,
                          magma_int_t* indxq, float rho, magma_int_t cutpnt,
                          float* work, magma_int_t* iwork, float* dwork,
                          char range, float vl, float vu,
                          magma_int_t il, magma_int_t iu, magma_int_t* info);
magma_int_t magma_claex3( magma_int_t k, magma_int_t n, magma_int_t n1, float* d,
                          float* q, magma_int_t ldq, float rho,
                          float* dlamda, float* q2, magma_int_t* indx,
                          magma_int_t* ctot, float* w, float* s, magma_int_t* indxq,
                          float* dwork,
                          char range, float vl, float vu, magma_int_t il, magma_int_t iu,
                          magma_int_t* info );
#endif

magma_int_t magma_chegst( magma_int_t itype, char uplo, magma_int_t n,
                          cuFloatComplex *a, magma_int_t lda,
                          cuFloatComplex *b, magma_int_t ldb, magma_int_t *info);

/* //////////////////////////////////////////////////////////////////////////// 
 -- MAGMA function definitions / Data on GPU
*/
magma_int_t magma_cgels_gpu(  char trans, magma_int_t m, magma_int_t n, magma_int_t nrhs,
                              cuFloatComplex *dA,    magma_int_t ldda, 
                              cuFloatComplex *dB,    magma_int_t lddb, 
                              cuFloatComplex *hwork, magma_int_t lwork, 
                              magma_int_t *info);
magma_int_t magma_cgels3_gpu( char trans, magma_int_t m, magma_int_t n, magma_int_t nrhs,
                              cuFloatComplex *dA,    magma_int_t ldda,
                              cuFloatComplex *dB,    magma_int_t lddb,
                              cuFloatComplex *hwork, magma_int_t lwork,
                              magma_int_t *info);
magma_int_t magma_cgelqf_gpu( magma_int_t m, magma_int_t n,
                              cuFloatComplex *dA,    magma_int_t ldda,   cuFloatComplex *tau,
                              cuFloatComplex *work, magma_int_t lwork, magma_int_t *info);
magma_int_t magma_cgeqrf_gpu( magma_int_t m, magma_int_t n, 
                              cuFloatComplex *dA,  magma_int_t ldda, 
                              cuFloatComplex *tau, cuFloatComplex *dT, 
                              magma_int_t *info);
magma_int_t magma_cgeqrf2_gpu(magma_int_t m, magma_int_t n, 
                              cuFloatComplex *dA,  magma_int_t ldda, 
                              cuFloatComplex *tau, magma_int_t *info);
magma_int_t magma_cgeqrf2_mgpu(magma_int_t num_gpus, magma_int_t m, magma_int_t n,
                               cuFloatComplex **dlA, magma_int_t ldda,
                               cuFloatComplex *tau, magma_int_t *info );
magma_int_t magma_cgeqrf3_gpu(magma_int_t m, magma_int_t n, 
                              cuFloatComplex *dA,  magma_int_t ldda, 
                              cuFloatComplex *tau, cuFloatComplex *dT, 
                              magma_int_t *info);
magma_int_t magma_cgeqrs_gpu( magma_int_t m, magma_int_t n, magma_int_t nrhs, 
                              cuFloatComplex *dA,     magma_int_t ldda, 
                              cuFloatComplex *tau,   cuFloatComplex *dT,
                              cuFloatComplex *dB,    magma_int_t lddb,
                              cuFloatComplex *hwork, magma_int_t lhwork, 
                              magma_int_t *info);
magma_int_t magma_cgeqrs3_gpu( magma_int_t m, magma_int_t n, magma_int_t nrhs, 
                              cuFloatComplex *dA,     magma_int_t ldda, 
                              cuFloatComplex *tau,   cuFloatComplex *dT,
                              cuFloatComplex *dB,    magma_int_t lddb,
                              cuFloatComplex *hwork, magma_int_t lhwork, 
                              magma_int_t *info);
magma_int_t magma_cgessm_gpu( char storev, magma_int_t m, magma_int_t n, magma_int_t k, magma_int_t ib, 
                              magma_int_t *ipiv, 
                              cuFloatComplex *dL1, magma_int_t lddl1, 
                              cuFloatComplex *dL,  magma_int_t lddl, 
                              cuFloatComplex *dA,  magma_int_t ldda, 
                              magma_int_t *info);
magma_int_t magma_cgesv_gpu(  magma_int_t n, magma_int_t nrhs, 
                              cuFloatComplex *dA, magma_int_t ldda, magma_int_t *ipiv, 
                              cuFloatComplex *dB, magma_int_t lddb, magma_int_t *info);
magma_int_t magma_cgetrf_incpiv_gpu( char storev, magma_int_t m, magma_int_t n, magma_int_t ib,
                              cuFloatComplex *hA, magma_int_t ldha, cuFloatComplex *dA, magma_int_t ldda,
                              cuFloatComplex *hL, magma_int_t ldhl, cuFloatComplex *dL, magma_int_t lddl,
                              magma_int_t *ipiv, 
                              cuFloatComplex *dwork, magma_int_t lddwork,
                              magma_int_t *info);
magma_int_t magma_cgetrf_gpu( magma_int_t m, magma_int_t n, 
                              cuFloatComplex *dA, magma_int_t ldda, 
                              magma_int_t *ipiv, magma_int_t *info);
magma_int_t 
magma_cgetrf_nopiv_gpu      ( magma_int_t m, magma_int_t n,
                              cuFloatComplex *dA, magma_int_t ldda,
                              magma_int_t *info);
magma_int_t magma_cgetri_gpu( magma_int_t n, 
                              cuFloatComplex *dA, magma_int_t ldda, magma_int_t *ipiv, 
                              cuFloatComplex *dwork, magma_int_t lwork, magma_int_t *info);
magma_int_t magma_cgetrs_gpu( char trans, magma_int_t n, magma_int_t nrhs, 
                              cuFloatComplex *dA, magma_int_t ldda, magma_int_t *ipiv, 
                              cuFloatComplex *dB, magma_int_t lddb, magma_int_t *info);
magma_int_t magma_clabrd_gpu( magma_int_t m, magma_int_t n, magma_int_t nb, 
                              cuFloatComplex *a, magma_int_t lda, cuFloatComplex *da, magma_int_t ldda,
                              float *d, float *e, cuFloatComplex *tauq, cuFloatComplex *taup,  
                              cuFloatComplex *x, magma_int_t ldx, cuFloatComplex *dx, magma_int_t lddx, 
                              cuFloatComplex *y, magma_int_t ldy, cuFloatComplex *dy, magma_int_t lddy);
magma_int_t magma_clarfb_gpu( char side, char trans, char direct, char storev, 
                              magma_int_t m, magma_int_t n, magma_int_t k,
                              const cuFloatComplex *dv, magma_int_t ldv,
                              const cuFloatComplex *dt, magma_int_t ldt, 
                              cuFloatComplex *dc,       magma_int_t ldc,
                              cuFloatComplex *dwork,    magma_int_t ldwork );
magma_int_t magma_cposv_gpu(  char uplo, magma_int_t n, magma_int_t nrhs, 
                              cuFloatComplex *dA, magma_int_t ldda, 
                              cuFloatComplex *dB, magma_int_t lddb, magma_int_t *info);
magma_int_t magma_cpotrf_gpu( char uplo,  magma_int_t n, 
                              cuFloatComplex *dA, magma_int_t ldda, magma_int_t *info);
magma_int_t magma_cpotri_gpu( char uplo,  magma_int_t n,
                              cuFloatComplex *dA, magma_int_t ldda, magma_int_t *info);
magma_int_t magma_clauum_gpu( char uplo,  magma_int_t n,
                              cuFloatComplex *dA, magma_int_t ldda, magma_int_t *info);
magma_int_t magma_ctrtri_gpu( char uplo,  char diag, magma_int_t n,
                              cuFloatComplex *dA, magma_int_t ldda, magma_int_t *info);
magma_int_t magma_chetrd_gpu( char uplo, magma_int_t n,
                              cuFloatComplex *da, magma_int_t ldda,
                              float *d, float *e, cuFloatComplex *tau,
                              cuFloatComplex *wa,  magma_int_t ldwa,
                              cuFloatComplex *work, magma_int_t lwork,
                              magma_int_t *info);
magma_int_t magma_chetrd2_gpu(char uplo, magma_int_t n,
                              cuFloatComplex *da, magma_int_t ldda,
                              float *d, float *e, cuFloatComplex *tau,
                              cuFloatComplex *wa,  magma_int_t ldwa,
                              cuFloatComplex *work, magma_int_t lwork,
                              cuFloatComplex *dwork, magma_int_t ldwork,
                              magma_int_t *info);
magma_int_t magma_cpotrs_gpu( char uplo,  magma_int_t n, magma_int_t nrhs, 
                              cuFloatComplex *dA, magma_int_t ldda, 
                              cuFloatComplex *dB, magma_int_t lddb, magma_int_t *info);
magma_int_t magma_cssssm_gpu( char storev, magma_int_t m1, magma_int_t n1, 
                              magma_int_t m2, magma_int_t n2, magma_int_t k, magma_int_t ib, 
                              cuFloatComplex *dA1, magma_int_t ldda1, 
                              cuFloatComplex *dA2, magma_int_t ldda2, 
                              cuFloatComplex *dL1, magma_int_t lddl1, 
                              cuFloatComplex *dL2, magma_int_t lddl2,
                              magma_int_t *IPIV, magma_int_t *info);
magma_int_t magma_ctstrf_gpu( char storev, magma_int_t m, magma_int_t n, magma_int_t ib, magma_int_t nb,
                              cuFloatComplex *hU, magma_int_t ldhu, cuFloatComplex *dU, magma_int_t lddu, 
                              cuFloatComplex *hA, magma_int_t ldha, cuFloatComplex *dA, magma_int_t ldda, 
                              cuFloatComplex *hL, magma_int_t ldhl, cuFloatComplex *dL, magma_int_t lddl,
                              magma_int_t *ipiv, 
                              cuFloatComplex *hwork, magma_int_t ldhwork, 
                              cuFloatComplex *dwork, magma_int_t lddwork,
                              magma_int_t *info);
magma_int_t magma_cungqr_gpu( magma_int_t m, magma_int_t n, magma_int_t k, 
                              cuFloatComplex *da, magma_int_t ldda, 
                              cuFloatComplex *tau, cuFloatComplex *dwork, 
                              magma_int_t nb, magma_int_t *info );
magma_int_t magma_cunmql2_gpu(const char side, const char trans,
                              magma_int_t m, magma_int_t n, magma_int_t k,
                              cuFloatComplex *da, magma_int_t ldda,
                              cuFloatComplex *tau,
                              cuFloatComplex *dc, magma_int_t lddc,
                              cuFloatComplex *wa, magma_int_t ldwa,
                              magma_int_t *info);
magma_int_t magma_cunmqr_gpu( char side, char trans, 
                              magma_int_t m, magma_int_t n, magma_int_t k,
                              cuFloatComplex *a,    magma_int_t lda, cuFloatComplex *tau, 
                              cuFloatComplex *c,    magma_int_t ldc,
                              cuFloatComplex *work, magma_int_t lwork, 
                              cuFloatComplex *td,   magma_int_t nb, magma_int_t *info);
magma_int_t magma_cunmqr2_gpu(const char side, const char trans,
                              magma_int_t m, magma_int_t n, magma_int_t k,
                              cuFloatComplex *da,   magma_int_t ldda,
                              cuFloatComplex *tau,
                              cuFloatComplex *dc,    magma_int_t lddc,
                              cuFloatComplex *wa,    magma_int_t ldwa,
                              magma_int_t *info);
magma_int_t magma_cunmtr_gpu( char side, char uplo, char trans,
                              magma_int_t m, magma_int_t n,
                              cuFloatComplex *da,    magma_int_t ldda,
                              cuFloatComplex *tau,
                              cuFloatComplex *dc,    magma_int_t lddc,
                              cuFloatComplex *wa,    magma_int_t ldwa,
                              magma_int_t *info);

#if defined(PRECISION_z) || defined(PRECISION_c)
magma_int_t magma_cheevd_gpu( char jobz, char uplo,
                              magma_int_t n,
                              cuFloatComplex *da, magma_int_t ldda,
                              float *w,
                              cuFloatComplex *wa,  magma_int_t ldwa,
                              cuFloatComplex *work, magma_int_t lwork,
                              float *rwork, magma_int_t lrwork,
                              magma_int_t *iwork, magma_int_t liwork,
                              magma_int_t *info);
magma_int_t magma_cheevdx_gpu(char jobz, char range, char uplo,
                              magma_int_t n, cuFloatComplex *da, 
                              magma_int_t ldda, float vl, float vu, 
                              magma_int_t il, magma_int_t iu,
                              magma_int_t *m, float *w,
                              cuFloatComplex *wa,  magma_int_t ldwa,
                              cuFloatComplex *work, magma_int_t lwork,
                              float *rwork, magma_int_t lrwork,
                              magma_int_t *iwork, magma_int_t liwork,
                              magma_int_t *info);
magma_int_t magma_cheevr_gpu( char jobz, char range, char uplo, magma_int_t n,
                              cuFloatComplex *da, magma_int_t ldda, float vl, float vu,
                              magma_int_t il, magma_int_t iu, float abstol, magma_int_t *m,
                              float *w, cuFloatComplex *dz, magma_int_t lddz,
                              magma_int_t *isuppz,
                              cuFloatComplex *wa, magma_int_t ldwa,
                              cuFloatComplex *wz, magma_int_t ldwz,
                              cuFloatComplex *work, magma_int_t lwork,
                              float *rwork, magma_int_t lrwork, magma_int_t *iwork,
                              magma_int_t liwork, magma_int_t *info);
#else
magma_int_t magma_cheevd_gpu( char jobz, char uplo,
                              magma_int_t n,
                              cuFloatComplex *da, magma_int_t ldda,
                              cuFloatComplex *w,
                              cuFloatComplex *wa,  magma_int_t ldwa,
                              cuFloatComplex *work, magma_int_t lwork,
                              magma_int_t *iwork, magma_int_t liwork,
                              magma_int_t *info);
#endif

magma_int_t magma_cheevx_gpu( char jobz, char range, char uplo, magma_int_t n,
                              cuFloatComplex *da, magma_int_t ldda, float vl, 
                              float vu, magma_int_t il, magma_int_t iu, 
                              float abstol, magma_int_t *m,
                              float *w, cuFloatComplex *dz, magma_int_t lddz,
                              cuFloatComplex *wa, magma_int_t ldwa,
                              cuFloatComplex *wz, magma_int_t ldwz,
                              cuFloatComplex *work, magma_int_t lwork,
                              float *rwork, magma_int_t *iwork, 
                              magma_int_t *ifail, magma_int_t *info);
magma_int_t magma_chegst_gpu(magma_int_t itype, char uplo, magma_int_t n,
                             cuFloatComplex *da, magma_int_t ldda,
                             cuFloatComplex *db, magma_int_t lddb, magma_int_t *info);


/* //////////////////////////////////////////////////////////////////////////// 
 -- MAGMA utility function definitions
*/

void magma_cprint    ( magma_int_t m, magma_int_t n, const cuFloatComplex  *A, magma_int_t lda  );
void magma_cprint_gpu( magma_int_t m, magma_int_t n, const cuFloatComplex *dA, magma_int_t ldda );

void cpanel_to_q( char uplo, magma_int_t ib, cuFloatComplex *A, magma_int_t lda, cuFloatComplex *work );
void cq_to_panel( char uplo, magma_int_t ib, cuFloatComplex *A, magma_int_t lda, cuFloatComplex *work );

#ifdef __cplusplus
}
#endif

#undef PRECISION_c
#endif /* _MAGMA_C_H_ */
