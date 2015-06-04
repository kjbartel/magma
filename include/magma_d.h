/*
 *   -- MAGMA (version 1.3.0) --
 *      Univ. of Tennessee, Knoxville
 *      Univ. of California, Berkeley
 *      Univ. of Colorado, Denver
 *      November 2012
 *
 * @generated d Wed Nov 14 22:52:26 2012
 */

#ifndef _MAGMA_D_H_
#define _MAGMA_D_H_
#define PRECISION_d

#ifdef __cplusplus
extern "C" {
#endif

/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA Auxiliary functions to get the NB used
*/
int magma_get_dpotrf_nb(int m);
int magma_get_dgetrf_nb(int m);
int magma_get_dgetri_nb(int m);
int magma_get_dgeqp3_nb(int m);
int magma_get_dgeqrf_nb(int m);
int magma_get_dgeqlf_nb(int m);
int magma_get_dgehrd_nb(int m);
int magma_get_dsytrd_nb(int m);
int magma_get_dgelqf_nb(int m);
int magma_get_dgebrd_nb(int m);
int magma_get_dsygst_nb(int m);
int magma_get_dgesvd_nb(int m);

/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA function definitions / Data on CPU
*/
magma_int_t magma_dgebrd( magma_int_t m, magma_int_t n, double *A, 
                          magma_int_t lda, double *d, double *e,
                          double *tauq,  double *taup, 
                          double *work, magma_int_t lwork, magma_int_t *info);
magma_int_t magma_dgehrd2(magma_int_t n, magma_int_t ilo, magma_int_t ihi,
                          double *A, magma_int_t lda, double *tau, 
                          double *work, magma_int_t *lwork, magma_int_t *info);
magma_int_t magma_dgehrd( magma_int_t n, magma_int_t ilo, magma_int_t ihi,
                          double *A, magma_int_t lda, double *tau,
                          double *work, magma_int_t lwork,
                          double *d_T, magma_int_t *info);
magma_int_t magma_dgelqf( magma_int_t m, magma_int_t n, 
                          double *A,    magma_int_t lda,   double *tau, 
                          double *work, magma_int_t lwork, magma_int_t *info);
magma_int_t magma_dgeqlf( magma_int_t m, magma_int_t n, 
                          double *A,    magma_int_t lda,   double *tau, 
                          double *work, magma_int_t lwork, magma_int_t *info);
magma_int_t magma_dgeqrf( magma_int_t m, magma_int_t n, double *A, 
                          magma_int_t lda, double *tau, double *work, 
                          magma_int_t lwork, magma_int_t *info);
magma_int_t magma_dgeqr2( magma_int_t *m, magma_int_t *n, double *a,
                          magma_int_t *lda, double *tau, double *work, 
                          magma_int_t *info);
magma_int_t magma_dgeqrf4(magma_int_t num_gpus, magma_int_t m, magma_int_t n,
                          double *a,    magma_int_t lda, double *tau,
                          double *work, magma_int_t lwork, magma_int_t *info );
magma_int_t magma_dgeqrf_ooc( magma_int_t m, magma_int_t n, double *A,
                          magma_int_t lda, double *tau, double *work,
                          magma_int_t lwork, magma_int_t *info);
magma_int_t magma_dgesv ( magma_int_t n, magma_int_t nrhs, 
                          double *A, magma_int_t lda, magma_int_t *ipiv, 
                          double *B, magma_int_t ldb, magma_int_t *info);
magma_int_t magma_dgetrf( magma_int_t m, magma_int_t n, double *A, 
                          magma_int_t lda, magma_int_t *ipiv, 
                          magma_int_t *info);
magma_int_t magma_dgetrf2(magma_int_t m, magma_int_t n, double *a, 
                          magma_int_t lda, magma_int_t *ipiv, magma_int_t *info);

magma_int_t magma_dlaqps( magma_int_t m, magma_int_t n, magma_int_t offset, 
                          magma_int_t nb, magma_int_t *kb, 
                          double *A,  magma_int_t lda,
                          double *dA, magma_int_t ldda,
                          magma_int_t *jpvt, double *tau, double *vn1, double *vn2, 
                          double *auxv, 
                          double *F,  magma_int_t ldf,
                          double *dF, magma_int_t lddf );
magma_int_t magma_dlarf(  char *side, magma_int_t *m, magma_int_t *n,
                          double *v, magma_int_t *incv, double *tau,
                          double *c__, magma_int_t *ldc, double *work);
magma_int_t magma_dlarfg( magma_int_t *n, double *alpha, double *x,
                          magma_int_t *incx, double *tau);
magma_int_t magma_dlatrd( char uplo, magma_int_t n, magma_int_t nb, double *a, 
                          magma_int_t lda, double *e, double *tau, 
                          double *w, magma_int_t ldw,
                          double *da, magma_int_t ldda, 
                          double *dw, magma_int_t lddw);
magma_int_t magma_dlatrd2(char uplo, magma_int_t n, magma_int_t nb,
                          double *a,  magma_int_t lda,
                          double *e, double *tau,
                          double *w,  magma_int_t ldw,
                          double *da, magma_int_t ldda,
                          double *dw, magma_int_t lddw,
                          double *dwork, magma_int_t ldwork);
magma_int_t magma_dlahr2( magma_int_t m, magma_int_t n, magma_int_t nb, 
                          double *da, double *dv, double *a, 
                          magma_int_t lda, double *tau, double *t, 
                          magma_int_t ldt, double *y, magma_int_t ldy);
magma_int_t magma_dlahru( magma_int_t n, magma_int_t ihi, magma_int_t k, magma_int_t nb, 
                          double *a, magma_int_t lda, 
                          double *da, double *y, 
                          double *v, double *t, 
                          double *dwork);
magma_int_t magma_dposv ( char uplo, magma_int_t n, magma_int_t nrhs, 
                          double *A, magma_int_t lda, 
                          double *B, magma_int_t ldb, magma_int_t *info);
magma_int_t magma_dpotrf( char uplo, magma_int_t n, double *A, 
                          magma_int_t lda, magma_int_t *info);
magma_int_t magma_dpotri( char uplo, magma_int_t n, double *A,
                          magma_int_t lda, magma_int_t *info);
magma_int_t magma_dlauum( char uplo, magma_int_t n, double *A,
                          magma_int_t lda, magma_int_t *info);
magma_int_t magma_dtrtri( char uplo, char diag, magma_int_t n, double *A, 
                          magma_int_t lda, magma_int_t *info);
magma_int_t magma_dsytrd( char uplo, magma_int_t n, double *A, 
                          magma_int_t lda, double *d, double *e, 
                          double *tau, double *work, magma_int_t lwork, 
                          magma_int_t *info);
magma_int_t magma_dorgqr( magma_int_t m, magma_int_t n, magma_int_t k,
                          double *a, magma_int_t lda,
                          double *tau, double *dwork,
                          magma_int_t nb, magma_int_t *info );
magma_int_t magma_dormql( const char side, const char trans,
                          magma_int_t m, magma_int_t n, magma_int_t k,
                          double *a, magma_int_t lda,
                          double *tau,
                          double *c, magma_int_t ldc,
                          double *work, magma_int_t lwork,
                          magma_int_t *info);
magma_int_t magma_dormqr( char side, char trans, 
                          magma_int_t m, magma_int_t n, magma_int_t k, 
                          double *a, magma_int_t lda, double *tau, 
                          double *c, magma_int_t ldc, 
                          double *work, magma_int_t lwork, magma_int_t *info);
magma_int_t magma_dormtr( char side, char uplo, char trans,
                          magma_int_t m, magma_int_t n,
                          double *a,    magma_int_t lda,
                          double *tau,
                          double *c,    magma_int_t ldc,
                          double *work, magma_int_t lwork,
                          magma_int_t *info);
magma_int_t magma_dorghr( magma_int_t n, magma_int_t ilo, magma_int_t ihi,
                          double *a, magma_int_t lda,
                          double *tau,
                          double *dT, magma_int_t nb,
                          magma_int_t *info);
magma_int_t magma_dsyev( char jobz, char uplo, magma_int_t n,
                         double *a, magma_int_t lda, double *w,
                         double *work, magma_int_t lwork,
                         double *rwork, magma_int_t *info);
magma_int_t magma_dsyevx(char jobz, char range, char uplo, magma_int_t n,
                         double *a, magma_int_t lda, double vl, double vu,
                         magma_int_t il, magma_int_t iu, double abstol, magma_int_t *m,
                         double *w, double *z, magma_int_t ldz, 
                         double *work, magma_int_t lwork,
                         double *rwork, magma_int_t *iwork, 
                         magma_int_t *ifail, magma_int_t *info);
#if defined(PRECISION_z) || defined(PRECISION_c)
magma_int_t  magma_dgeev( char jobvl, char jobvr, magma_int_t n,
                          double *a, magma_int_t lda,
                          double *w,
                          double *vl, magma_int_t ldvl,
                          double *vr, magma_int_t ldvr,
                          double *work, magma_int_t lwork,
                          double *rwork, magma_int_t *info);
magma_int_t magma_dgeqp3( magma_int_t m, magma_int_t n,
                          double *a, magma_int_t lda, 
                          magma_int_t *jpvt, double *tau,
                          double *work, magma_int_t lwork, 
                          double *rwork, magma_int_t *info);
magma_int_t magma_dgesvd( char jobu, char jobvt, magma_int_t m, magma_int_t n,
                          double *a,    magma_int_t lda, double *s, 
                          double *u,    magma_int_t ldu, 
                          double *vt,   magma_int_t ldvt,
                          double *work, magma_int_t lwork,
                          double *rwork, magma_int_t *info );
magma_int_t magma_dsyevd( char jobz, char uplo, magma_int_t n,
                          double *a, magma_int_t lda, double *w,
                          double *work, magma_int_t lwork,
                          double *rwork, magma_int_t lrwork,
                          magma_int_t *iwork, magma_int_t liwork, magma_int_t *info);
magma_int_t magma_dsyevr( char jobz, char range, char uplo, magma_int_t n,
                          double *a, magma_int_t lda, double vl, double vu,
                          magma_int_t il, magma_int_t iu, double abstol, magma_int_t *m,
                          double *w, double *z, magma_int_t ldz, 
                          magma_int_t *isuppz,
                          double *work, magma_int_t lwork,
                          double *rwork, magma_int_t lrwork, magma_int_t *iwork,
                          magma_int_t liwork, magma_int_t *info);
magma_int_t magma_dsygvd( magma_int_t itype, char jobz, char uplo, magma_int_t n,
                          double *a, magma_int_t lda,
                          double *b, magma_int_t ldb,
                          double *w, double *work, magma_int_t lwork,
                          double *rwork, magma_int_t lrwork, magma_int_t *iwork,
                          magma_int_t liwork, magma_int_t *info);
magma_int_t magma_dsygvdx(magma_int_t itype, char jobz, char range, char uplo, 
                          magma_int_t n, double *a, magma_int_t lda,
                          double *b, magma_int_t ldb,
                          double vl, double vu, magma_int_t il, magma_int_t iu,
                          magma_int_t *m, double *w, double *work, 
                          magma_int_t lwork, double *rwork,
                          magma_int_t lrwork, magma_int_t *iwork,
                          magma_int_t liwork, magma_int_t *info);
magma_int_t magma_dsygvx( magma_int_t itype, char jobz, char range, char uplo, 
                          magma_int_t n, double *a, magma_int_t lda, 
                          double *b, magma_int_t ldb,
                          double vl, double vu, magma_int_t il, magma_int_t iu,
                          double abstol, magma_int_t *m, double *w, 
                          double *z, magma_int_t ldz,
                          double *work, magma_int_t lwork, double *rwork,
                          magma_int_t *iwork, magma_int_t *ifail, magma_int_t *info);
magma_int_t magma_dsygvr( magma_int_t itype, char jobz, char range, char uplo, 
                          magma_int_t n, double *a, magma_int_t lda,
                          double *b, magma_int_t ldb,
                          double vl, double vu, magma_int_t il, magma_int_t iu,
                          double abstol, magma_int_t *m, double *w, 
                          double *z, magma_int_t ldz,
                          magma_int_t *isuppz, double *work, magma_int_t lwork,
                          double *rwork, magma_int_t lrwork, magma_int_t *iwork,
                          magma_int_t liwork, magma_int_t *info);
magma_int_t magma_dstedx( char range, magma_int_t n, double vl, double vu,
                          magma_int_t il, magma_int_t iu, double *D, double *E,
                          double *Z, magma_int_t ldz,
                          double *rwork, magma_int_t ldrwork, magma_int_t *iwork,
                          magma_int_t liwork, double* dwork, magma_int_t *info);
#else
magma_int_t  magma_dgeev( char jobvl, char jobvr, magma_int_t n,
                          double *a,    magma_int_t lda,
                          double *wr, double *wi,
                          double *vl,   magma_int_t ldvl,
                          double *vr,   magma_int_t ldvr,
                          double *work, magma_int_t lwork,
                          magma_int_t *info);
magma_int_t magma_dgeqp3( magma_int_t m, magma_int_t n,
                          double *a, magma_int_t lda, 
                          magma_int_t *jpvt, double *tau,
                          double *work, magma_int_t lwork,
                          magma_int_t *info);
magma_int_t magma_dgesvd( char jobu, char jobvt, magma_int_t m, magma_int_t n,
                          double *a,    magma_int_t lda, double *s, 
                          double *u,    magma_int_t ldu, 
                          double *vt,   magma_int_t ldvt,
                          double *work, magma_int_t lwork,
                          magma_int_t *info );
magma_int_t magma_dsyevd( char jobz, char uplo, magma_int_t n,
                          double *a, magma_int_t lda, double *w,
                          double *work, magma_int_t lwork,
                          magma_int_t *iwork, magma_int_t liwork, magma_int_t *info);
magma_int_t magma_dsygvd( magma_int_t itype, char jobz, char uplo, magma_int_t n,
                          double *a, magma_int_t lda,
                          double *b, magma_int_t ldb,
                          double *w, double *work, magma_int_t lwork,
                          magma_int_t *iwork, magma_int_t liwork, magma_int_t *info);
magma_int_t magma_dstedx( char range, magma_int_t n, double vl, double vu,
                          magma_int_t il, magma_int_t iu, double* d, double* e,
                          double* z, magma_int_t ldz,
                          double* work, magma_int_t lwork,
                          magma_int_t* iwork, magma_int_t liwork,
                          double* dwork, magma_int_t* info);
magma_int_t magma_dlaex0( magma_int_t n, double* d, double* e, double* q, magma_int_t ldq,
                          double* work, magma_int_t* iwork, double* dwork,
                          char range, double vl, double vu,
                          magma_int_t il, magma_int_t iu, magma_int_t* info);
magma_int_t magma_dlaex1( magma_int_t n, double* d, double* q, magma_int_t ldq,
                          magma_int_t* indxq, double rho, magma_int_t cutpnt,
                          double* work, magma_int_t* iwork, double* dwork,
                          char range, double vl, double vu,
                          magma_int_t il, magma_int_t iu, magma_int_t* info);
magma_int_t magma_dlaex3( magma_int_t k, magma_int_t n, magma_int_t n1, double* d,
                          double* q, magma_int_t ldq, double rho,
                          double* dlamda, double* q2, magma_int_t* indx,
                          magma_int_t* ctot, double* w, double* s, magma_int_t* indxq,
                          double* dwork,
                          char range, double vl, double vu, magma_int_t il, magma_int_t iu,
                          magma_int_t* info );
#endif

magma_int_t magma_dsygst( magma_int_t itype, char uplo, magma_int_t n,
                          double *a, magma_int_t lda,
                          double *b, magma_int_t ldb, magma_int_t *info);

/* //////////////////////////////////////////////////////////////////////////// 
 -- MAGMA function definitions / Data on GPU
*/
magma_int_t magma_dgels_gpu(  char trans, magma_int_t m, magma_int_t n, magma_int_t nrhs,
                              double *dA,    magma_int_t ldda, 
                              double *dB,    magma_int_t lddb, 
                              double *hwork, magma_int_t lwork, 
                              magma_int_t *info);
magma_int_t magma_dgels3_gpu( char trans, magma_int_t m, magma_int_t n, magma_int_t nrhs,
                              double *dA,    magma_int_t ldda,
                              double *dB,    magma_int_t lddb,
                              double *hwork, magma_int_t lwork,
                              magma_int_t *info);
magma_int_t magma_dgelqf_gpu( magma_int_t m, magma_int_t n,
                              double *dA,    magma_int_t ldda,   double *tau,
                              double *work, magma_int_t lwork, magma_int_t *info);
magma_int_t magma_dgeqrf_gpu( magma_int_t m, magma_int_t n, 
                              double *dA,  magma_int_t ldda, 
                              double *tau, double *dT, 
                              magma_int_t *info);
magma_int_t magma_dgeqrf2_gpu(magma_int_t m, magma_int_t n, 
                              double *dA,  magma_int_t ldda, 
                              double *tau, magma_int_t *info);
magma_int_t magma_dgeqrf2_mgpu(magma_int_t num_gpus, magma_int_t m, magma_int_t n,
                               double **dlA, magma_int_t ldda,
                               double *tau, magma_int_t *info );
magma_int_t magma_dgeqrf3_gpu(magma_int_t m, magma_int_t n, 
                              double *dA,  magma_int_t ldda, 
                              double *tau, double *dT, 
                              magma_int_t *info);
magma_int_t magma_dgeqrs_gpu( magma_int_t m, magma_int_t n, magma_int_t nrhs, 
                              double *dA,     magma_int_t ldda, 
                              double *tau,   double *dT,
                              double *dB,    magma_int_t lddb,
                              double *hwork, magma_int_t lhwork, 
                              magma_int_t *info);
magma_int_t magma_dgeqrs3_gpu( magma_int_t m, magma_int_t n, magma_int_t nrhs, 
                              double *dA,     magma_int_t ldda, 
                              double *tau,   double *dT,
                              double *dB,    magma_int_t lddb,
                              double *hwork, magma_int_t lhwork, 
                              magma_int_t *info);
magma_int_t magma_dgessm_gpu( char storev, magma_int_t m, magma_int_t n, magma_int_t k, magma_int_t ib, 
                              magma_int_t *ipiv, 
                              double *dL1, magma_int_t lddl1, 
                              double *dL,  magma_int_t lddl, 
                              double *dA,  magma_int_t ldda, 
                              magma_int_t *info);
magma_int_t magma_dgesv_gpu(  magma_int_t n, magma_int_t nrhs, 
                              double *dA, magma_int_t ldda, magma_int_t *ipiv, 
                              double *dB, magma_int_t lddb, magma_int_t *info);
magma_int_t magma_dgetrf_incpiv_gpu( char storev, magma_int_t m, magma_int_t n, magma_int_t ib,
                              double *hA, magma_int_t ldha, double *dA, magma_int_t ldda,
                              double *hL, magma_int_t ldhl, double *dL, magma_int_t lddl,
                              magma_int_t *ipiv, 
                              double *dwork, magma_int_t lddwork,
                              magma_int_t *info);
magma_int_t magma_dgetrf_gpu( magma_int_t m, magma_int_t n, 
                              double *dA, magma_int_t ldda, 
                              magma_int_t *ipiv, magma_int_t *info);
magma_int_t 
magma_dgetrf_nopiv_gpu      ( magma_int_t m, magma_int_t n,
                              double *dA, magma_int_t ldda,
                              magma_int_t *info);
magma_int_t magma_dgetri_gpu( magma_int_t n, 
                              double *dA, magma_int_t ldda, magma_int_t *ipiv, 
                              double *dwork, magma_int_t lwork, magma_int_t *info);
magma_int_t magma_dgetrs_gpu( char trans, magma_int_t n, magma_int_t nrhs, 
                              double *dA, magma_int_t ldda, magma_int_t *ipiv, 
                              double *dB, magma_int_t lddb, magma_int_t *info);
magma_int_t magma_dlabrd_gpu( magma_int_t m, magma_int_t n, magma_int_t nb, 
                              double *a, magma_int_t lda, double *da, magma_int_t ldda,
                              double *d, double *e, double *tauq, double *taup,  
                              double *x, magma_int_t ldx, double *dx, magma_int_t lddx, 
                              double *y, magma_int_t ldy, double *dy, magma_int_t lddy);
magma_int_t magma_dlarfb_gpu( char side, char trans, char direct, char storev, 
                              magma_int_t m, magma_int_t n, magma_int_t k,
                              const double *dv, magma_int_t ldv,
                              const double *dt, magma_int_t ldt, 
                              double *dc,       magma_int_t ldc,
                              double *dwork,    magma_int_t ldwork );
magma_int_t magma_dposv_gpu(  char uplo, magma_int_t n, magma_int_t nrhs, 
                              double *dA, magma_int_t ldda, 
                              double *dB, magma_int_t lddb, magma_int_t *info);
magma_int_t magma_dpotrf_gpu( char uplo,  magma_int_t n, 
                              double *dA, magma_int_t ldda, magma_int_t *info);
magma_int_t magma_dpotri_gpu( char uplo,  magma_int_t n,
                              double *dA, magma_int_t ldda, magma_int_t *info);
magma_int_t magma_dlauum_gpu( char uplo,  magma_int_t n,
                              double *dA, magma_int_t ldda, magma_int_t *info);
magma_int_t magma_dtrtri_gpu( char uplo,  char diag, magma_int_t n,
                              double *dA, magma_int_t ldda, magma_int_t *info);
magma_int_t magma_dsytrd_gpu( char uplo, magma_int_t n,
                              double *da, magma_int_t ldda,
                              double *d, double *e, double *tau,
                              double *wa,  magma_int_t ldwa,
                              double *work, magma_int_t lwork,
                              magma_int_t *info);
magma_int_t magma_dsytrd2_gpu(char uplo, magma_int_t n,
                              double *da, magma_int_t ldda,
                              double *d, double *e, double *tau,
                              double *wa,  magma_int_t ldwa,
                              double *work, magma_int_t lwork,
                              double *dwork, magma_int_t ldwork,
                              magma_int_t *info);
magma_int_t magma_dpotrs_gpu( char uplo,  magma_int_t n, magma_int_t nrhs, 
                              double *dA, magma_int_t ldda, 
                              double *dB, magma_int_t lddb, magma_int_t *info);
magma_int_t magma_dssssm_gpu( char storev, magma_int_t m1, magma_int_t n1, 
                              magma_int_t m2, magma_int_t n2, magma_int_t k, magma_int_t ib, 
                              double *dA1, magma_int_t ldda1, 
                              double *dA2, magma_int_t ldda2, 
                              double *dL1, magma_int_t lddl1, 
                              double *dL2, magma_int_t lddl2,
                              magma_int_t *IPIV, magma_int_t *info);
magma_int_t magma_dtstrf_gpu( char storev, magma_int_t m, magma_int_t n, magma_int_t ib, magma_int_t nb,
                              double *hU, magma_int_t ldhu, double *dU, magma_int_t lddu, 
                              double *hA, magma_int_t ldha, double *dA, magma_int_t ldda, 
                              double *hL, magma_int_t ldhl, double *dL, magma_int_t lddl,
                              magma_int_t *ipiv, 
                              double *hwork, magma_int_t ldhwork, 
                              double *dwork, magma_int_t lddwork,
                              magma_int_t *info);
magma_int_t magma_dorgqr_gpu( magma_int_t m, magma_int_t n, magma_int_t k, 
                              double *da, magma_int_t ldda, 
                              double *tau, double *dwork, 
                              magma_int_t nb, magma_int_t *info );
magma_int_t magma_dormql2_gpu(const char side, const char trans,
                              magma_int_t m, magma_int_t n, magma_int_t k,
                              double *da, magma_int_t ldda,
                              double *tau,
                              double *dc, magma_int_t lddc,
                              double *wa, magma_int_t ldwa,
                              magma_int_t *info);
magma_int_t magma_dormqr_gpu( char side, char trans, 
                              magma_int_t m, magma_int_t n, magma_int_t k,
                              double *a,    magma_int_t lda, double *tau, 
                              double *c,    magma_int_t ldc,
                              double *work, magma_int_t lwork, 
                              double *td,   magma_int_t nb, magma_int_t *info);
magma_int_t magma_dormqr2_gpu(const char side, const char trans,
                              magma_int_t m, magma_int_t n, magma_int_t k,
                              double *da,   magma_int_t ldda,
                              double *tau,
                              double *dc,    magma_int_t lddc,
                              double *wa,    magma_int_t ldwa,
                              magma_int_t *info);
magma_int_t magma_dormtr_gpu( char side, char uplo, char trans,
                              magma_int_t m, magma_int_t n,
                              double *da,    magma_int_t ldda,
                              double *tau,
                              double *dc,    magma_int_t lddc,
                              double *wa,    magma_int_t ldwa,
                              magma_int_t *info);

#if defined(PRECISION_z) || defined(PRECISION_c)
magma_int_t magma_dsyevd_gpu( char jobz, char uplo,
                              magma_int_t n,
                              double *da, magma_int_t ldda,
                              double *w,
                              double *wa,  magma_int_t ldwa,
                              double *work, magma_int_t lwork,
                              double *rwork, magma_int_t lrwork,
                              magma_int_t *iwork, magma_int_t liwork,
                              magma_int_t *info);
magma_int_t magma_dsyevdx_gpu(char jobz, char range, char uplo,
                              magma_int_t n, double *da, 
                              magma_int_t ldda, double vl, double vu, 
                              magma_int_t il, magma_int_t iu,
                              magma_int_t *m, double *w,
                              double *wa,  magma_int_t ldwa,
                              double *work, magma_int_t lwork,
                              double *rwork, magma_int_t lrwork,
                              magma_int_t *iwork, magma_int_t liwork,
                              magma_int_t *info);
magma_int_t magma_dsyevr_gpu( char jobz, char range, char uplo, magma_int_t n,
                              double *da, magma_int_t ldda, double vl, double vu,
                              magma_int_t il, magma_int_t iu, double abstol, magma_int_t *m,
                              double *w, double *dz, magma_int_t lddz,
                              magma_int_t *isuppz,
                              double *wa, magma_int_t ldwa,
                              double *wz, magma_int_t ldwz,
                              double *work, magma_int_t lwork,
                              double *rwork, magma_int_t lrwork, magma_int_t *iwork,
                              magma_int_t liwork, magma_int_t *info);
#else
magma_int_t magma_dsyevd_gpu( char jobz, char uplo,
                              magma_int_t n,
                              double *da, magma_int_t ldda,
                              double *w,
                              double *wa,  magma_int_t ldwa,
                              double *work, magma_int_t lwork,
                              magma_int_t *iwork, magma_int_t liwork,
                              magma_int_t *info);
#endif

magma_int_t magma_dsyevx_gpu( char jobz, char range, char uplo, magma_int_t n,
                              double *da, magma_int_t ldda, double vl, 
                              double vu, magma_int_t il, magma_int_t iu, 
                              double abstol, magma_int_t *m,
                              double *w, double *dz, magma_int_t lddz,
                              double *wa, magma_int_t ldwa,
                              double *wz, magma_int_t ldwz,
                              double *work, magma_int_t lwork,
                              double *rwork, magma_int_t *iwork, 
                              magma_int_t *ifail, magma_int_t *info);
magma_int_t magma_dsygst_gpu(magma_int_t itype, char uplo, magma_int_t n,
                             double *da, magma_int_t ldda,
                             double *db, magma_int_t lddb, magma_int_t *info);


/* //////////////////////////////////////////////////////////////////////////// 
 -- MAGMA utility function definitions
*/

void magma_dprint    ( magma_int_t m, magma_int_t n, const double  *A, magma_int_t lda  );
void magma_dprint_gpu( magma_int_t m, magma_int_t n, const double *dA, magma_int_t ldda );

void dpanel_to_q( char uplo, magma_int_t ib, double *A, magma_int_t lda, double *work );
void dq_to_panel( char uplo, magma_int_t ib, double *A, magma_int_t lda, double *work );

#ifdef __cplusplus
}
#endif

#undef PRECISION_d
#endif /* _MAGMA_D_H_ */
