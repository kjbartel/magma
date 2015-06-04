/*
 *   -- MAGMA (version 1.2.0) --
 *      Univ. of Tennessee, Knoxville
 *      Univ. of California, Berkeley
 *      Univ. of Colorado, Denver
 *      May 2012
 *
 * @generated d Tue May 15 18:17:06 2012
 */

#ifndef MAGMA_DLAPACK_H
#define MAGMA_DLAPACK_H

#define PRECISION_d
#ifdef __cplusplus
extern "C" {
#endif

/* ////////////////////////////////////////////////////////////////////////////
   -- LAPACK Externs used in MAGMA
*/

#define blasf77_daxpy      FORTRAN_NAME( daxpy,  ZAXPY  )
#define blasf77_dcopy      FORTRAN_NAME( dcopy,  ZCOPY  )

/* real versions use C wrapper to return value; no name mangling. */
#if  defined(PRECISION_z) || defined(PRECISION_c)    
#define blasf77_ddot      ddot
#else
#define blasf77_ddot      FORTRAN_NAME( ddot,  DDOT  )
#endif

#define blasf77_dgemm      FORTRAN_NAME( dgemm,  DGEMM  )
#define blasf77_dgemv      FORTRAN_NAME( dgemv,  ZGEMV  )
#define blasf77_dsymm      FORTRAN_NAME( dsymm,  DSYMM  )
#define blasf77_dsymv      FORTRAN_NAME( dsymv,  DSYMV  )
#define blasf77_dsyr2k     FORTRAN_NAME( dsyr2k, DSYR2K )
#define blasf77_dher2      FORTRAN_NAME( dsyr2,  DSYR2  )
#define blasf77_dsyrk      FORTRAN_NAME( dsyrk,  DSYRK  )
#define blasf77_dscal      FORTRAN_NAME( dscal,  ZSCAL  )
#define blasf77_ddscal     FORTRAN_NAME( dscal, ZDSCAL ) 
#define blasf77_dsymm      FORTRAN_NAME( dsymm,  DSYMM  )
#define blasf77_dsyr2k     FORTRAN_NAME( dsyr2k, DSYR2K )
#define blasf77_dsyrk      FORTRAN_NAME( dsyrk,  DSYRK  )
#define blasf77_dswap      FORTRAN_NAME( dswap,  ZSWAP  )
#define blasf77_dtrmm      FORTRAN_NAME( dtrmm,  DTRMM  )
#define blasf77_dtrmv      FORTRAN_NAME( dtrmv,  ZTRMV  )
#define blasf77_dtrsm      FORTRAN_NAME( dtrsm,  DTRSM  )
#define blasf77_dtrsv      FORTRAN_NAME( dtrsv,  DTRSV  )
#define blasf77_dger      FORTRAN_NAME( dger,  ZGERU  )

#define lapackf77_dbdsqr   FORTRAN_NAME( dbdsqr, ZBDSQR )
#define lapackf77_dgebak   FORTRAN_NAME( dgebak, ZGEBAK )
#define lapackf77_dgebal   FORTRAN_NAME( dgebal, ZGEBAL )
#define lapackf77_dgebd2   FORTRAN_NAME( dgebd2, ZGEBD2 )
#define lapackf77_dgebrd   FORTRAN_NAME( dgebrd, DGEBRD )
#define lapackf77_dgeev    FORTRAN_NAME( dgeev,  DGEEV  )
#define lapackf77_dgehd2   FORTRAN_NAME( dgehd2, ZGEHD2 )
#define lapackf77_dgehrd   FORTRAN_NAME( dgehrd, DGEHRD )
#define lapackf77_dgelqf   FORTRAN_NAME( dgelqf, DGELQF )
#define lapackf77_dgels    FORTRAN_NAME( dgels,  DGELS  )
#define lapackf77_dgeqlf   FORTRAN_NAME( dgeqlf, ZGEQLF )
#define lapackf77_dgeqp3   FORTRAN_NAME( dgeqp3, DGEQP3 )
#define lapackf77_dgeqrf   FORTRAN_NAME( dgeqrf, DGEQRF )
#define lapackf77_dgesvd   FORTRAN_NAME( dgesvd, DGESVD )
#define lapackf77_dgetrf   FORTRAN_NAME( dgetrf, DGETRF )
#define lapackf77_dgetri   FORTRAN_NAME( dgetri, DGETRI )
#define lapackf77_dgetrs   FORTRAN_NAME( dgetrs, DGETRS )
#define lapackf77_dsyev    FORTRAN_NAME( dsyev,  DSYEV  )
#define lapackf77_dsyevd   FORTRAN_NAME( dsyevd, DSYEVD )
#define lapackf77_dhegs2   FORTRAN_NAME( dsygs2, ZHEGS2 )
#define lapackf77_dhegvd   FORTRAN_NAME( dsygvd, DSYGVD )
#define lapackf77_dsytd2   FORTRAN_NAME( dsytd2, ZHETD2 )
#define lapackf77_dsytrd   FORTRAN_NAME( dsytrd, DSYTRD )
#define lapackf77_dsbtrd   FORTRAN_NAME( dsbtrd, DSBTRD )
#define lapackf77_dhseqr   FORTRAN_NAME( dhseqr, ZHSEQR )
#define lapackf77_dlacpy   FORTRAN_NAME( dlacpy, ZLACPY )
#define lapackf77_dlacgv   FORTRAN_NAME( dlacgv, ZLACGV )
#define lapackf77_dlange   FORTRAN_NAME( dlange, DLANGE )
#define lapackf77_dlansy   FORTRAN_NAME( dlansy, DLANSY )
#define lapackf77_dlansy   FORTRAN_NAME( dlansy, DLANSY )
#define lapackf77_dlarfb   FORTRAN_NAME( dlarfb, DLARFB )
#define lapackf77_dlarfg   FORTRAN_NAME( dlarfg, ZLARFG )
#define lapackf77_dlarft   FORTRAN_NAME( dlarft, ZLARFT )
#define lapackf77_dlarnv   FORTRAN_NAME( dlarnv, ZLARNV )
#define lapackf77_dlartg   FORTRAN_NAME( dlartg, ZLARTG )
#define lapackf77_dlascl   FORTRAN_NAME( dlascl, ZLASCL )
#define lapackf77_dlaset   FORTRAN_NAME( dlaset, ZLASET )
#define lapackf77_dlaswp   FORTRAN_NAME( dlaswp, ZLASWP )
#define lapackf77_dlatrd   FORTRAN_NAME( dlatrd, DLATRD )
#define lapackf77_dlabrd   FORTRAN_NAME( dlabrd, DLABRD )
#define lapackf77_dlauum   FORTRAN_NAME( dlauum, ZLAUUM )
#define lapackf77_dpotrf   FORTRAN_NAME( dpotrf, DPOTRF )
#define lapackf77_dpotrs   FORTRAN_NAME( dpotrs, DPOTRS )
#define lapackf77_dpotri   FORTRAN_NAME( dpotri, ZPOTRI )
#define lapackf77_dtrevc   FORTRAN_NAME( dtrevc, ZTREVC )
#define lapackf77_dstebz   FORTRAN_NAME( dstebz, DSTEBZ )
#define lapackf77_dlamc3   FORTRAN_NAME( dlamc3, DLAMC3 )
#define lapackf77_dlaed4   FORTRAN_NAME( dlaed4, DLAED4 )
#define lapackf77_dlamrg   FORTRAN_NAME( dlamrg, DLAMRG )
#define lapackf77_dtrtri   FORTRAN_NAME( dtrtri, ZTRTRI )
#define lapackf77_dsteqr   FORTRAN_NAME( dsteqr, ZSTEQR )
#define lapackf77_dstedc   FORTRAN_NAME( dstedc, ZSTEDC )
#define lapackf77_dstein   FORTRAN_NAME( dstein, ZSTEIN )
#define lapackf77_dstemr   FORTRAN_NAME( dstemr, ZSTEMR )
#define lapackf77_dsymv    FORTRAN_NAME( dsymv,  ZSYMV  )
#define lapackf77_dorg2r   FORTRAN_NAME( dorg2r, ZUNG2R )
#define lapackf77_dorgbr   FORTRAN_NAME( dorgbr, ZUNGBR )
#define lapackf77_dorghr   FORTRAN_NAME( dorghr, DORGHR )
#define lapackf77_dorglq   FORTRAN_NAME( dorglq, DORGLQ )
#define lapackf77_dungql   FORTRAN_NAME( dorgql, ZUNGQL )
#define lapackf77_dorgqr   FORTRAN_NAME( dorgqr, DORGQR )
#define lapackf77_dorgtr   FORTRAN_NAME( dorgtr, ZUNGTR )
#define lapackf77_dorm2r   FORTRAN_NAME( dorm2r, ZUNM2R )
#define lapackf77_dormbr   FORTRAN_NAME( dormbr, ZUNMBR )
#define lapackf77_dormlq   FORTRAN_NAME( dormlq, DORMLQ )
#define lapackf77_dormql   FORTRAN_NAME( dormql, DORMQL )
#define lapackf77_dormqr   FORTRAN_NAME( dormqr, DORMQR )
#define lapackf77_dormtr   FORTRAN_NAME( dormtr, DORMTR )

/* testing functions */
#define lapackf77_dbdt01   FORTRAN_NAME( dbdt01, ZBDT01 )
#define lapackf77_dget22   FORTRAN_NAME( dget22, ZGET22 )
#define lapackf77_dqpt01   FORTRAN_NAME( dqpt01, ZQPT01 )
#define lapackf77_dsyt21   FORTRAN_NAME( dsyt21, ZHET21 )
#define lapackf77_dhst01   FORTRAN_NAME( dhst01, ZHST01 )
#define lapackf77_dqrt02   FORTRAN_NAME( dqrt02, ZQRT02 )
#define lapackf77_dort01   FORTRAN_NAME( dort01, ZUNT01 )
#define lapackf77_dlarfy   FORTRAN_NAME( dlarfy, ZLARFY )
#define lapackf77_dlarfx   FORTRAN_NAME( dlarfx, ZLARFX )
#define lapackf77_dstt21   FORTRAN_NAME( dstt21, ZSTT21 )


#if defined(PRECISION_z) || defined(PRECISION_c)
#define DWORKFORZ        double *rwork,
#define DWORKFORZ_AND_LD double *rwork, magma_int_t *ldrwork,
#define WSPLIT           double *w
#else
#define DWORKFORZ 
#define DWORKFORZ_AND_LD
#define WSPLIT           double *wr, double *wi
#endif

  /*
   * BLAS functions (Alphabetical order)
   */
void     blasf77_daxpy(const int *, double *, double *, 
                       const int *, double *, const int *);
void     blasf77_dcopy(const int *, double *, const int *,
                       double *, const int *);
#if defined(PRECISION_z) || defined(PRECISION_c)
void     blasf77_ddot(double *, int *, double *, int *, 
                       double *, int *);
#endif
void     blasf77_dgemm(const char *, const char *, const int *, const int *, const int *,
                       double *, double *, const int *, 
                       double *, const int *, double *,
                       double *, const int *);
void     blasf77_dgemv(const char *, const int  *, const int *, double *, 
                       double *, const int *, double *, const int *, 
                       double *, double *, const int *);
void     blasf77_dger(int *, int *, double *, double *, int *, 
                       double *, int *, double *, int *);
void     blasf77_dsymm(const char *, const char *, const int *, const int *, 
                       double *, double *, const int *, 
                       double *, const int *, double *,
                       double *, const int *);
void     blasf77_dsymv(const char *, const int  *, double *, double *,
                       const int *, double *, const int *, double *,
                       double *, const int *);
void    blasf77_dsyr2k(const char *, const char *, const int *, const int *, 
                       double *, double *, const int *, 
                       double *, const int *, double *, 
                       double *, const int *);
void     blasf77_dher2(const char *, int *, double *, 
                       double *, int *, double *, int *, 
                       double *, int *);
void    blasf77_dsyrk( const char *, const char *, const int *, const int *, double *, 
                       double *, const int *, double *, double *, 
                       const int *);
void    blasf77_dscal( const int *, double *, double *, const int *);
#if defined(PRECISION_z) || defined(PRECISION_c)
void    blasf77_ddscal( const int *, double *, double *, const int *);
#endif
void    blasf77_dsymm( const char *, const char *, const int *, const int *, 
                       double *, double *, const int *, 
                       double *, const int *, double *,
                       double *, const int *);
void    blasf77_dsyr2k(const char *, const char *, const int *, const int *, 
                       double *, double *, const int *, 
                       double *, const int *, double *, 
                       double *, const int *);
void    blasf77_dsyrk( const char *, const char *, const int *, const int *, 
                       double *, double *, const int *, 
                       double *, double *, const int *);
void    blasf77_dswap( int *, double *, int *, double *, int *);
void    blasf77_dtrmm( const char *, const char *, const char *, const char *, 
                       const int *, const int *, double *,
                       double *, const int *, double *,const int *);
void    blasf77_dtrmv( const char *, const char *, const char *, const int *, 
                       double*,  const int *, double *, const int*);
void    blasf77_dtrsm( const char *, const char *, const char *, const char *, 
                       const int *, const int *, double *, 
                       double *, const int *, double *,const int*);
void    blasf77_dtrsv( const char *, const char *, const char *, const int *, 
                       double *, const int *, double *, const int*);

  /*
   * Lapack functions (Alphabetical order)
   */
void    lapackf77_dbdsqr(const char *uplo, magma_int_t *n, magma_int_t *nvct, 
                         magma_int_t *nru,  magma_int_t *ncc, double *D, double *E, 
                         double *VT, magma_int_t *ldvt, 
                         double *U, magma_int_t *ldu, 
                         double *C, magma_int_t *ldc, 
                         double *work, magma_int_t *info);
void    lapackf77_dgebak(const char *job, const char *side, magma_int_t *n, 
                         magma_int_t *ilo, magma_int_t *ihi, 
                         double *scale, magma_int_t *m,
                         double *v, magma_int_t *ldv, magma_int_t *info);
void    lapackf77_dgebal(const char *job, magma_int_t *n, double *A, magma_int_t *lda, 
                         magma_int_t *ilo, magma_int_t *ihi, double *scale, magma_int_t *info);
void    lapackf77_dgebd2(magma_int_t *m, magma_int_t *n, 
                         double *a, magma_int_t *lda, double *d, double *e,
                         double *tauq, double *taup,
                         double *work, magma_int_t *info);
void    lapackf77_dgebrd(magma_int_t *m, magma_int_t *n, 
                         double *a, magma_int_t *lda, double *d, double *e,
                         double *tauq, double *taup, 
                         double *work, magma_int_t *lwork, magma_int_t *info);
void     lapackf77_dgeev(const char *jobl, const char *jobr, magma_int_t *n, 
                         double *a, magma_int_t *lda, WSPLIT, 
                         double *vl, magma_int_t *ldvl, 
                         double *vr, magma_int_t *ldvr, 
                         double *work, magma_int_t *lwork, 
                         DWORKFORZ magma_int_t *info);
void    lapackf77_dgehd2(magma_int_t *n, magma_int_t *ilo, magma_int_t *ihi, 
                         double *a, magma_int_t *lda, double *tau, 
                         double *work, magma_int_t *info);
void    lapackf77_dgehrd(magma_int_t *n, magma_int_t *ilo, magma_int_t *ihi, 
                         double *a, magma_int_t *lda, double *tau,
                         double *work, magma_int_t *lwork, magma_int_t *info);
void    lapackf77_dgelqf(magma_int_t *m, magma_int_t *n, 
                         double *a, magma_int_t *lda, double *tau,
                         double *work, magma_int_t *lwork, magma_int_t *info);
void     lapackf77_dgels(const char *trans, 
                         magma_int_t *m, magma_int_t *n, magma_int_t *nrhs, 
                         double *a, magma_int_t *lda, 
                         double *b, magma_int_t *ldb,
                         double *work, magma_int_t *lwork, magma_int_t *info);
void    lapackf77_dgeqlf(magma_int_t *m, magma_int_t *n,
                         double *a, magma_int_t *lda, double *tau, 
                         double *work, magma_int_t *lwork, magma_int_t *info);
void    lapackf77_dgeqp3(magma_int_t *m, magma_int_t *n, double *a, magma_int_t *lda,
                         magma_int_t *jpvt, double *tau,
                         double *work, magma_int_t *lwork, 
                         DWORKFORZ magma_int_t *info);
void    lapackf77_dgeqrf(magma_int_t *m, magma_int_t *n,
                         double *a, magma_int_t *lda, double *tau,
                         double *work, magma_int_t *lwork, magma_int_t *info);
void    lapackf77_dgetrf(magma_int_t *m, magma_int_t *n, 
                         double *a, magma_int_t *lda, 
                         magma_int_t *ipiv, magma_int_t *info);
void    lapackf77_dgetri(magma_int_t *n,
                         double *a, magma_int_t *lda, magma_int_t *ipiv,
                         double *work, magma_int_t *lwork, magma_int_t *info);
void    lapackf77_dgetrs(const char* trans,
                         magma_int_t *n, magma_int_t *nrhs,
                         double *a, magma_int_t *lda, magma_int_t *ipiv,
                         double *b, magma_int_t *ldb, magma_int_t *info);
void    lapackf77_dgesvd(const char *jobu, const char *jobvt, 
                         magma_int_t *m, magma_int_t *n, 
                         double *a, magma_int_t *lda, 
                         double *s, double *u, magma_int_t *ldu, 
                         double *vt, magma_int_t *ldvt, 
                         double *work, magma_int_t *lwork, 
                         DWORKFORZ magma_int_t *info );
void    lapackf77_dsyev(const char *jobz, const char *uplo, magma_int_t *n, 
                         double *a, magma_int_t *lda, double *w, 
                         double *work, magma_int_t *lwork,
                         DWORKFORZ magma_int_t *info);
void    lapackf77_dsyevd(const char *jobz, const char *uplo, magma_int_t *n, 
                         double *a, magma_int_t *lda, double *w, 
                         double *work, magma_int_t *lwork,
                         DWORKFORZ_AND_LD magma_int_t *iwork, 
                         magma_int_t *liwork, magma_int_t *info);
void    lapackf77_dhegs2(int *itype, const char *uplo, int *n, 
                         double *a, int *lda, 
                         double *b, int *ldb, int *info);
void    lapackf77_dhegvd(magma_int_t *itype, const char *jobz, const char *uplo, 
                         magma_int_t *n, double *a, magma_int_t *lda,
                         double *b, magma_int_t *ldb, double *w,
                         double *work, magma_int_t *lwork, 
                         DWORKFORZ_AND_LD magma_int_t *iwork, magma_int_t *liwork,
                         magma_int_t *info);
void    lapackf77_dsytd2(const char *uplo, magma_int_t *n, 
                         double *a, magma_int_t *lda, 
                         double *d, double *e, double *tau, magma_int_t *info);
void    lapackf77_dsytrd(const char *uplo, magma_int_t *n, 
                         double *a, magma_int_t *lda, 
                         double *d, double *e, double *tau, 
                         double *work, magma_int_t *lwork, magma_int_t *info);
void    lapackf77_dsbtrd(const char *vect, const char *uplo, magma_int_t *n, magma_int_t *kd, 
                         double *ab, magma_int_t *ldab, double *d__, double *e, 
                         double *q, magma_int_t *ldq, double *work, 
                         magma_int_t *info);
void    lapackf77_dhseqr(const char *job, const char *compz, magma_int_t *n, 
                         magma_int_t *ilo, magma_int_t *ihi, 
                         double *H, magma_int_t *ldh, WSPLIT, 
                         double *Z, magma_int_t *ldz, 
                         double *work, magma_int_t *lwork, magma_int_t *info);
void    lapackf77_dlacpy(const char *uplo, magma_int_t *m, magma_int_t *n, 
                         const double *a, magma_int_t *lda, 
                         double *b, magma_int_t *ldb);
void    lapackf77_dlacgv(magma_int_t *n, double *x, magma_int_t *incx);
double  lapackf77_dlange(const char *norm, magma_int_t *m, magma_int_t *n, 
                         const double *a, magma_int_t *lda, double *work);
double  lapackf77_dlansy(const char *norm, const char *uplo, magma_int_t *n, 
                         const double *a, magma_int_t *lda, double * work);
double  lapackf77_dlansy(const char *norm, const char *uplo, magma_int_t *n, 
                         const double *a, magma_int_t *lda, double * work);
void    lapackf77_dlarfb(const char *side, const char *trans, const char *direct, 
                         const char *storev, magma_int_t *m, magma_int_t *n, magma_int_t *k, 
                         const double *v, magma_int_t *ldv, 
                         const double *t, magma_int_t *ldt, 
                         double *c, magma_int_t *ldc, 
                         double *work, magma_int_t *ldwork);
void    lapackf77_dlarfg(magma_int_t *n, double *alpha, 
                         double *x, magma_int_t *incx, double *tau);
void    lapackf77_dlarft(const char *direct, const char *storev, magma_int_t *n, magma_int_t *k, 
                         double *v, magma_int_t *ldv, const double *tau, 
                         double *t, magma_int_t *ldt);
void    lapackf77_dlarnv(magma_int_t *idist, magma_int_t *iseed, magma_int_t *n, 
                         double *x);
void    lapackf77_dlartg(double *F, double *G, double *cs, 
                         double *SN, double *R);
void    lapackf77_dlascl(const char *type, magma_int_t *kl, magma_int_t *ku, 
                         double *cfrom, double *cto, 
                         magma_int_t *m, magma_int_t *n, 
                         double *A, magma_int_t *lda, magma_int_t *info);
void    lapackf77_dlaset(const char *uplo, magma_int_t *m, magma_int_t *n, 
                         double *alpha, double *beta,
                         double *A, magma_int_t *lda);
void    lapackf77_dlaswp(magma_int_t *n, double *a, magma_int_t *lda, 
                         magma_int_t *k1, magma_int_t *k2, magma_int_t *ipiv,
                         magma_int_t *incx);
void    lapackf77_dlatrd(const char *uplo, magma_int_t *n, magma_int_t *nb, 
                         double *a, magma_int_t *lda, double *e,
                         double *tau, double *work, magma_int_t *ldwork);
void    lapackf77_dlabrd(magma_int_t *m, magma_int_t *n, magma_int_t *nb, 
                         double *a, magma_int_t *lda, double *d__, double *e, 
                         double *tauq, double *taup,
                         double *x, magma_int_t *ldx,
                         double *y, magma_int_t *ldy);
void    lapackf77_dpotrf(const char *uplo, magma_int_t *n, 
                         double *a, magma_int_t *lda, magma_int_t *info);
void    lapackf77_dpotrs(const char *uplo, magma_int_t *n, magma_int_t *nrhs,
                         double *a, magma_int_t *lda,
                         double *b, magma_int_t *ldb, magma_int_t *info);
void    lapackf77_dpotri(const char *uplo, magma_int_t *n, 
                         double *a, magma_int_t *lda, magma_int_t *info);
void    lapackf77_dlauum(const char *uplo, magma_int_t *n, 
                         double *a, magma_int_t *lda, magma_int_t *info);
void    lapackf77_dtrevc(const char *side, const char *howmny, magma_int_t *select, magma_int_t *n, 
                         double *T,  magma_int_t *ldt,  double *VL, magma_int_t *ldvl,
                         double *VR, magma_int_t *ldvr, magma_int_t *MM, magma_int_t *M, 
                         double *work, DWORKFORZ magma_int_t *info);
void    lapackf77_dstebz(const char *range, const char *order, magma_int_t *n, double *vl, double *vu,
                         magma_int_t *il, magma_int_t *iu, double *abstol,
                         double *d__, double *e, magma_int_t *m, magma_int_t *nsplit,
                         double *w, magma_int_t *iblock, magma_int_t *isplit, double *work,
                         magma_int_t *iwork, magma_int_t *info);
double  lapackf77_dlamc3(double* a, double* b);
void    lapackf77_dlamrg(magma_int_t* n1, magma_int_t* n2, double* a, 
                         magma_int_t* dtrd1, magma_int_t* dtrd2, magma_int_t* index);
void    lapackf77_dlaed4(magma_int_t* n, magma_int_t* i, double* d, double* z,
                         double* delta, double* rho, double* dlam, magma_int_t* info);
void    lapackf77_dsteqr(const char *compz, magma_int_t *n, double *D, double *E, 
                         double *Z, magma_int_t *ldz, 
                         double *work, magma_int_t *info);
void    lapackf77_dstedc(const char *compz, magma_int_t *n, double *D, double *E, 
                         double *Z, magma_int_t *ldz, 
                         double *work, magma_int_t *ldwork, 
                         DWORKFORZ_AND_LD magma_int_t *iwork, magma_int_t *liwork,
                         magma_int_t *info);
void    lapackf77_dstein(magma_int_t *n, double *d__, double *e, 
                         magma_int_t *m, double *w, magma_int_t *iblock, magma_int_t *isplit, 
                         double *z__, magma_int_t *ldz, double *work, magma_int_t *iwork, 
                         magma_int_t *ifail, magma_int_t *info);
void    lapackf77_dstemr(const char *jobz, const char *range, magma_int_t *n, double *d__, double *e, 
                         double *vl, double *vu, magma_int_t *il, magma_int_t *iu, magma_int_t *m,
                         double *w, double *z__, magma_int_t *ldz, magma_int_t *nzc, 
                         magma_int_t *isuppz, magma_int_t *tryrac, double *work, magma_int_t *lwork, 
                         magma_int_t *iwork, magma_int_t *liwork, magma_int_t *info);
void    lapackf77_dtrtri(const char *uplo, const char *diag, magma_int_t *n,
                         double *a, magma_int_t *lda, magma_int_t *info);
#if defined(PRECISION_z) || defined(PRECISION_c)
void    lapackf77_dsymv(const char *uplo, const magma_int_t *N, const double *alpha, 
                        const double *A, const magma_int_t *lda, 
                        const double *X, const magma_int_t *incX,
                        const double *beta, 
                        double *Y, const magma_int_t *incY);
#endif
void    lapackf77_dorg2r(magma_int_t *m, magma_int_t *n, magma_int_t *k, 
                         double *a, magma_int_t *lda,
                         const double *tau, double *work,
                         magma_int_t *info);
void    lapackf77_dorgbr(const char *vect, magma_int_t *m, magma_int_t *n, magma_int_t *k,
                         double *a, magma_int_t *lda, const double *tau,
                         double *work, magma_int_t *lwork, magma_int_t *info);
void    lapackf77_dorghr(magma_int_t *n, magma_int_t *ilo, magma_int_t *ihi, 
                         double *a, magma_int_t *lda, const double *tau,
                         double *work, magma_int_t *lwork, magma_int_t *info);
void    lapackf77_dorglq(magma_int_t *m, magma_int_t *n, magma_int_t *k, 
                         double *a, magma_int_t *lda, const double *tau, 
                         double *work, magma_int_t *ldwork, magma_int_t *info);
void    lapackf77_dungql(magma_int_t *, magma_int_t *, magma_int_t *,
                         double *, magma_int_t *, double *, 
                         double *, magma_int_t *, magma_int_t *);
void    lapackf77_dorgqr(magma_int_t *m, magma_int_t *n, magma_int_t *k, 
                         double *a, magma_int_t *lda, const double *tau, 
                         double *work, magma_int_t *ldwork, magma_int_t *info);
void    lapackf77_dorgtr(const char *uplo, magma_int_t *n, 
                         double *a, magma_int_t *lda, const double *tau, 
                         double *work, magma_int_t *lwork, magma_int_t *info);
void    lapackf77_dorm2r(const char *side, const char *trans, 
                         magma_int_t *m, magma_int_t *n, magma_int_t *k, 
                         const double *a, magma_int_t *lda, 
                         const double *tau, double *c, magma_int_t *ldc,
                         double *work, magma_int_t *info);
void    lapackf77_dormbr(const char *vect, const char *side, const char *trans,
                         magma_int_t *M, magma_int_t *N, magma_int_t *K, 
                         double *A, magma_int_t *lda, double *Tau,
                         double *C, magma_int_t *ldc, 
                         double *work, magma_int_t *ldwork, magma_int_t *info);
void    lapackf77_dormlq(const char *side, const char *trans, 
                         magma_int_t *m, magma_int_t *n, magma_int_t *k,
                         const double *a, magma_int_t *lda, 
                         const double *tau, double *c, magma_int_t *ldc, 
                         double *work, magma_int_t *lwork, magma_int_t *info);
void    lapackf77_dormql(const char *side, const char *trans, 
                         magma_int_t *m, magma_int_t *n, magma_int_t *k,
                         const double *a, magma_int_t *lda, 
                         const double *tau, double *c, magma_int_t *ldc,
                         double *work, magma_int_t *lwork, magma_int_t *info);
void    lapackf77_dormqr(const char *side, const char *trans, 
                         magma_int_t *m, magma_int_t *n, magma_int_t *k, 
                         const double *a, magma_int_t *lda, 
                         const double *tau, double *c, magma_int_t *ldc, 
                         double *work, magma_int_t *lwork, magma_int_t *info);
void    lapackf77_dormtr(const char *side, const char *uplo, const char *trans,
                         magma_int_t *M, magma_int_t *N,
                         double *A, magma_int_t *lda, double *Tau,
                         double *C, magma_int_t *ldc, 
                         double *work, magma_int_t *ldwork, magma_int_t *info);


  /*
   * Testing functions
   */

#if defined(PRECISION_z) || defined(PRECISION_c)

void    lapackf77_dbdt01(int *m, int *n, int *kd, double *A, int *lda, 
                         double *Q, int *ldq, double *D, double *E, 
                         double *PT, int *ldpt, double *work, 
                         double *rwork, double *resid);
void    lapackf77_dget22(const char *transa, const char *transe, const char *transw, int *n,
                         double *a, int *lda, double *e, int *lde,
                         double *w, double *work,
                         double *rwork, double *result);
void    lapackf77_dsyt21(int *itype, const char *uplo, int *n, int *kband, 
                         double *A, int *lda, double *D, double *E, 
                         double *U, int *ldu, double *V, int *ldv, 
                         double *TAU, double *work,
                         double *rwork, double *result);
void    lapackf77_dhst01(int *n, int *ilo, int *ihi, double *A, int *lda, 
                         double *H, int *ldh, double *Q, int *ldq,
                         double *work, int *lwork, double *rwork, double *result);
void    lapackf77_dstt21(int *n, int *kband, double *AD, double *AE, double *SD,
                         double *SE, double *U, int *ldu, 
                         double *work, double *rwork, double *result);
void    lapackf77_dort01(const char *rowcol, int *m, int *n, double *U, int *ldu,
                         double *work, int *lwork, double *rwork, double *resid);

#else

void    lapackf77_dbdt01(int *m, int *n, int *kd, double *A, int *lda, 
                         double *Q, int *ldq, double *D, double *E, 
                         double *PT, int *ldpt, 
                         double *work, double *resid);
void    lapackf77_dget22(const char *transa, const char *transe, const char *transw, int *n,
                         double *a, int *lda, double *e, int *lde,
                         double *wr, double *wi, 
                         double *work, double *result);
void    lapackf77_dsyt21(int *itype, const char *uplo, int *n, int *kband, 
                         double *A, int *lda, double *D, double *E,
                         double *U, int *ldu, double *V, int *ldv, 
                         double *TAU, double *work, double *result);
void    lapackf77_dhst01(int *n, int *ilo, int *ihi, double *A, int *lda, 
                         double *H, int *ldh, double *Q, int *ldq, 
                         double *work, int *lwork, double *result);
void    lapackf77_dstt21(int *n, int *kband, double *AD, double *AE, double *SD, 
                         double *SE, double *U, int *ldu, 
                         double *work, double *result);
void    lapackf77_dort01(const char *rowcol, int *m, int *n, double *U, int *ldu,
                         double *work, int *lwork, double *resid);
#endif

void    lapackf77_dlarfy(const char *uplo, int *N, double *V, int *incv, 
                         double *tau, double *C, int *ldc, 
                         double *work);
void    lapackf77_dlarfx(const char *, int *, int *, 
                         double *, double *, 
                         double *, int *, double *);
double  lapackf77_dqpt01(int *m, int *n, int *k, double *a,
                         double *af, int *lda, double *tau, int *jpvt,
                         double *work, int *lwork);
void    lapackf77_dqrt02(int *m, int *n, int *k, double *A, double *AF,
                         double *Q, double *R, int *lda, 
                         double *TAU, double *work, int *lwork,
                         double *rwork, double *result);

#ifdef __cplusplus
}
#endif

#undef DWORKFORZ 
#undef DWORKFORZ_AND_LD
#undef WSPLIT
#undef PRECISION_d
#endif /* MAGMA ZLAPACK */
