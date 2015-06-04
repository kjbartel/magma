/*
 *   -- MAGMA (version 1.2.0) --
 *      Univ. of Tennessee, Knoxville
 *      Univ. of California, Berkeley
 *      Univ. of Colorado, Denver
 *      May 2012
 *
 * @generated c Tue May 15 18:17:06 2012
 */

#ifndef MAGMA_CLAPACK_H
#define MAGMA_CLAPACK_H

#define PRECISION_c
#ifdef __cplusplus
extern "C" {
#endif

/* ////////////////////////////////////////////////////////////////////////////
   -- LAPACK Externs used in MAGMA
*/

#define blasf77_caxpy      FORTRAN_NAME( caxpy,  ZAXPY  )
#define blasf77_ccopy      FORTRAN_NAME( ccopy,  ZCOPY  )

/* complex versions use C wrapper to return value; no name mangling. */
#if  defined(PRECISION_z) || defined(PRECISION_c)    
#define blasf77_cdotc      cdotc
#else
#define blasf77_cdotc      FORTRAN_NAME( cdotc,  CDOTC  )
#endif

#define blasf77_cgemm      FORTRAN_NAME( cgemm,  CGEMM  )
#define blasf77_cgemv      FORTRAN_NAME( cgemv,  ZGEMV  )
#define blasf77_chemm      FORTRAN_NAME( chemm,  CHEMM  )
#define blasf77_chemv      FORTRAN_NAME( chemv,  CHEMV  )
#define blasf77_cher2k     FORTRAN_NAME( cher2k, CHER2K )
#define blasf77_cher2      FORTRAN_NAME( cher2,  CHER2  )
#define blasf77_cherk      FORTRAN_NAME( cherk,  CHERK  )
#define blasf77_cscal      FORTRAN_NAME( cscal,  ZSCAL  )
#define blasf77_csscal     FORTRAN_NAME( csscal, ZDSCAL ) 
#define blasf77_csymm      FORTRAN_NAME( csymm,  CSYMM  )
#define blasf77_csyr2k     FORTRAN_NAME( csyr2k, CSYR2K )
#define blasf77_csyrk      FORTRAN_NAME( csyrk,  CSYRK  )
#define blasf77_cswap      FORTRAN_NAME( cswap,  ZSWAP  )
#define blasf77_ctrmm      FORTRAN_NAME( ctrmm,  CTRMM  )
#define blasf77_ctrmv      FORTRAN_NAME( ctrmv,  ZTRMV  )
#define blasf77_ctrsm      FORTRAN_NAME( ctrsm,  CTRSM  )
#define blasf77_ctrsv      FORTRAN_NAME( ctrsv,  CTRSV  )
#define blasf77_cgeru      FORTRAN_NAME( cgeru,  ZGERU  )

#define lapackf77_cbdsqr   FORTRAN_NAME( cbdsqr, ZBDSQR )
#define lapackf77_cgebak   FORTRAN_NAME( cgebak, ZGEBAK )
#define lapackf77_cgebal   FORTRAN_NAME( cgebal, ZGEBAL )
#define lapackf77_cgebd2   FORTRAN_NAME( cgebd2, ZGEBD2 )
#define lapackf77_cgebrd   FORTRAN_NAME( cgebrd, CGEBRD )
#define lapackf77_cgeev    FORTRAN_NAME( cgeev,  CGEEV  )
#define lapackf77_cgehd2   FORTRAN_NAME( cgehd2, ZGEHD2 )
#define lapackf77_cgehrd   FORTRAN_NAME( cgehrd, CGEHRD )
#define lapackf77_cgelqf   FORTRAN_NAME( cgelqf, CGELQF )
#define lapackf77_cgels    FORTRAN_NAME( cgels,  CGELS  )
#define lapackf77_cgeqlf   FORTRAN_NAME( cgeqlf, ZGEQLF )
#define lapackf77_cgeqp3   FORTRAN_NAME( cgeqp3, CGEQP3 )
#define lapackf77_cgeqrf   FORTRAN_NAME( cgeqrf, CGEQRF )
#define lapackf77_cgesvd   FORTRAN_NAME( cgesvd, CGESVD )
#define lapackf77_cgetrf   FORTRAN_NAME( cgetrf, CGETRF )
#define lapackf77_cgetri   FORTRAN_NAME( cgetri, CGETRI )
#define lapackf77_cgetrs   FORTRAN_NAME( cgetrs, CGETRS )
#define lapackf77_cheev    FORTRAN_NAME( cheev,  CHEEV  )
#define lapackf77_cheevd   FORTRAN_NAME( cheevd, CHEEVD )
#define lapackf77_chegs2   FORTRAN_NAME( chegs2, ZHEGS2 )
#define lapackf77_chegvd   FORTRAN_NAME( chegvd, CHEGVD )
#define lapackf77_chetd2   FORTRAN_NAME( chetd2, ZHETD2 )
#define lapackf77_chetrd   FORTRAN_NAME( chetrd, CHETRD )
#define lapackf77_chbtrd   FORTRAN_NAME( chbtrd, CHBTRD )
#define lapackf77_chseqr   FORTRAN_NAME( chseqr, ZHSEQR )
#define lapackf77_clacpy   FORTRAN_NAME( clacpy, ZLACPY )
#define lapackf77_clacgv   FORTRAN_NAME( clacgv, ZLACGV )
#define lapackf77_clange   FORTRAN_NAME( clange, CLANGE )
#define lapackf77_clanhe   FORTRAN_NAME( clanhe, CLANHE )
#define lapackf77_clansy   FORTRAN_NAME( clansy, CLANSY )
#define lapackf77_clarfb   FORTRAN_NAME( clarfb, CLARFB )
#define lapackf77_clarfg   FORTRAN_NAME( clarfg, ZLARFG )
#define lapackf77_clarft   FORTRAN_NAME( clarft, ZLARFT )
#define lapackf77_clarnv   FORTRAN_NAME( clarnv, ZLARNV )
#define lapackf77_clartg   FORTRAN_NAME( clartg, ZLARTG )
#define lapackf77_clascl   FORTRAN_NAME( clascl, ZLASCL )
#define lapackf77_claset   FORTRAN_NAME( claset, ZLASET )
#define lapackf77_claswp   FORTRAN_NAME( claswp, ZLASWP )
#define lapackf77_clatrd   FORTRAN_NAME( clatrd, CLATRD )
#define lapackf77_clabrd   FORTRAN_NAME( clabrd, CLABRD )
#define lapackf77_clauum   FORTRAN_NAME( clauum, ZLAUUM )
#define lapackf77_cpotrf   FORTRAN_NAME( cpotrf, CPOTRF )
#define lapackf77_cpotrs   FORTRAN_NAME( cpotrs, CPOTRS )
#define lapackf77_cpotri   FORTRAN_NAME( cpotri, ZPOTRI )
#define lapackf77_ctrevc   FORTRAN_NAME( ctrevc, ZTREVC )
#define lapackf77_sstebz   FORTRAN_NAME( sstebz, DSTEBZ )
#define lapackf77_slamc3   FORTRAN_NAME( slamc3, DLAMC3 )
#define lapackf77_slaed4   FORTRAN_NAME( slaed4, DLAED4 )
#define lapackf77_slamrg   FORTRAN_NAME( slamrg, DLAMRG )
#define lapackf77_ctrtri   FORTRAN_NAME( ctrtri, ZTRTRI )
#define lapackf77_csteqr   FORTRAN_NAME( csteqr, ZSTEQR )
#define lapackf77_cstedc   FORTRAN_NAME( cstedc, ZSTEDC )
#define lapackf77_cstein   FORTRAN_NAME( cstein, ZSTEIN )
#define lapackf77_cstemr   FORTRAN_NAME( cstemr, ZSTEMR )
#define lapackf77_csymv    FORTRAN_NAME( csymv,  ZSYMV  )
#define lapackf77_cung2r   FORTRAN_NAME( cung2r, ZUNG2R )
#define lapackf77_cungbr   FORTRAN_NAME( cungbr, ZUNGBR )
#define lapackf77_cunghr   FORTRAN_NAME( cunghr, CUNGHR )
#define lapackf77_cunglq   FORTRAN_NAME( cunglq, CUNGLQ )
#define lapackf77_cungql   FORTRAN_NAME( cungql, ZUNGQL )
#define lapackf77_cungqr   FORTRAN_NAME( cungqr, CUNGQR )
#define lapackf77_cungtr   FORTRAN_NAME( cungtr, ZUNGTR )
#define lapackf77_cunm2r   FORTRAN_NAME( cunm2r, ZUNM2R )
#define lapackf77_cunmbr   FORTRAN_NAME( cunmbr, ZUNMBR )
#define lapackf77_cunmlq   FORTRAN_NAME( cunmlq, CUNMLQ )
#define lapackf77_cunmql   FORTRAN_NAME( cunmql, CUNMQL )
#define lapackf77_cunmqr   FORTRAN_NAME( cunmqr, CUNMQR )
#define lapackf77_cunmtr   FORTRAN_NAME( cunmtr, CUNMTR )

/* testing functions */
#define lapackf77_cbdt01   FORTRAN_NAME( cbdt01, ZBDT01 )
#define lapackf77_cget22   FORTRAN_NAME( cget22, ZGET22 )
#define lapackf77_cqpt01   FORTRAN_NAME( cqpt01, ZQPT01 )
#define lapackf77_chet21   FORTRAN_NAME( chet21, ZHET21 )
#define lapackf77_chst01   FORTRAN_NAME( chst01, ZHST01 )
#define lapackf77_cqrt02   FORTRAN_NAME( cqrt02, ZQRT02 )
#define lapackf77_cunt01   FORTRAN_NAME( cunt01, ZUNT01 )
#define lapackf77_clarfy   FORTRAN_NAME( clarfy, ZLARFY )
#define lapackf77_clarfx   FORTRAN_NAME( clarfx, ZLARFX )
#define lapackf77_cstt21   FORTRAN_NAME( cstt21, ZSTT21 )


#if defined(PRECISION_z) || defined(PRECISION_c)
#define DWORKFORZ        float *rwork,
#define DWORKFORZ_AND_LD float *rwork, magma_int_t *ldrwork,
#define WSPLIT           cuFloatComplex *w
#else
#define DWORKFORZ 
#define DWORKFORZ_AND_LD
#define WSPLIT           float *wr, float *wi
#endif

  /*
   * BLAS functions (Alphabetical order)
   */
void     blasf77_caxpy(const int *, cuFloatComplex *, cuFloatComplex *, 
                       const int *, cuFloatComplex *, const int *);
void     blasf77_ccopy(const int *, cuFloatComplex *, const int *,
                       cuFloatComplex *, const int *);
#if defined(PRECISION_z) || defined(PRECISION_c)
void     blasf77_cdotc(cuFloatComplex *, int *, cuFloatComplex *, int *, 
                       cuFloatComplex *, int *);
#endif
void     blasf77_cgemm(const char *, const char *, const int *, const int *, const int *,
                       cuFloatComplex *, cuFloatComplex *, const int *, 
                       cuFloatComplex *, const int *, cuFloatComplex *,
                       cuFloatComplex *, const int *);
void     blasf77_cgemv(const char *, const int  *, const int *, cuFloatComplex *, 
                       cuFloatComplex *, const int *, cuFloatComplex *, const int *, 
                       cuFloatComplex *, cuFloatComplex *, const int *);
void     blasf77_cgeru(int *, int *, cuFloatComplex *, cuFloatComplex *, int *, 
                       cuFloatComplex *, int *, cuFloatComplex *, int *);
void     blasf77_chemm(const char *, const char *, const int *, const int *, 
                       cuFloatComplex *, cuFloatComplex *, const int *, 
                       cuFloatComplex *, const int *, cuFloatComplex *,
                       cuFloatComplex *, const int *);
void     blasf77_chemv(const char *, const int  *, cuFloatComplex *, cuFloatComplex *,
                       const int *, cuFloatComplex *, const int *, cuFloatComplex *,
                       cuFloatComplex *, const int *);
void    blasf77_cher2k(const char *, const char *, const int *, const int *, 
                       cuFloatComplex *, cuFloatComplex *, const int *, 
                       cuFloatComplex *, const int *, float *, 
                       cuFloatComplex *, const int *);
void     blasf77_cher2(const char *, int *, cuFloatComplex *, 
                       cuFloatComplex *, int *, cuFloatComplex *, int *, 
                       cuFloatComplex *, int *);
void    blasf77_cherk( const char *, const char *, const int *, const int *, float *, 
                       cuFloatComplex *, const int *, float *, cuFloatComplex *, 
                       const int *);
void    blasf77_cscal( const int *, cuFloatComplex *, cuFloatComplex *, const int *);
#if defined(PRECISION_z) || defined(PRECISION_c)
void    blasf77_csscal( const int *, float *, cuFloatComplex *, const int *);
#endif
void    blasf77_csymm( const char *, const char *, const int *, const int *, 
                       cuFloatComplex *, cuFloatComplex *, const int *, 
                       cuFloatComplex *, const int *, cuFloatComplex *,
                       cuFloatComplex *, const int *);
void    blasf77_csyr2k(const char *, const char *, const int *, const int *, 
                       cuFloatComplex *, cuFloatComplex *, const int *, 
                       cuFloatComplex *, const int *, cuFloatComplex *, 
                       cuFloatComplex *, const int *);
void    blasf77_csyrk( const char *, const char *, const int *, const int *, 
                       cuFloatComplex *, cuFloatComplex *, const int *, 
                       cuFloatComplex *, cuFloatComplex *, const int *);
void    blasf77_cswap( int *, cuFloatComplex *, int *, cuFloatComplex *, int *);
void    blasf77_ctrmm( const char *, const char *, const char *, const char *, 
                       const int *, const int *, cuFloatComplex *,
                       cuFloatComplex *, const int *, cuFloatComplex *,const int *);
void    blasf77_ctrmv( const char *, const char *, const char *, const int *, 
                       cuFloatComplex*,  const int *, cuFloatComplex *, const int*);
void    blasf77_ctrsm( const char *, const char *, const char *, const char *, 
                       const int *, const int *, cuFloatComplex *, 
                       cuFloatComplex *, const int *, cuFloatComplex *,const int*);
void    blasf77_ctrsv( const char *, const char *, const char *, const int *, 
                       cuFloatComplex *, const int *, cuFloatComplex *, const int*);

  /*
   * Lapack functions (Alphabetical order)
   */
void    lapackf77_cbdsqr(const char *uplo, magma_int_t *n, magma_int_t *nvct, 
                         magma_int_t *nru,  magma_int_t *ncc, float *D, float *E, 
                         cuFloatComplex *VT, magma_int_t *ldvt, 
                         cuFloatComplex *U, magma_int_t *ldu, 
                         cuFloatComplex *C, magma_int_t *ldc, 
                         float *work, magma_int_t *info);
void    lapackf77_cgebak(const char *job, const char *side, magma_int_t *n, 
                         magma_int_t *ilo, magma_int_t *ihi, 
                         float *scale, magma_int_t *m,
                         cuFloatComplex *v, magma_int_t *ldv, magma_int_t *info);
void    lapackf77_cgebal(const char *job, magma_int_t *n, cuFloatComplex *A, magma_int_t *lda, 
                         magma_int_t *ilo, magma_int_t *ihi, float *scale, magma_int_t *info);
void    lapackf77_cgebd2(magma_int_t *m, magma_int_t *n, 
                         cuFloatComplex *a, magma_int_t *lda, float *d, float *e,
                         cuFloatComplex *tauq, cuFloatComplex *taup,
                         cuFloatComplex *work, magma_int_t *info);
void    lapackf77_cgebrd(magma_int_t *m, magma_int_t *n, 
                         cuFloatComplex *a, magma_int_t *lda, float *d, float *e,
                         cuFloatComplex *tauq, cuFloatComplex *taup, 
                         cuFloatComplex *work, magma_int_t *lwork, magma_int_t *info);
void     lapackf77_cgeev(const char *jobl, const char *jobr, magma_int_t *n, 
                         cuFloatComplex *a, magma_int_t *lda, WSPLIT, 
                         cuFloatComplex *vl, magma_int_t *ldvl, 
                         cuFloatComplex *vr, magma_int_t *ldvr, 
                         cuFloatComplex *work, magma_int_t *lwork, 
                         DWORKFORZ magma_int_t *info);
void    lapackf77_cgehd2(magma_int_t *n, magma_int_t *ilo, magma_int_t *ihi, 
                         cuFloatComplex *a, magma_int_t *lda, cuFloatComplex *tau, 
                         cuFloatComplex *work, magma_int_t *info);
void    lapackf77_cgehrd(magma_int_t *n, magma_int_t *ilo, magma_int_t *ihi, 
                         cuFloatComplex *a, magma_int_t *lda, cuFloatComplex *tau,
                         cuFloatComplex *work, magma_int_t *lwork, magma_int_t *info);
void    lapackf77_cgelqf(magma_int_t *m, magma_int_t *n, 
                         cuFloatComplex *a, magma_int_t *lda, cuFloatComplex *tau,
                         cuFloatComplex *work, magma_int_t *lwork, magma_int_t *info);
void     lapackf77_cgels(const char *trans, 
                         magma_int_t *m, magma_int_t *n, magma_int_t *nrhs, 
                         cuFloatComplex *a, magma_int_t *lda, 
                         cuFloatComplex *b, magma_int_t *ldb,
                         cuFloatComplex *work, magma_int_t *lwork, magma_int_t *info);
void    lapackf77_cgeqlf(magma_int_t *m, magma_int_t *n,
                         cuFloatComplex *a, magma_int_t *lda, cuFloatComplex *tau, 
                         cuFloatComplex *work, magma_int_t *lwork, magma_int_t *info);
void    lapackf77_cgeqp3(magma_int_t *m, magma_int_t *n, cuFloatComplex *a, magma_int_t *lda,
                         magma_int_t *jpvt, cuFloatComplex *tau,
                         cuFloatComplex *work, magma_int_t *lwork, 
                         DWORKFORZ magma_int_t *info);
void    lapackf77_cgeqrf(magma_int_t *m, magma_int_t *n,
                         cuFloatComplex *a, magma_int_t *lda, cuFloatComplex *tau,
                         cuFloatComplex *work, magma_int_t *lwork, magma_int_t *info);
void    lapackf77_cgetrf(magma_int_t *m, magma_int_t *n, 
                         cuFloatComplex *a, magma_int_t *lda, 
                         magma_int_t *ipiv, magma_int_t *info);
void    lapackf77_cgetri(magma_int_t *n,
                         cuFloatComplex *a, magma_int_t *lda, magma_int_t *ipiv,
                         cuFloatComplex *work, magma_int_t *lwork, magma_int_t *info);
void    lapackf77_cgetrs(const char* trans,
                         magma_int_t *n, magma_int_t *nrhs,
                         cuFloatComplex *a, magma_int_t *lda, magma_int_t *ipiv,
                         cuFloatComplex *b, magma_int_t *ldb, magma_int_t *info);
void    lapackf77_cgesvd(const char *jobu, const char *jobvt, 
                         magma_int_t *m, magma_int_t *n, 
                         cuFloatComplex *a, magma_int_t *lda, 
                         float *s, cuFloatComplex *u, magma_int_t *ldu, 
                         cuFloatComplex *vt, magma_int_t *ldvt, 
                         cuFloatComplex *work, magma_int_t *lwork, 
                         DWORKFORZ magma_int_t *info );
void    lapackf77_cheev(const char *jobz, const char *uplo, magma_int_t *n, 
                         cuFloatComplex *a, magma_int_t *lda, float *w, 
                         cuFloatComplex *work, magma_int_t *lwork,
                         DWORKFORZ magma_int_t *info);
void    lapackf77_cheevd(const char *jobz, const char *uplo, magma_int_t *n, 
                         cuFloatComplex *a, magma_int_t *lda, float *w, 
                         cuFloatComplex *work, magma_int_t *lwork,
                         DWORKFORZ_AND_LD magma_int_t *iwork, 
                         magma_int_t *liwork, magma_int_t *info);
void    lapackf77_chegs2(int *itype, const char *uplo, int *n, 
                         cuFloatComplex *a, int *lda, 
                         cuFloatComplex *b, int *ldb, int *info);
void    lapackf77_chegvd(magma_int_t *itype, const char *jobz, const char *uplo, 
                         magma_int_t *n, cuFloatComplex *a, magma_int_t *lda,
                         cuFloatComplex *b, magma_int_t *ldb, float *w,
                         cuFloatComplex *work, magma_int_t *lwork, 
                         DWORKFORZ_AND_LD magma_int_t *iwork, magma_int_t *liwork,
                         magma_int_t *info);
void    lapackf77_chetd2(const char *uplo, magma_int_t *n, 
                         cuFloatComplex *a, magma_int_t *lda, 
                         float *d, float *e, cuFloatComplex *tau, magma_int_t *info);
void    lapackf77_chetrd(const char *uplo, magma_int_t *n, 
                         cuFloatComplex *a, magma_int_t *lda, 
                         float *d, float *e, cuFloatComplex *tau, 
                         cuFloatComplex *work, magma_int_t *lwork, magma_int_t *info);
void    lapackf77_chbtrd(const char *vect, const char *uplo, magma_int_t *n, magma_int_t *kd, 
                         cuFloatComplex *ab, magma_int_t *ldab, float *d__, float *e, 
                         cuFloatComplex *q, magma_int_t *ldq, cuFloatComplex *work, 
                         magma_int_t *info);
void    lapackf77_chseqr(const char *job, const char *compz, magma_int_t *n, 
                         magma_int_t *ilo, magma_int_t *ihi, 
                         cuFloatComplex *H, magma_int_t *ldh, WSPLIT, 
                         cuFloatComplex *Z, magma_int_t *ldz, 
                         cuFloatComplex *work, magma_int_t *lwork, magma_int_t *info);
void    lapackf77_clacpy(const char *uplo, magma_int_t *m, magma_int_t *n, 
                         const cuFloatComplex *a, magma_int_t *lda, 
                         cuFloatComplex *b, magma_int_t *ldb);
void    lapackf77_clacgv(magma_int_t *n, cuFloatComplex *x, magma_int_t *incx);
float  lapackf77_clange(const char *norm, magma_int_t *m, magma_int_t *n, 
                         const cuFloatComplex *a, magma_int_t *lda, float *work);
float  lapackf77_clanhe(const char *norm, const char *uplo, magma_int_t *n, 
                         const cuFloatComplex *a, magma_int_t *lda, float * work);
float  lapackf77_clansy(const char *norm, const char *uplo, magma_int_t *n, 
                         const cuFloatComplex *a, magma_int_t *lda, float * work);
void    lapackf77_clarfb(const char *side, const char *trans, const char *direct, 
                         const char *storev, magma_int_t *m, magma_int_t *n, magma_int_t *k, 
                         const cuFloatComplex *v, magma_int_t *ldv, 
                         const cuFloatComplex *t, magma_int_t *ldt, 
                         cuFloatComplex *c, magma_int_t *ldc, 
                         cuFloatComplex *work, magma_int_t *ldwork);
void    lapackf77_clarfg(magma_int_t *n, cuFloatComplex *alpha, 
                         cuFloatComplex *x, magma_int_t *incx, cuFloatComplex *tau);
void    lapackf77_clarft(const char *direct, const char *storev, magma_int_t *n, magma_int_t *k, 
                         cuFloatComplex *v, magma_int_t *ldv, const cuFloatComplex *tau, 
                         cuFloatComplex *t, magma_int_t *ldt);
void    lapackf77_clarnv(magma_int_t *idist, magma_int_t *iseed, magma_int_t *n, 
                         cuFloatComplex *x);
void    lapackf77_clartg(cuFloatComplex *F, cuFloatComplex *G, float *cs, 
                         cuFloatComplex *SN, cuFloatComplex *R);
void    lapackf77_clascl(const char *type, magma_int_t *kl, magma_int_t *ku, 
                         float *cfrom, float *cto, 
                         magma_int_t *m, magma_int_t *n, 
                         cuFloatComplex *A, magma_int_t *lda, magma_int_t *info);
void    lapackf77_claset(const char *uplo, magma_int_t *m, magma_int_t *n, 
                         cuFloatComplex *alpha, cuFloatComplex *beta,
                         cuFloatComplex *A, magma_int_t *lda);
void    lapackf77_claswp(magma_int_t *n, cuFloatComplex *a, magma_int_t *lda, 
                         magma_int_t *k1, magma_int_t *k2, magma_int_t *ipiv,
                         magma_int_t *incx);
void    lapackf77_clatrd(const char *uplo, magma_int_t *n, magma_int_t *nb, 
                         cuFloatComplex *a, magma_int_t *lda, float *e,
                         cuFloatComplex *tau, cuFloatComplex *work, magma_int_t *ldwork);
void    lapackf77_clabrd(magma_int_t *m, magma_int_t *n, magma_int_t *nb, 
                         cuFloatComplex *a, magma_int_t *lda, float *d__, float *e, 
                         cuFloatComplex *tauq, cuFloatComplex *taup,
                         cuFloatComplex *x, magma_int_t *ldx,
                         cuFloatComplex *y, magma_int_t *ldy);
void    lapackf77_cpotrf(const char *uplo, magma_int_t *n, 
                         cuFloatComplex *a, magma_int_t *lda, magma_int_t *info);
void    lapackf77_cpotrs(const char *uplo, magma_int_t *n, magma_int_t *nrhs,
                         cuFloatComplex *a, magma_int_t *lda,
                         cuFloatComplex *b, magma_int_t *ldb, magma_int_t *info);
void    lapackf77_cpotri(const char *uplo, magma_int_t *n, 
                         cuFloatComplex *a, magma_int_t *lda, magma_int_t *info);
void    lapackf77_clauum(const char *uplo, magma_int_t *n, 
                         cuFloatComplex *a, magma_int_t *lda, magma_int_t *info);
void    lapackf77_ctrevc(const char *side, const char *howmny, magma_int_t *select, magma_int_t *n, 
                         cuFloatComplex *T,  magma_int_t *ldt,  cuFloatComplex *VL, magma_int_t *ldvl,
                         cuFloatComplex *VR, magma_int_t *ldvr, magma_int_t *MM, magma_int_t *M, 
                         cuFloatComplex *work, DWORKFORZ magma_int_t *info);
void    lapackf77_sstebz(const char *range, const char *order, magma_int_t *n, float *vl, float *vu,
                         magma_int_t *il, magma_int_t *iu, float *abstol,
                         float *d__, float *e, magma_int_t *m, magma_int_t *nsplit,
                         float *w, magma_int_t *iblock, magma_int_t *isplit, float *work,
                         magma_int_t *iwork, magma_int_t *info);
float  lapackf77_slamc3(float* a, float* b);
void    lapackf77_slamrg(magma_int_t* n1, magma_int_t* n2, float* a, 
                         magma_int_t* dtrd1, magma_int_t* dtrd2, magma_int_t* index);
void    lapackf77_slaed4(magma_int_t* n, magma_int_t* i, float* d, float* z,
                         float* delta, float* rho, float* dlam, magma_int_t* info);
void    lapackf77_csteqr(const char *compz, magma_int_t *n, float *D, float *E, 
                         cuFloatComplex *Z, magma_int_t *ldz, 
                         float *work, magma_int_t *info);
void    lapackf77_cstedc(const char *compz, magma_int_t *n, float *D, float *E, 
                         cuFloatComplex *Z, magma_int_t *ldz, 
                         cuFloatComplex *work, magma_int_t *ldwork, 
                         DWORKFORZ_AND_LD magma_int_t *iwork, magma_int_t *liwork,
                         magma_int_t *info);
void    lapackf77_cstein(magma_int_t *n, float *d__, float *e, 
                         magma_int_t *m, float *w, magma_int_t *iblock, magma_int_t *isplit, 
                         cuFloatComplex *z__, magma_int_t *ldz, float *work, magma_int_t *iwork, 
                         magma_int_t *ifail, magma_int_t *info);
void    lapackf77_cstemr(const char *jobz, const char *range, magma_int_t *n, float *d__, float *e, 
                         float *vl, float *vu, magma_int_t *il, magma_int_t *iu, magma_int_t *m,
                         float *w, cuFloatComplex *z__, magma_int_t *ldz, magma_int_t *nzc, 
                         magma_int_t *isuppz, magma_int_t *tryrac, float *work, magma_int_t *lwork, 
                         magma_int_t *iwork, magma_int_t *liwork, magma_int_t *info);
void    lapackf77_ctrtri(const char *uplo, const char *diag, magma_int_t *n,
                         cuFloatComplex *a, magma_int_t *lda, magma_int_t *info);
#if defined(PRECISION_z) || defined(PRECISION_c)
void    lapackf77_csymv(const char *uplo, const magma_int_t *N, const cuFloatComplex *alpha, 
                        const cuFloatComplex *A, const magma_int_t *lda, 
                        const cuFloatComplex *X, const magma_int_t *incX,
                        const cuFloatComplex *beta, 
                        cuFloatComplex *Y, const magma_int_t *incY);
#endif
void    lapackf77_cung2r(magma_int_t *m, magma_int_t *n, magma_int_t *k, 
                         cuFloatComplex *a, magma_int_t *lda,
                         const cuFloatComplex *tau, cuFloatComplex *work,
                         magma_int_t *info);
void    lapackf77_cungbr(const char *vect, magma_int_t *m, magma_int_t *n, magma_int_t *k,
                         cuFloatComplex *a, magma_int_t *lda, const cuFloatComplex *tau,
                         cuFloatComplex *work, magma_int_t *lwork, magma_int_t *info);
void    lapackf77_cunghr(magma_int_t *n, magma_int_t *ilo, magma_int_t *ihi, 
                         cuFloatComplex *a, magma_int_t *lda, const cuFloatComplex *tau,
                         cuFloatComplex *work, magma_int_t *lwork, magma_int_t *info);
void    lapackf77_cunglq(magma_int_t *m, magma_int_t *n, magma_int_t *k, 
                         cuFloatComplex *a, magma_int_t *lda, const cuFloatComplex *tau, 
                         cuFloatComplex *work, magma_int_t *ldwork, magma_int_t *info);
void    lapackf77_cungql(magma_int_t *, magma_int_t *, magma_int_t *,
                         cuFloatComplex *, magma_int_t *, cuFloatComplex *, 
                         cuFloatComplex *, magma_int_t *, magma_int_t *);
void    lapackf77_cungqr(magma_int_t *m, magma_int_t *n, magma_int_t *k, 
                         cuFloatComplex *a, magma_int_t *lda, const cuFloatComplex *tau, 
                         cuFloatComplex *work, magma_int_t *ldwork, magma_int_t *info);
void    lapackf77_cungtr(const char *uplo, magma_int_t *n, 
                         cuFloatComplex *a, magma_int_t *lda, const cuFloatComplex *tau, 
                         cuFloatComplex *work, magma_int_t *lwork, magma_int_t *info);
void    lapackf77_cunm2r(const char *side, const char *trans, 
                         magma_int_t *m, magma_int_t *n, magma_int_t *k, 
                         const cuFloatComplex *a, magma_int_t *lda, 
                         const cuFloatComplex *tau, cuFloatComplex *c, magma_int_t *ldc,
                         cuFloatComplex *work, magma_int_t *info);
void    lapackf77_cunmbr(const char *vect, const char *side, const char *trans,
                         magma_int_t *M, magma_int_t *N, magma_int_t *K, 
                         cuFloatComplex *A, magma_int_t *lda, cuFloatComplex *Tau,
                         cuFloatComplex *C, magma_int_t *ldc, 
                         cuFloatComplex *work, magma_int_t *ldwork, magma_int_t *info);
void    lapackf77_cunmlq(const char *side, const char *trans, 
                         magma_int_t *m, magma_int_t *n, magma_int_t *k,
                         const cuFloatComplex *a, magma_int_t *lda, 
                         const cuFloatComplex *tau, cuFloatComplex *c, magma_int_t *ldc, 
                         cuFloatComplex *work, magma_int_t *lwork, magma_int_t *info);
void    lapackf77_cunmql(const char *side, const char *trans, 
                         magma_int_t *m, magma_int_t *n, magma_int_t *k,
                         const cuFloatComplex *a, magma_int_t *lda, 
                         const cuFloatComplex *tau, cuFloatComplex *c, magma_int_t *ldc,
                         cuFloatComplex *work, magma_int_t *lwork, magma_int_t *info);
void    lapackf77_cunmqr(const char *side, const char *trans, 
                         magma_int_t *m, magma_int_t *n, magma_int_t *k, 
                         const cuFloatComplex *a, magma_int_t *lda, 
                         const cuFloatComplex *tau, cuFloatComplex *c, magma_int_t *ldc, 
                         cuFloatComplex *work, magma_int_t *lwork, magma_int_t *info);
void    lapackf77_cunmtr(const char *side, const char *uplo, const char *trans,
                         magma_int_t *M, magma_int_t *N,
                         cuFloatComplex *A, magma_int_t *lda, cuFloatComplex *Tau,
                         cuFloatComplex *C, magma_int_t *ldc, 
                         cuFloatComplex *work, magma_int_t *ldwork, magma_int_t *info);


  /*
   * Testing functions
   */

#if defined(PRECISION_z) || defined(PRECISION_c)

void    lapackf77_cbdt01(int *m, int *n, int *kd, cuFloatComplex *A, int *lda, 
                         cuFloatComplex *Q, int *ldq, float *D, float *E, 
                         cuFloatComplex *PT, int *ldpt, cuFloatComplex *work, 
                         float *rwork, float *resid);
void    lapackf77_cget22(const char *transa, const char *transe, const char *transw, int *n,
                         cuFloatComplex *a, int *lda, cuFloatComplex *e, int *lde,
                         cuFloatComplex *w, cuFloatComplex *work,
                         float *rwork, float *result);
void    lapackf77_chet21(int *itype, const char *uplo, int *n, int *kband, 
                         cuFloatComplex *A, int *lda, float *D, float *E, 
                         cuFloatComplex *U, int *ldu, cuFloatComplex *V, int *ldv, 
                         cuFloatComplex *TAU, cuFloatComplex *work,
                         float *rwork, float *result);
void    lapackf77_chst01(int *n, int *ilo, int *ihi, cuFloatComplex *A, int *lda, 
                         cuFloatComplex *H, int *ldh, cuFloatComplex *Q, int *ldq,
                         cuFloatComplex *work, int *lwork, float *rwork, float *result);
void    lapackf77_cstt21(int *n, int *kband, float *AD, float *AE, float *SD,
                         float *SE, cuFloatComplex *U, int *ldu, 
                         cuFloatComplex *work, float *rwork, float *result);
void    lapackf77_cunt01(const char *rowcol, int *m, int *n, cuFloatComplex *U, int *ldu,
                         cuFloatComplex *work, int *lwork, float *rwork, float *resid);

#else

void    lapackf77_cbdt01(int *m, int *n, int *kd, cuFloatComplex *A, int *lda, 
                         cuFloatComplex *Q, int *ldq, float *D, float *E, 
                         cuFloatComplex *PT, int *ldpt, 
                         cuFloatComplex *work, float *resid);
void    lapackf77_cget22(const char *transa, const char *transe, const char *transw, int *n,
                         cuFloatComplex *a, int *lda, cuFloatComplex *e, int *lde,
                         cuFloatComplex *wr, cuFloatComplex *wi, 
                         float *work, float *result);
void    lapackf77_chet21(int *itype, const char *uplo, int *n, int *kband, 
                         cuFloatComplex *A, int *lda, float *D, float *E,
                         cuFloatComplex *U, int *ldu, cuFloatComplex *V, int *ldv, 
                         cuFloatComplex *TAU, cuFloatComplex *work, float *result);
void    lapackf77_chst01(int *n, int *ilo, int *ihi, cuFloatComplex *A, int *lda, 
                         cuFloatComplex *H, int *ldh, cuFloatComplex *Q, int *ldq, 
                         cuFloatComplex *work, int *lwork, float *result);
void    lapackf77_cstt21(int *n, int *kband, float *AD, float *AE, float *SD, 
                         float *SE, cuFloatComplex *U, int *ldu, 
                         cuFloatComplex *work, float *result);
void    lapackf77_cunt01(const char *rowcol, int *m, int *n, cuFloatComplex *U, int *ldu,
                         cuFloatComplex *work, int *lwork, float *resid);
#endif

void    lapackf77_clarfy(const char *uplo, int *N, cuFloatComplex *V, int *incv, 
                         cuFloatComplex *tau, cuFloatComplex *C, int *ldc, 
                         cuFloatComplex *work);
void    lapackf77_clarfx(const char *, int *, int *, 
                         cuFloatComplex *, cuFloatComplex *, 
                         cuFloatComplex *, int *, cuFloatComplex *);
float  lapackf77_cqpt01(int *m, int *n, int *k, cuFloatComplex *a,
                         cuFloatComplex *af, int *lda, cuFloatComplex *tau, int *jpvt,
                         cuFloatComplex *work, int *lwork);
void    lapackf77_cqrt02(int *m, int *n, int *k, cuFloatComplex *A, cuFloatComplex *AF,
                         cuFloatComplex *Q, cuFloatComplex *R, int *lda, 
                         cuFloatComplex *TAU, cuFloatComplex *work, int *lwork,
                         float *rwork, float *result);

#ifdef __cplusplus
}
#endif

#undef DWORKFORZ 
#undef DWORKFORZ_AND_LD
#undef WSPLIT
#undef PRECISION_c
#endif /* MAGMA ZLAPACK */
