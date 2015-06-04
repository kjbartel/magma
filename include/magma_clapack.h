/*
 *   -- MAGMA (version 1.3.0) --
 *      Univ. of Tennessee, Knoxville
 *      Univ. of California, Berkeley
 *      Univ. of Colorado, Denver
 *      November 2012
 *
 * @generated c Wed Nov 14 22:52:27 2012
 */

#ifndef MAGMA_CLAPACK_H
#define MAGMA_CLAPACK_H

#define PRECISION_c

#ifdef __cplusplus
extern "C" {
#endif

/* ////////////////////////////////////////////////////////////////////////////
   -- BLAS and LAPACK functions (alphabetical order)
*/
#define blasf77_icamax     FORTRAN_NAME( icamax, ICAMAX )
#define blasf77_caxpy      FORTRAN_NAME( caxpy,  CAXPY  )
#define blasf77_ccopy      FORTRAN_NAME( ccopy,  CCOPY  )
#define blasf77_cgemm      FORTRAN_NAME( cgemm,  CGEMM  )
#define blasf77_cgemv      FORTRAN_NAME( cgemv,  CGEMV  )
#define blasf77_cgerc      FORTRAN_NAME( cgerc,  CGERC  )
#define blasf77_cgeru      FORTRAN_NAME( cgeru,  CGERU  )
#define blasf77_chemm      FORTRAN_NAME( chemm,  CHEMM  )
#define blasf77_chemv      FORTRAN_NAME( chemv,  CHEMV  )
#define blasf77_cher2      FORTRAN_NAME( cher2,  CHER2  )
#define blasf77_cher2k     FORTRAN_NAME( cher2k, CHER2K )
#define blasf77_cherk      FORTRAN_NAME( cherk,  CHERK  )
#define blasf77_cscal      FORTRAN_NAME( cscal,  CSCAL  )
#define blasf77_csscal     FORTRAN_NAME( csscal, CSSCAL )
#define blasf77_cswap      FORTRAN_NAME( cswap,  CSWAP  )
#define blasf77_csymm      FORTRAN_NAME( csymm,  CSYMM  )
#define blasf77_csyr2k     FORTRAN_NAME( csyr2k, CSYR2K )
#define blasf77_csyrk      FORTRAN_NAME( csyrk,  CSYRK  )
#define blasf77_ctrmm      FORTRAN_NAME( ctrmm,  CTRMM  )
#define blasf77_ctrmv      FORTRAN_NAME( ctrmv,  CTRMV  )
#define blasf77_ctrsm      FORTRAN_NAME( ctrsm,  CTRSM  )
#define blasf77_ctrsv      FORTRAN_NAME( ctrsv,  CTRSV  )

#define lapackf77_slaed4   FORTRAN_NAME( slaed4, SLAED4 )
#define lapackf77_slamc3   FORTRAN_NAME( slamc3, SLAMC3 )
#define lapackf77_slamrg   FORTRAN_NAME( slamrg, SLAMRG )
#define lapackf77_sstebz   FORTRAN_NAME( sstebz, SSTEBZ )

#define lapackf77_cbdsqr   FORTRAN_NAME( cbdsqr, CBDSQR )
#define lapackf77_cgebak   FORTRAN_NAME( cgebak, CGEBAK )
#define lapackf77_cgebal   FORTRAN_NAME( cgebal, CGEBAL )
#define lapackf77_cgebd2   FORTRAN_NAME( cgebd2, CGEBD2 )
#define lapackf77_cgebrd   FORTRAN_NAME( cgebrd, CGEBRD )
#define lapackf77_cgeev    FORTRAN_NAME( cgeev,  CGEEV  )
#define lapackf77_cgehd2   FORTRAN_NAME( cgehd2, CGEHD2 )
#define lapackf77_cgehrd   FORTRAN_NAME( cgehrd, CGEHRD )
#define lapackf77_cgelqf   FORTRAN_NAME( cgelqf, CGELQF )
#define lapackf77_cgels    FORTRAN_NAME( cgels,  CGELS  )
#define lapackf77_cgeqlf   FORTRAN_NAME( cgeqlf, CGEQLF )
#define lapackf77_cgeqp3   FORTRAN_NAME( cgeqp3, CGEQP3 )
#define lapackf77_cgeqrf   FORTRAN_NAME( cgeqrf, CGEQRF )
#define lapackf77_cgesvd   FORTRAN_NAME( cgesvd, CGESVD )
#define lapackf77_cgetrf   FORTRAN_NAME( cgetrf, CGETRF )
#define lapackf77_cgetri   FORTRAN_NAME( cgetri, CGETRI )
#define lapackf77_cgetrs   FORTRAN_NAME( cgetrs, CGETRS )
#define lapackf77_chbtrd   FORTRAN_NAME( chbtrd, CHBTRD )
#define lapackf77_cheev    FORTRAN_NAME( cheev,  CHEEV  )
#define lapackf77_cheevd   FORTRAN_NAME( cheevd, CHEEVD )
#define lapackf77_chegs2   FORTRAN_NAME( chegs2, CHEGS2 )
#define lapackf77_chegvd   FORTRAN_NAME( chegvd, CHEGVD )
#define lapackf77_chetd2   FORTRAN_NAME( chetd2, CHETD2 )
#define lapackf77_chetrd   FORTRAN_NAME( chetrd, CHETRD )
#define lapackf77_chetrf   FORTRAN_NAME( chetrf, CHETRF )
#define lapackf77_chseqr   FORTRAN_NAME( chseqr, CHSEQR )
#define lapackf77_clabrd   FORTRAN_NAME( clabrd, CLABRD )
#define lapackf77_cladiv   FORTRAN_NAME( cladiv, ZLADIV )
#define lapackf77_clacgv   FORTRAN_NAME( clacgv, CLACGV )
#define lapackf77_clacpy   FORTRAN_NAME( clacpy, CLACPY )
#define lapackf77_clahef   FORTRAN_NAME( clahef, CLAHEF )
#define lapackf77_clange   FORTRAN_NAME( clange, CLANGE )
#define lapackf77_clanhe   FORTRAN_NAME( clanhe, CLANHE )
#define lapackf77_clanht   FORTRAN_NAME( clanht, CLANHT )
#define lapackf77_clansy   FORTRAN_NAME( clansy, CLANSY )
#define lapackf77_slapy3   FORTRAN_NAME( slapy3, DLAPY3 )
#define lapackf77_claqp2   FORTRAN_NAME( claqp2, CLAQP2 )
#define lapackf77_clarf    FORTRAN_NAME( clarf,  ZLARF  )
#define lapackf77_clarfb   FORTRAN_NAME( clarfb, CLARFB )
#define lapackf77_clarfg   FORTRAN_NAME( clarfg, CLARFG )
#define lapackf77_clarft   FORTRAN_NAME( clarft, CLARFT )
#define lapackf77_clarnv   FORTRAN_NAME( clarnv, CLARNV )
#define lapackf77_clartg   FORTRAN_NAME( clartg, CLARTG )
#define lapackf77_clascl   FORTRAN_NAME( clascl, CLASCL )
#define lapackf77_claset   FORTRAN_NAME( claset, CLASET )
#define lapackf77_claswp   FORTRAN_NAME( claswp, CLASWP )
#define lapackf77_clatrd   FORTRAN_NAME( clatrd, CLATRD )
#define lapackf77_clauum   FORTRAN_NAME( clauum, CLAUUM )
#define lapackf77_clavhe   FORTRAN_NAME( clavhe, CLAVHE )
#define lapackf77_cpotrf   FORTRAN_NAME( cpotrf, CPOTRF )
#define lapackf77_cpotri   FORTRAN_NAME( cpotri, CPOTRI )
#define lapackf77_cpotrs   FORTRAN_NAME( cpotrs, CPOTRS )
#define lapackf77_cstedc   FORTRAN_NAME( cstedc, CSTEDC )
#define lapackf77_cstein   FORTRAN_NAME( cstein, CSTEIN )
#define lapackf77_cstemr   FORTRAN_NAME( cstemr, CSTEMR )
#define lapackf77_csteqr   FORTRAN_NAME( csteqr, CSTEQR )
#define lapackf77_csymv    FORTRAN_NAME( csymv,  CSYMV  )
#define lapackf77_ctrevc   FORTRAN_NAME( ctrevc, CTREVC )
#define lapackf77_ctrtri   FORTRAN_NAME( ctrtri, CTRTRI )
#define lapackf77_cung2r   FORTRAN_NAME( cung2r, CUNG2R )
#define lapackf77_cungbr   FORTRAN_NAME( cungbr, CUNGBR )
#define lapackf77_cunghr   FORTRAN_NAME( cunghr, CUNGHR )
#define lapackf77_cunglq   FORTRAN_NAME( cunglq, CUNGLQ )
#define lapackf77_cungql   FORTRAN_NAME( cungql, CUNGQL )
#define lapackf77_cungqr   FORTRAN_NAME( cungqr, CUNGQR )
#define lapackf77_cungtr   FORTRAN_NAME( cungtr, CUNGTR )
#define lapackf77_cunm2r   FORTRAN_NAME( cunm2r, CUNM2R )
#define lapackf77_cunmbr   FORTRAN_NAME( cunmbr, CUNMBR )
#define lapackf77_cunmlq   FORTRAN_NAME( cunmlq, CUNMLQ )
#define lapackf77_cunmql   FORTRAN_NAME( cunmql, CUNMQL )
#define lapackf77_cunmqr   FORTRAN_NAME( cunmqr, CUNMQR )
#define lapackf77_cunmtr   FORTRAN_NAME( cunmtr, CUNMTR )

/* testing functions */
#define lapackf77_cbdt01   FORTRAN_NAME( cbdt01, CBDT01 )
#define lapackf77_cget22   FORTRAN_NAME( cget22, CGET22 )
#define lapackf77_cqpt01   FORTRAN_NAME( cqpt01, CQPT01 )
#define lapackf77_chet21   FORTRAN_NAME( chet21, CHET21 )
#define lapackf77_chst01   FORTRAN_NAME( chst01, CHST01 )
#define lapackf77_cqrt02   FORTRAN_NAME( cqrt02, CQRT02 )
#define lapackf77_cunt01   FORTRAN_NAME( cunt01, CUNT01 )
#define lapackf77_clarfy   FORTRAN_NAME( clarfy, CLARFY )
#define lapackf77_clarfx   FORTRAN_NAME( clarfx, CLARFX )
#define lapackf77_cstt21   FORTRAN_NAME( cstt21, CSTT21 )

// macros to handle differences in arguments between complex and real versions of routines.
#if defined(PRECISION_z) || defined(PRECISION_c)
#define DWORKFORZ        float *rwork,
#define DWORKFORZ_AND_LD float *rwork, const magma_int_t *ldrwork,
#define WSPLIT           cuFloatComplex *w
#else
#define DWORKFORZ
#define DWORKFORZ_AND_LD
#define WSPLIT           float *wr, float *wi
#endif

/*
 * BLAS functions (alphabetical order)
 */
magma_int_t blasf77_icamax(
                     const magma_int_t *n,
                     const cuFloatComplex *x, const magma_int_t *incx);

void blasf77_caxpy(  const magma_int_t *n,
                     const cuFloatComplex *alpha,
                     const cuFloatComplex *x, const magma_int_t *incx,
                           cuFloatComplex *y, const magma_int_t *incy );

void blasf77_ccopy(  const magma_int_t *n,
                     const cuFloatComplex *x, const magma_int_t *incx,
                           cuFloatComplex *y, const magma_int_t *incy );

void blasf77_cgemm(  const char *transa, const char *transb,
                     const magma_int_t *m, const magma_int_t *n, const magma_int_t *k,
                     const cuFloatComplex *alpha,
                     const cuFloatComplex *A, const magma_int_t *lda,
                     const cuFloatComplex *B, const magma_int_t *ldb,
                     const cuFloatComplex *beta,
                           cuFloatComplex *C, const magma_int_t *ldc );

void blasf77_cgemv(  const char *transa,
                     const magma_int_t *m, const magma_int_t *n,
                     const cuFloatComplex *alpha,
                     const cuFloatComplex *A, const magma_int_t *lda,
                     const cuFloatComplex *x, const magma_int_t *incx,
                     const cuFloatComplex *beta,
                           cuFloatComplex *y, const magma_int_t *incy );

void blasf77_cgerc(  const magma_int_t *m, const magma_int_t *n,
                     const cuFloatComplex *alpha,
                     const cuFloatComplex *x, const magma_int_t *incx,
                     const cuFloatComplex *y, const magma_int_t *incy,
                           cuFloatComplex *A, const magma_int_t *lda );

#if defined(PRECISION_z) || defined(PRECISION_c)
void blasf77_cgeru(  const magma_int_t *m, const magma_int_t *n,
                     const cuFloatComplex *alpha,
                     const cuFloatComplex *x, const magma_int_t *incx,
                     const cuFloatComplex *y, const magma_int_t *incy,
                           cuFloatComplex *A, const magma_int_t *lda );
#endif

void blasf77_chemm(  const char *side, const char *uplo,
                     const magma_int_t *m, const magma_int_t *n,
                     const cuFloatComplex *alpha,
                     const cuFloatComplex *A, const magma_int_t *lda,
                     const cuFloatComplex *B, const magma_int_t *ldb,
                     const cuFloatComplex *beta,
                           cuFloatComplex *C, const magma_int_t *ldc );

void blasf77_chemv(  const char *uplo,
                     const magma_int_t *n,
                     const cuFloatComplex *alpha,
                     const cuFloatComplex *A, const magma_int_t *lda,
                     const cuFloatComplex *x, const magma_int_t *incx,
                     const cuFloatComplex *beta,
                           cuFloatComplex *y, const magma_int_t *incy );

void blasf77_cher2(  const char *uplo,
                     const magma_int_t *n,
                     const cuFloatComplex *alpha,
                     const cuFloatComplex *x, const magma_int_t *incx,
                     const cuFloatComplex *y, const magma_int_t *incy,
                           cuFloatComplex *A, const magma_int_t *lda );

void blasf77_cher2k(  const char *uplo, const char *trans,
                     const magma_int_t *n, const magma_int_t *k,
                     const cuFloatComplex *alpha,
                     const cuFloatComplex *A, const magma_int_t *lda,
                     const cuFloatComplex *B, const magma_int_t *ldb,
                     const float *beta,
                           cuFloatComplex *C, const magma_int_t *ldc );

void blasf77_cherk(  const char *uplo, const char *trans,
                     const magma_int_t *n, const magma_int_t *k,
                     const float *alpha,
                     const cuFloatComplex *A, const magma_int_t *lda,
                     const float *beta,
                           cuFloatComplex *C, const magma_int_t *ldc );

void blasf77_cscal(  const magma_int_t *n,
                     const cuFloatComplex *alpha,
                           cuFloatComplex *x, const magma_int_t *incx );

#if defined(PRECISION_z) || defined(PRECISION_c)
void blasf77_csscal( const magma_int_t *n,
                     const float *alpha,
                           cuFloatComplex *x, const magma_int_t *incx );
#endif

void blasf77_cswap(  const magma_int_t *n,
                     cuFloatComplex *x, const magma_int_t *incx,
                     cuFloatComplex *y, const magma_int_t *incy );

void blasf77_csymm(  const char *side, const char *uplo,
                     const magma_int_t *m, const magma_int_t *n,
                     const cuFloatComplex *alpha,
                     const cuFloatComplex *A, const magma_int_t *lda,
                     const cuFloatComplex *B, const magma_int_t *ldb,
                     const cuFloatComplex *beta,
                           cuFloatComplex *C, const magma_int_t *ldc );

void blasf77_csyr2k( const char *uplo, const char *trans,
                     const magma_int_t *n, const magma_int_t *k,
                     const cuFloatComplex *alpha,
                     const cuFloatComplex *A, const magma_int_t *lda,
                     const cuFloatComplex *B, const magma_int_t *ldb,
                     const cuFloatComplex *beta,
                           cuFloatComplex *C, const magma_int_t *ldc );

void blasf77_csyrk(  const char *uplo, const char *trans,
                     const magma_int_t *n, const magma_int_t *k,
                     const cuFloatComplex *alpha,
                     const cuFloatComplex *A, const magma_int_t *lda,
                     const cuFloatComplex *beta,
                           cuFloatComplex *C, const magma_int_t *ldc );

void blasf77_ctrmm(  const char *side, const char *uplo, const char *transa, const char *diag,
                     const magma_int_t *m, const magma_int_t *n,
                     const cuFloatComplex *alpha,
                     const cuFloatComplex *A, const magma_int_t *lda,
                           cuFloatComplex *B, const magma_int_t *ldb );

void blasf77_ctrmv(  const char *uplo, const char *transa, const char *diag,
                     const magma_int_t *n,
                     const cuFloatComplex *A, const magma_int_t *lda,
                           cuFloatComplex *x, const magma_int_t *incx );

void blasf77_ctrsm(  const char *side, const char *uplo, const char *transa, const char *diag,
                     const magma_int_t *m, const magma_int_t *n,
                     const cuFloatComplex *alpha,
                     const cuFloatComplex *A, const magma_int_t *lda,
                           cuFloatComplex *B, const magma_int_t *ldb );

void blasf77_ctrsv(  const char *uplo, const char *transa, const char *diag,
                     const magma_int_t *n,
                     const cuFloatComplex *A, const magma_int_t *lda,
                           cuFloatComplex *x, const magma_int_t *incx );

/*
 * LAPACK functions (alphabetical order)
 */
void   lapackf77_cbdsqr( const char *uplo,
                         const magma_int_t *n, const magma_int_t *ncvt, const magma_int_t *nru,  const magma_int_t *ncc,
                         float *d, float *e,
                         cuFloatComplex *Vt, const magma_int_t *ldvt,
                         cuFloatComplex *U, const magma_int_t *ldu,
                         cuFloatComplex *C, const magma_int_t *ldc,
                         float *work,
                         magma_int_t *info );

void   lapackf77_cgebak( const char *job, const char *side,
                         const magma_int_t *n,
                         const magma_int_t *ilo, const magma_int_t *ihi,
                         const float *scale, const magma_int_t *m,
                         cuFloatComplex *V, const magma_int_t *ldv,
                         magma_int_t *info );

void   lapackf77_cgebal( const char *job,
                         const magma_int_t *n,
                         cuFloatComplex *A, const magma_int_t *lda,
                         magma_int_t *ilo, magma_int_t *ihi,
                         float *scale,
                         magma_int_t *info );

void   lapackf77_cgebd2( const magma_int_t *m, const magma_int_t *n,
                         cuFloatComplex *A, const magma_int_t *lda,
                         float *d, float *e,
                         cuFloatComplex *tauq,
                         cuFloatComplex *taup,
                         cuFloatComplex *work,
                         magma_int_t *info );

void   lapackf77_cgebrd( const magma_int_t *m, const magma_int_t *n,
                         cuFloatComplex *A, const magma_int_t *lda,
                         float *d, float *e,
                         cuFloatComplex *tauq,
                         cuFloatComplex *taup,
                         cuFloatComplex *work, const magma_int_t *lwork,
                         magma_int_t *info );

void   lapackf77_cgeev(  const char *jobvl, const char *jobvr,
                         const magma_int_t *n,
                         cuFloatComplex *A,    const magma_int_t *lda,
                         WSPLIT,
                         cuFloatComplex *Vl,   const magma_int_t *ldvl,
                         cuFloatComplex *Vr,   const magma_int_t *ldvr,
                         cuFloatComplex *work, const magma_int_t *lwork,
                         DWORKFORZ
                         magma_int_t *info );

void   lapackf77_cgehd2( const magma_int_t *n,
                         const magma_int_t *ilo, const magma_int_t *ihi,
                         cuFloatComplex *A, const magma_int_t *lda,
                         cuFloatComplex *tau,
                         cuFloatComplex *work,
                         magma_int_t *info );

void   lapackf77_cgehrd( const magma_int_t *n,
                         const magma_int_t *ilo, const magma_int_t *ihi,
                         cuFloatComplex *A, const magma_int_t *lda,
                         cuFloatComplex *tau,
                         cuFloatComplex *work, const magma_int_t *lwork,
                         magma_int_t *info );

void   lapackf77_cgelqf( const magma_int_t *m, const magma_int_t *n,
                         cuFloatComplex *A, const magma_int_t *lda,
                         cuFloatComplex *tau,
                         cuFloatComplex *work, const magma_int_t *lwork,
                         magma_int_t *info );

void   lapackf77_cgels(  const char *trans,
                         const magma_int_t *m, const magma_int_t *n, const magma_int_t *nrhs,
                         cuFloatComplex *A, const magma_int_t *lda,
                         cuFloatComplex *B, const magma_int_t *ldb,
                         cuFloatComplex *work, const magma_int_t *lwork,
                         magma_int_t *info );

void   lapackf77_cgeqlf( const magma_int_t *m, const magma_int_t *n,
                         cuFloatComplex *A, const magma_int_t *lda,
                         cuFloatComplex *tau,
                         cuFloatComplex *work, const magma_int_t *lwork,
                         magma_int_t *info );

void   lapackf77_cgeqp3( const magma_int_t *m, const magma_int_t *n,
                         cuFloatComplex *A, const magma_int_t *lda,
                         magma_int_t *jpvt,
                         cuFloatComplex *tau,
                         cuFloatComplex *work, const magma_int_t *lwork,
                         DWORKFORZ
                         magma_int_t *info );

void   lapackf77_cgeqrf( const magma_int_t *m, const magma_int_t *n,
                         cuFloatComplex *A, const magma_int_t *lda,
                         cuFloatComplex *tau,
                         cuFloatComplex *work, const magma_int_t *lwork,
                         magma_int_t *info );

void   lapackf77_cgesvd( const char *jobu, const char *jobvt,
                         const magma_int_t *m, const magma_int_t *n,
                         cuFloatComplex *A, const magma_int_t *lda,
                         float *s,
                         cuFloatComplex *U,  const magma_int_t *ldu,
                         cuFloatComplex *Vt, const magma_int_t *ldvt,
                         cuFloatComplex *work, const magma_int_t *lwork,
                         DWORKFORZ
                         magma_int_t *info );

void   lapackf77_cgetrf( const magma_int_t *m, const magma_int_t *n,
                         cuFloatComplex *A, const magma_int_t *lda,
                         magma_int_t *ipiv,
                         magma_int_t *info );

void   lapackf77_cgetri( const magma_int_t *n,
                         cuFloatComplex *A, const magma_int_t *lda,
                         const magma_int_t *ipiv,
                         cuFloatComplex *work, const magma_int_t *lwork,
                         magma_int_t *info );

void   lapackf77_cgetrs( const char* trans,
                         const magma_int_t *n, const magma_int_t *nrhs,
                         const cuFloatComplex *A, const magma_int_t *lda,
                         const magma_int_t *ipiv,
                         cuFloatComplex *B, const magma_int_t *ldb,
                         magma_int_t *info );

void   lapackf77_chbtrd( const char *vect, const char *uplo,
                         const magma_int_t *n, const magma_int_t *kd,
                         cuFloatComplex *Ab, const magma_int_t *ldab,
                         float *d, float *e,
                         cuFloatComplex *Q, const magma_int_t *ldq,
                         cuFloatComplex *work,
                         magma_int_t *info );

void   lapackf77_cheev(  const char *jobz, const char *uplo,
                         const magma_int_t *n,
                         cuFloatComplex *A, const magma_int_t *lda,
                         float *w,
                         cuFloatComplex *work, const magma_int_t *lwork,
                         DWORKFORZ
                         magma_int_t *info );

void   lapackf77_cheevd( const char *jobz, const char *uplo,
                         const magma_int_t *n,
                         cuFloatComplex *A, const magma_int_t *lda,
                         float *w,
                         cuFloatComplex *work, const magma_int_t *lwork,
                         DWORKFORZ_AND_LD
                         magma_int_t *iwork, const magma_int_t *liwork,
                         magma_int_t *info );

void   lapackf77_chegs2( const magma_int_t *itype, const char *uplo,
                         const magma_int_t *n,
                         cuFloatComplex *A, const magma_int_t *lda,
                         cuFloatComplex *B, const magma_int_t *ldb,
                         magma_int_t *info );

void   lapackf77_chegvd( const magma_int_t *itype, const char *jobz, const char *uplo,
                         const magma_int_t *n,
                         cuFloatComplex *A, const magma_int_t *lda,
                         cuFloatComplex *B, const magma_int_t *ldb,
                         float *w,
                         cuFloatComplex *work, const magma_int_t *lwork,
                         DWORKFORZ_AND_LD
                         magma_int_t *iwork, const magma_int_t *liwork,
                         magma_int_t *info );

void   lapackf77_chetd2( const char *uplo,
                         const magma_int_t *n,
                         cuFloatComplex *A, const magma_int_t *lda,
                         float *d, float *e,
                         cuFloatComplex *tau,
                         magma_int_t *info );

void   lapackf77_chetrd( const char *uplo,
                         const magma_int_t *n,
                         cuFloatComplex *A, const magma_int_t *lda,
                         float *d, float *e,
                         cuFloatComplex *tau,
                         cuFloatComplex *work, const magma_int_t *lwork,
                         magma_int_t *info );

void   lapackf77_chetrf( const char *uplo,
                         const magma_int_t *n,
                         cuFloatComplex *A, const magma_int_t *lda,
                         magma_int_t *ipiv,
                         cuFloatComplex *work, const magma_int_t *lwork,
                         magma_int_t *info );

void   lapackf77_chseqr( const char *job, const char *compz,
                         const magma_int_t *n,
                         const magma_int_t *ilo, const magma_int_t *ihi,
                         cuFloatComplex *H, const magma_int_t *ldh,
                         WSPLIT,
                         cuFloatComplex *Z, const magma_int_t *ldz,
                         cuFloatComplex *work, const magma_int_t *lwork,
                         magma_int_t *info );

void   lapackf77_clabrd( const magma_int_t *m, const magma_int_t *n, const magma_int_t *nb,
                         cuFloatComplex *A, const magma_int_t *lda,
                         float *d, float *e,
                         cuFloatComplex *tauq,
                         cuFloatComplex *taup,
                         cuFloatComplex *X, const magma_int_t *ldx,
                         cuFloatComplex *Y, const magma_int_t *ldy );

void   lapackf77_cladiv( cuFloatComplex *ret_val, cuFloatComplex *x, 
                         cuFloatComplex *y );

void   lapackf77_clacgv( const magma_int_t *n,
                         cuFloatComplex *x, const magma_int_t *incx );

void   lapackf77_clacpy( const char *uplo,
                         const magma_int_t *m, const magma_int_t *n,
                         const cuFloatComplex *A, const magma_int_t *lda,
                         cuFloatComplex *B, const magma_int_t *ldb );

void   lapackf77_clahef( const char *uplo,
                         const magma_int_t *n, const magma_int_t *kn,
                         magma_int_t *kb,
                         cuFloatComplex *A, const magma_int_t lda,
                         magma_int_t *ipiv,
                         cuFloatComplex *work, const magma_int_t *ldwork,
                         magma_int_t *info );

float lapackf77_clange( const char *norm,
                         const magma_int_t *m, const magma_int_t *n,
                         const cuFloatComplex *A, const magma_int_t *lda,
                         float *work );

float lapackf77_clanhe( const char *norm, const char *uplo,
                         const magma_int_t *n,
                         const cuFloatComplex *A, const magma_int_t *lda,
                         float * work );

float lapackf77_clanht( const char* norm, const magma_int_t* n,
                         const float* d, const cuFloatComplex* e );

float lapackf77_clansy( const char *norm, const char *uplo,
                         const magma_int_t *n,
                         const cuFloatComplex *A, const magma_int_t *lda,
                         float * work );

void lapackf77_claqp2 (  magma_int_t *m, magma_int_t *n, magma_int_t *offset,
                         cuFloatComplex *a, magma_int_t *lda, magma_int_t *jpvt, 
                         cuFloatComplex *tau,
                         float *vn1, float *vn2, cuFloatComplex *work);

void lapackf77_clarf  (  char *, magma_int_t *, magma_int_t *,
                         cuFloatComplex *, magma_int_t *, cuFloatComplex *, cuFloatComplex *,
                         magma_int_t *, cuFloatComplex *);

void   lapackf77_clarfb( const char *side, const char *trans, const char *direct, const char *storev,
                         const magma_int_t *m, const magma_int_t *n, const magma_int_t *k,
                         const cuFloatComplex *V, const magma_int_t *ldv,
                         const cuFloatComplex *T, const magma_int_t *ldt,
                         cuFloatComplex *C, const magma_int_t *ldc,
                         cuFloatComplex *work, const magma_int_t *ldwork );

void   lapackf77_clarfg( const magma_int_t *n,
                         cuFloatComplex *alpha,
                         cuFloatComplex *x, const magma_int_t *incx,
                         cuFloatComplex *tau );

void   lapackf77_clarft( const char *direct, const char *storev,
                         const magma_int_t *n, const magma_int_t *k,
                         cuFloatComplex *V, const magma_int_t *ldv,
                         const cuFloatComplex *tau,
                         cuFloatComplex *T, const magma_int_t *ldt );

void   lapackf77_clarnv( const magma_int_t *idist, magma_int_t *iseed, const magma_int_t *n,
                         cuFloatComplex *x );

void   lapackf77_clartg( cuFloatComplex *F,
                         cuFloatComplex *G,
                         float *cs,
                         cuFloatComplex *SN,
                         cuFloatComplex *R );

void   lapackf77_clascl( const char *type,
                         const magma_int_t *kl, const magma_int_t *ku,
                         float *cfrom,
                         float *cto,
                         const magma_int_t *m, const magma_int_t *n,
                         cuFloatComplex *A, const magma_int_t *lda,
                         magma_int_t *info );

void   lapackf77_claset( const char *uplo,
                         const magma_int_t *m, const magma_int_t *n,
                         const cuFloatComplex *alpha,
                         const cuFloatComplex *beta,
                         cuFloatComplex *A, const magma_int_t *lda );

void   lapackf77_claswp( const magma_int_t *n,
                         cuFloatComplex *A, const magma_int_t *lda,
                         const magma_int_t *k1, const magma_int_t *k2,
                         magma_int_t *ipiv,
                         const magma_int_t *incx );

void   lapackf77_clatrd( const char *uplo,
                         const magma_int_t *n, const magma_int_t *nb,
                         cuFloatComplex *A, const magma_int_t *lda,
                         float *e,
                         cuFloatComplex *tau,
                         cuFloatComplex *work, const magma_int_t *ldwork );

void   lapackf77_clauum( const char *uplo,
                         const magma_int_t *n,
                         cuFloatComplex *A, const magma_int_t *lda,
                         magma_int_t *info );

void   lapackf77_clavhe( const char *uplo, const char *trans, const char *diag,
                         magma_int_t *n, magma_int_t *nrhs,
                         cuFloatComplex *A, magma_int_t *lda,
                         magma_int_t *ipiv,
                         cuFloatComplex *B, magma_int_t *ldb,
                         magma_int_t *info );

void   lapackf77_cpotrf( const char *uplo,
                         const magma_int_t *n,
                         cuFloatComplex *A, const magma_int_t *lda,
                         magma_int_t *info );

void   lapackf77_cpotri( const char *uplo,
                         const magma_int_t *n,
                         cuFloatComplex *A, const magma_int_t *lda,
                         magma_int_t *info );

void   lapackf77_cpotrs( const char *uplo,
                         const magma_int_t *n, const magma_int_t *nrhs,
                         const cuFloatComplex *A, const magma_int_t *lda,
                         cuFloatComplex *B, const magma_int_t *ldb,
                         magma_int_t *info );

void   lapackf77_cstedc( const char *compz,
                         const magma_int_t *n,
                         float *d, float *e,
                         cuFloatComplex *Z, const magma_int_t *ldz,
                         cuFloatComplex *work, const magma_int_t *lwork,
                         DWORKFORZ_AND_LD
                         magma_int_t *iwork, const magma_int_t *liwork,
                         magma_int_t *info );

void   lapackf77_cstein( const magma_int_t *n,
                         const float *d, const float *e,
                         const magma_int_t *m,
                         const float *w,
                         const magma_int_t *iblock,
                         const magma_int_t *isplit,
                         cuFloatComplex *Z, const magma_int_t *ldz,
                         float *work, magma_int_t *iwork, magma_int_t *ifailv,
                         magma_int_t *info );

void   lapackf77_cstemr( const char *jobz, const char *range,
                         const magma_int_t *n,
                         float *d, float *e,
                         const float *vl, const float *vu,
                         const magma_int_t *il, const magma_int_t *iu,
                         magma_int_t *m,
                         float *w,
                         cuFloatComplex *Z, const magma_int_t *ldz,
                         const magma_int_t *nzc, magma_int_t *isuppz, magma_int_t *tryrac,
                         float *work, const magma_int_t *lwork,
                         magma_int_t *iwork, const magma_int_t *liwork,
                         magma_int_t *info );

void   lapackf77_csteqr( const char *compz,
                         const magma_int_t *n,
                         float *d, float *e,
                         cuFloatComplex *Z, const magma_int_t *ldz,
                         float *work,
                         magma_int_t *info );

#if defined(PRECISION_z) || defined(PRECISION_c)
void   lapackf77_csymv(  const char *uplo,
                         const magma_int_t *n,
                         const cuFloatComplex *alpha,
                         const cuFloatComplex *A, const magma_int_t *lda,
                         const cuFloatComplex *x, const magma_int_t *incx,
                         const cuFloatComplex *beta,
                               cuFloatComplex *y, const magma_int_t *incy );
#endif

void   lapackf77_ctrevc( const char *side, const char *howmny,
                         magma_int_t *select, const magma_int_t *n,
                         cuFloatComplex *T,  const magma_int_t *ldt,
                         cuFloatComplex *Vl, const magma_int_t *ldvl,
                         cuFloatComplex *Vr, const magma_int_t *ldvr,
                         const magma_int_t *mm, magma_int_t *m,
                         cuFloatComplex *work,
                         DWORKFORZ
                         magma_int_t *info );

void   lapackf77_ctrtri( const char *uplo, const char *diag,
                         const magma_int_t *n,
                         cuFloatComplex *A, const magma_int_t *lda,
                         magma_int_t *info );

void   lapackf77_cung2r( const magma_int_t *m, const magma_int_t *n, const magma_int_t *k,
                         cuFloatComplex *A, const magma_int_t *lda,
                         const cuFloatComplex *tau,
                         cuFloatComplex *work,
                         magma_int_t *info );

void   lapackf77_cungbr( const char *vect,
                         const magma_int_t *m, const magma_int_t *n, const magma_int_t *k,
                         cuFloatComplex *A, const magma_int_t *lda,
                         const cuFloatComplex *tau,
                         cuFloatComplex *work, const magma_int_t *lwork,
                         magma_int_t *info );

void   lapackf77_cunghr( const magma_int_t *n,
                         const magma_int_t *ilo, const magma_int_t *ihi,
                         cuFloatComplex *A, const magma_int_t *lda,
                         const cuFloatComplex *tau,
                         cuFloatComplex *work, const magma_int_t *lwork,
                         magma_int_t *info );

void   lapackf77_cunglq( const magma_int_t *m, const magma_int_t *n, const magma_int_t *k,
                         cuFloatComplex *A, const magma_int_t *lda,
                         const cuFloatComplex *tau,
                         cuFloatComplex *work, const magma_int_t *lwork,
                         magma_int_t *info );

void   lapackf77_cungql( const magma_int_t *m, const magma_int_t *n, const magma_int_t *k,
                         cuFloatComplex *A, const magma_int_t *lda,
                         const cuFloatComplex *tau,
                         cuFloatComplex *work, const magma_int_t *lwork,
                         magma_int_t *info );

void   lapackf77_cungqr( const magma_int_t *m, const magma_int_t *n, const magma_int_t *k,
                         cuFloatComplex *A, const magma_int_t *lda,
                         const cuFloatComplex *tau,
                         cuFloatComplex *work, const magma_int_t *lwork,
                         magma_int_t *info );

void   lapackf77_cungtr( const char *uplo,
                         const magma_int_t *n,
                         cuFloatComplex *A, const magma_int_t *lda,
                         const cuFloatComplex *tau,
                         cuFloatComplex *work, const magma_int_t *lwork,
                         magma_int_t *info );

void   lapackf77_cunm2r( const char *side, const char *trans,
                         const magma_int_t *m, const magma_int_t *n, const magma_int_t *k,
                         const cuFloatComplex *A, const magma_int_t *lda,
                         const cuFloatComplex *tau,
                         cuFloatComplex *C, const magma_int_t *ldc,
                         cuFloatComplex *work,
                         magma_int_t *info );

void   lapackf77_cunmbr( const char *vect, const char *side, const char *trans,
                         const magma_int_t *m, const magma_int_t *n, const magma_int_t *k,
                         const cuFloatComplex *A, const magma_int_t *lda,
                         const cuFloatComplex *tau,
                         cuFloatComplex *C, const magma_int_t *ldc,
                         cuFloatComplex *work, const magma_int_t *lwork,
                         magma_int_t *info );

void   lapackf77_cunmlq( const char *side, const char *trans,
                         const magma_int_t *m, const magma_int_t *n, const magma_int_t *k,
                         const cuFloatComplex *A, const magma_int_t *lda,
                         const cuFloatComplex *tau,
                         cuFloatComplex *C, const magma_int_t *ldc,
                         cuFloatComplex *work, const magma_int_t *lwork,
                         magma_int_t *info );

void   lapackf77_cunmql( const char *side, const char *trans,
                         const magma_int_t *m, const magma_int_t *n, const magma_int_t *k,
                         const cuFloatComplex *A, const magma_int_t *lda,
                         const cuFloatComplex *tau,
                         cuFloatComplex *C, const magma_int_t *ldc,
                         cuFloatComplex *work, const magma_int_t *lwork,
                         magma_int_t *info );

void   lapackf77_cunmqr( const char *side, const char *trans,
                         const magma_int_t *m, const magma_int_t *n, const magma_int_t *k,
                         const cuFloatComplex *A, const magma_int_t *lda,
                         const cuFloatComplex *tau,
                         cuFloatComplex *C, const magma_int_t *ldc,
                         cuFloatComplex *work, const magma_int_t *lwork,
                         magma_int_t *info );

void   lapackf77_cunmtr( const char *side, const char *uplo, const char *trans,
                         const magma_int_t *m, const magma_int_t *n,
                         const cuFloatComplex *A, const magma_int_t *lda,
                         const cuFloatComplex *tau,
                         cuFloatComplex *C, const magma_int_t *ldc,
                         cuFloatComplex *work, const magma_int_t *lwork,
                         magma_int_t *info );

/*
 * Real precision extras
 */
void   lapackf77_sstebz( const char *range, const char *order,
                         const magma_int_t *n,
                         float *vl, float *vu,
                         magma_int_t *il, magma_int_t *iu,
                         float *abstol,
                         float *d, float *e,
                         const magma_int_t *m, const magma_int_t *nsplit,
                         float *w,
                         magma_int_t *iblock, magma_int_t *isplit,
                         float *work,
                         magma_int_t *iwork,
                         magma_int_t *info );

float lapackf77_slamc3( float* a, float* b );

void   lapackf77_slamrg( magma_int_t* n1, magma_int_t* n2,
                         float* a,
                         magma_int_t* dtrd1, magma_int_t* dtrd2, magma_int_t* index );

float lapackf77_slapy3(float *, float *, float *);

void   lapackf77_slaed4( magma_int_t* n, magma_int_t* i,
                         float* d,
                         float* z,
                         float* delta,
                         float* rho,
                         float* dlam, magma_int_t* info );

/*
 * Testing functions
 */
#if defined(PRECISION_z) || defined(PRECISION_c)
void   lapackf77_cbdt01( const magma_int_t *m, const magma_int_t *n, const magma_int_t *kd,
                         cuFloatComplex *A, const magma_int_t *lda,
                         cuFloatComplex *Q, const magma_int_t *ldq,
                         float *d, float *e,
                         cuFloatComplex *Pt, const magma_int_t *ldpt,
                         cuFloatComplex *work,
                         float *rwork,
                         float *resid );

void   lapackf77_cget22( const char *transa, const char *transe, const char *transw, const magma_int_t *n,
                         cuFloatComplex *A, const magma_int_t *lda,
                         cuFloatComplex *E, const magma_int_t *lde,
                         cuFloatComplex *w,
                         cuFloatComplex *work,
                         float *rwork,
                         float *result );

void   lapackf77_chet21( const magma_int_t *itype, const char *uplo,
                         const magma_int_t *n, const magma_int_t *kband,
                         cuFloatComplex *A, const magma_int_t *lda,
                         float *d, float *e,
                         cuFloatComplex *U, const magma_int_t *ldu,
                         cuFloatComplex *V, const magma_int_t *ldv,
                         cuFloatComplex *tau,
                         cuFloatComplex *work,
                         float *rwork,
                         float *result );

void   lapackf77_chst01( const magma_int_t *n, const magma_int_t *ilo, const magma_int_t *ihi,
                         cuFloatComplex *A, const magma_int_t *lda,
                         cuFloatComplex *H, const magma_int_t *ldh,
                         cuFloatComplex *Q, const magma_int_t *ldq,
                         cuFloatComplex *work, const magma_int_t *lwork,
                         float *rwork,
                         float *result );

void   lapackf77_cstt21( const magma_int_t *n, const magma_int_t *kband,
                         float *AD,
                         float *AE,
                         float *SD,
                         float *SE,
                         cuFloatComplex *U, const magma_int_t *ldu,
                         cuFloatComplex *work,
                         float *rwork,
                         float *result );

void   lapackf77_cunt01( const char *rowcol, const magma_int_t *m, const magma_int_t *n,
                         cuFloatComplex *U, const magma_int_t *ldu,
                         cuFloatComplex *work, const magma_int_t *lwork,
                         float *rwork,
                         float *resid );
#else
void   lapackf77_cbdt01( const magma_int_t *m, const magma_int_t *n, const magma_int_t *kd,
                         cuFloatComplex *A, const magma_int_t *lda,
                         cuFloatComplex *Q, const magma_int_t *ldq,
                         float *d, float *e,
                         cuFloatComplex *Pt, const magma_int_t *ldpt,
                         cuFloatComplex *work,
                         float *resid );

void   lapackf77_cget22( const char *transa, const char *transe, const char *transw, const magma_int_t *n,
                         cuFloatComplex *A, const magma_int_t *lda,
                         cuFloatComplex *E, const magma_int_t *lde,
                         cuFloatComplex *wr,
                         cuFloatComplex *wi,
                         float *work,
                         float *result );

void   lapackf77_chet21( magma_int_t *itype, const char *uplo, const magma_int_t *n, const magma_int_t *kband,
                         cuFloatComplex *A, const magma_int_t *lda,
                         float *d, float *e,
                         cuFloatComplex *U, const magma_int_t *ldu,
                         cuFloatComplex *V, const magma_int_t *ldv,
                         cuFloatComplex *tau,
                         cuFloatComplex *work,
                         float *result );

void   lapackf77_chst01( const magma_int_t *n, const magma_int_t *ilo, const magma_int_t *ihi,
                         cuFloatComplex *A, const magma_int_t *lda,
                         cuFloatComplex *H, const magma_int_t *ldh,
                         cuFloatComplex *Q, const magma_int_t *ldq,
                         cuFloatComplex *work, const magma_int_t *lwork,
                         float *result );

void   lapackf77_cstt21( const magma_int_t *n, const magma_int_t *kband,
                         float *AD,
                         float *AE,
                         float *SD,
                         float *SE,
                         cuFloatComplex *U, const magma_int_t *ldu,
                         cuFloatComplex *work,
                         float *result );

void   lapackf77_cunt01( const char *rowcol, const magma_int_t *m, const magma_int_t *n,
                         cuFloatComplex *U, const magma_int_t *ldu,
                         cuFloatComplex *work, const magma_int_t *lwork,
                         float *resid );
#endif

void   lapackf77_clarfy( const char *uplo, const magma_int_t *n,
                         cuFloatComplex *V, const magma_int_t *incv,
                         cuFloatComplex *tau,
                         cuFloatComplex *C, const magma_int_t *ldc,
                         cuFloatComplex *work );

void   lapackf77_clarfx( const char *side, const magma_int_t *m, const magma_int_t *n,
                         cuFloatComplex *V,
                         cuFloatComplex *tau,
                         cuFloatComplex *C, const magma_int_t *ldc,
                         cuFloatComplex *work );

float lapackf77_cqpt01( const magma_int_t *m, const magma_int_t *n, const magma_int_t *k,
                         cuFloatComplex *A,
                         cuFloatComplex *Af, const magma_int_t *lda,
                         cuFloatComplex *tau, magma_int_t *jpvt,
                         cuFloatComplex *work, const magma_int_t *lwork );

void   lapackf77_cqrt02( const magma_int_t *m, const magma_int_t *n, const magma_int_t *k,
                         cuFloatComplex *A,
                         cuFloatComplex *AF,
                         cuFloatComplex *Q,
                         cuFloatComplex *R, const magma_int_t *lda,
                         cuFloatComplex *tau,
                         cuFloatComplex *work, const magma_int_t *lwork,
                         float *rwork,
                         float *result );

#ifdef __cplusplus
}
#endif

#undef DWORKFORZ
#undef DWORKFORZ_AND_LD
#undef WSPLIT
#undef PRECISION_c

#endif /* MAGMA_CLAPACK_H */
