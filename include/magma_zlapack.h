/*
 *   -- MAGMA (version 1.3.0) --
 *      Univ. of Tennessee, Knoxville
 *      Univ. of California, Berkeley
 *      Univ. of Colorado, Denver
 *      November 2012
 *
 * @precisions normal z -> s d c
 */

#ifndef MAGMA_ZLAPACK_H
#define MAGMA_ZLAPACK_H

#define PRECISION_z

#ifdef __cplusplus
extern "C" {
#endif

/* ////////////////////////////////////////////////////////////////////////////
   -- BLAS and LAPACK functions (alphabetical order)
*/
#define blasf77_izamax     FORTRAN_NAME( izamax, IZAMAX )
#define blasf77_zaxpy      FORTRAN_NAME( zaxpy,  ZAXPY  )
#define blasf77_zcopy      FORTRAN_NAME( zcopy,  ZCOPY  )
#define blasf77_zgemm      FORTRAN_NAME( zgemm,  ZGEMM  )
#define blasf77_zgemv      FORTRAN_NAME( zgemv,  ZGEMV  )
#define blasf77_zgerc      FORTRAN_NAME( zgerc,  ZGERC  )
#define blasf77_zgeru      FORTRAN_NAME( zgeru,  ZGERU  )
#define blasf77_zhemm      FORTRAN_NAME( zhemm,  ZHEMM  )
#define blasf77_zhemv      FORTRAN_NAME( zhemv,  ZHEMV  )
#define blasf77_zher2      FORTRAN_NAME( zher2,  ZHER2  )
#define blasf77_zher2k     FORTRAN_NAME( zher2k, ZHER2K )
#define blasf77_zherk      FORTRAN_NAME( zherk,  ZHERK  )
#define blasf77_zscal      FORTRAN_NAME( zscal,  ZSCAL  )
#define blasf77_zdscal     FORTRAN_NAME( zdscal, ZDSCAL )
#define blasf77_zswap      FORTRAN_NAME( zswap,  ZSWAP  )
#define blasf77_zsymm      FORTRAN_NAME( zsymm,  ZSYMM  )
#define blasf77_zsyr2k     FORTRAN_NAME( zsyr2k, ZSYR2K )
#define blasf77_zsyrk      FORTRAN_NAME( zsyrk,  ZSYRK  )
#define blasf77_ztrmm      FORTRAN_NAME( ztrmm,  ZTRMM  )
#define blasf77_ztrmv      FORTRAN_NAME( ztrmv,  ZTRMV  )
#define blasf77_ztrsm      FORTRAN_NAME( ztrsm,  ZTRSM  )
#define blasf77_ztrsv      FORTRAN_NAME( ztrsv,  ZTRSV  )

#define lapackf77_dlaed4   FORTRAN_NAME( dlaed4, DLAED4 )
#define lapackf77_dlamc3   FORTRAN_NAME( dlamc3, DLAMC3 )
#define lapackf77_dlamrg   FORTRAN_NAME( dlamrg, DLAMRG )
#define lapackf77_dstebz   FORTRAN_NAME( dstebz, DSTEBZ )

#define lapackf77_zbdsqr   FORTRAN_NAME( zbdsqr, ZBDSQR )
#define lapackf77_zgebak   FORTRAN_NAME( zgebak, ZGEBAK )
#define lapackf77_zgebal   FORTRAN_NAME( zgebal, ZGEBAL )
#define lapackf77_zgebd2   FORTRAN_NAME( zgebd2, ZGEBD2 )
#define lapackf77_zgebrd   FORTRAN_NAME( zgebrd, ZGEBRD )
#define lapackf77_zgeev    FORTRAN_NAME( zgeev,  ZGEEV  )
#define lapackf77_zgehd2   FORTRAN_NAME( zgehd2, ZGEHD2 )
#define lapackf77_zgehrd   FORTRAN_NAME( zgehrd, ZGEHRD )
#define lapackf77_zgelqf   FORTRAN_NAME( zgelqf, ZGELQF )
#define lapackf77_zgels    FORTRAN_NAME( zgels,  ZGELS  )
#define lapackf77_zgeqlf   FORTRAN_NAME( zgeqlf, ZGEQLF )
#define lapackf77_zgeqp3   FORTRAN_NAME( zgeqp3, ZGEQP3 )
#define lapackf77_zgeqrf   FORTRAN_NAME( zgeqrf, ZGEQRF )
#define lapackf77_zgesvd   FORTRAN_NAME( zgesvd, ZGESVD )
#define lapackf77_zgetrf   FORTRAN_NAME( zgetrf, ZGETRF )
#define lapackf77_zgetri   FORTRAN_NAME( zgetri, ZGETRI )
#define lapackf77_zgetrs   FORTRAN_NAME( zgetrs, ZGETRS )
#define lapackf77_zhbtrd   FORTRAN_NAME( zhbtrd, ZHBTRD )
#define lapackf77_zheev    FORTRAN_NAME( zheev,  ZHEEV  )
#define lapackf77_zheevd   FORTRAN_NAME( zheevd, ZHEEVD )
#define lapackf77_zhegs2   FORTRAN_NAME( zhegs2, ZHEGS2 )
#define lapackf77_zhegvd   FORTRAN_NAME( zhegvd, ZHEGVD )
#define lapackf77_zhetd2   FORTRAN_NAME( zhetd2, ZHETD2 )
#define lapackf77_zhetrd   FORTRAN_NAME( zhetrd, ZHETRD )
#define lapackf77_zhetrf   FORTRAN_NAME( zhetrf, ZHETRF )
#define lapackf77_zhseqr   FORTRAN_NAME( zhseqr, ZHSEQR )
#define lapackf77_zlabrd   FORTRAN_NAME( zlabrd, ZLABRD )
#define lapackf77_zladiv   FORTRAN_NAME( zladiv, ZLADIV )
#define lapackf77_zlacgv   FORTRAN_NAME( zlacgv, ZLACGV )
#define lapackf77_zlacpy   FORTRAN_NAME( zlacpy, ZLACPY )
#define lapackf77_zlahef   FORTRAN_NAME( zlahef, ZLAHEF )
#define lapackf77_zlange   FORTRAN_NAME( zlange, ZLANGE )
#define lapackf77_zlanhe   FORTRAN_NAME( zlanhe, ZLANHE )
#define lapackf77_zlanht   FORTRAN_NAME( zlanht, ZLANHT )
#define lapackf77_zlansy   FORTRAN_NAME( zlansy, ZLANSY )
#define lapackf77_dlapy3   FORTRAN_NAME( dlapy3, DLAPY3 )
#define lapackf77_zlaqp2   FORTRAN_NAME( zlaqp2, ZLAQP2 )
#define lapackf77_zlarf    FORTRAN_NAME( zlarf,  ZLARF  )
#define lapackf77_zlarfb   FORTRAN_NAME( zlarfb, ZLARFB )
#define lapackf77_zlarfg   FORTRAN_NAME( zlarfg, ZLARFG )
#define lapackf77_zlarft   FORTRAN_NAME( zlarft, ZLARFT )
#define lapackf77_zlarnv   FORTRAN_NAME( zlarnv, ZLARNV )
#define lapackf77_zlartg   FORTRAN_NAME( zlartg, ZLARTG )
#define lapackf77_zlascl   FORTRAN_NAME( zlascl, ZLASCL )
#define lapackf77_zlaset   FORTRAN_NAME( zlaset, ZLASET )
#define lapackf77_zlaswp   FORTRAN_NAME( zlaswp, ZLASWP )
#define lapackf77_zlatrd   FORTRAN_NAME( zlatrd, ZLATRD )
#define lapackf77_zlauum   FORTRAN_NAME( zlauum, ZLAUUM )
#define lapackf77_zlavhe   FORTRAN_NAME( zlavhe, ZLAVHE )
#define lapackf77_zpotrf   FORTRAN_NAME( zpotrf, ZPOTRF )
#define lapackf77_zpotri   FORTRAN_NAME( zpotri, ZPOTRI )
#define lapackf77_zpotrs   FORTRAN_NAME( zpotrs, ZPOTRS )
#define lapackf77_zstedc   FORTRAN_NAME( zstedc, ZSTEDC )
#define lapackf77_zstein   FORTRAN_NAME( zstein, ZSTEIN )
#define lapackf77_zstemr   FORTRAN_NAME( zstemr, ZSTEMR )
#define lapackf77_zsteqr   FORTRAN_NAME( zsteqr, ZSTEQR )
#define lapackf77_zsymv    FORTRAN_NAME( zsymv,  ZSYMV  )
#define lapackf77_ztrevc   FORTRAN_NAME( ztrevc, ZTREVC )
#define lapackf77_ztrtri   FORTRAN_NAME( ztrtri, ZTRTRI )
#define lapackf77_zung2r   FORTRAN_NAME( zung2r, ZUNG2R )
#define lapackf77_zungbr   FORTRAN_NAME( zungbr, ZUNGBR )
#define lapackf77_zunghr   FORTRAN_NAME( zunghr, ZUNGHR )
#define lapackf77_zunglq   FORTRAN_NAME( zunglq, ZUNGLQ )
#define lapackf77_zungql   FORTRAN_NAME( zungql, ZUNGQL )
#define lapackf77_zungqr   FORTRAN_NAME( zungqr, ZUNGQR )
#define lapackf77_zungtr   FORTRAN_NAME( zungtr, ZUNGTR )
#define lapackf77_zunm2r   FORTRAN_NAME( zunm2r, ZUNM2R )
#define lapackf77_zunmbr   FORTRAN_NAME( zunmbr, ZUNMBR )
#define lapackf77_zunmlq   FORTRAN_NAME( zunmlq, ZUNMLQ )
#define lapackf77_zunmql   FORTRAN_NAME( zunmql, ZUNMQL )
#define lapackf77_zunmqr   FORTRAN_NAME( zunmqr, ZUNMQR )
#define lapackf77_zunmtr   FORTRAN_NAME( zunmtr, ZUNMTR )

/* testing functions */
#define lapackf77_zbdt01   FORTRAN_NAME( zbdt01, ZBDT01 )
#define lapackf77_zget22   FORTRAN_NAME( zget22, ZGET22 )
#define lapackf77_zqpt01   FORTRAN_NAME( zqpt01, ZQPT01 )
#define lapackf77_zhet21   FORTRAN_NAME( zhet21, ZHET21 )
#define lapackf77_zhst01   FORTRAN_NAME( zhst01, ZHST01 )
#define lapackf77_zqrt02   FORTRAN_NAME( zqrt02, ZQRT02 )
#define lapackf77_zunt01   FORTRAN_NAME( zunt01, ZUNT01 )
#define lapackf77_zlarfy   FORTRAN_NAME( zlarfy, ZLARFY )
#define lapackf77_zlarfx   FORTRAN_NAME( zlarfx, ZLARFX )
#define lapackf77_zstt21   FORTRAN_NAME( zstt21, ZSTT21 )

// macros to handle differences in arguments between complex and real versions of routines.
#if defined(PRECISION_z) || defined(PRECISION_c)
#define DWORKFORZ        double *rwork,
#define DWORKFORZ_AND_LD double *rwork, const magma_int_t *ldrwork,
#define WSPLIT           cuDoubleComplex *w
#else
#define DWORKFORZ
#define DWORKFORZ_AND_LD
#define WSPLIT           double *wr, double *wi
#endif

/*
 * BLAS functions (alphabetical order)
 */
magma_int_t blasf77_izamax(
                     const magma_int_t *n,
                     const cuDoubleComplex *x, const magma_int_t *incx);

void blasf77_zaxpy(  const magma_int_t *n,
                     const cuDoubleComplex *alpha,
                     const cuDoubleComplex *x, const magma_int_t *incx,
                           cuDoubleComplex *y, const magma_int_t *incy );

void blasf77_zcopy(  const magma_int_t *n,
                     const cuDoubleComplex *x, const magma_int_t *incx,
                           cuDoubleComplex *y, const magma_int_t *incy );

void blasf77_zgemm(  const char *transa, const char *transb,
                     const magma_int_t *m, const magma_int_t *n, const magma_int_t *k,
                     const cuDoubleComplex *alpha,
                     const cuDoubleComplex *A, const magma_int_t *lda,
                     const cuDoubleComplex *B, const magma_int_t *ldb,
                     const cuDoubleComplex *beta,
                           cuDoubleComplex *C, const magma_int_t *ldc );

void blasf77_zgemv(  const char *transa,
                     const magma_int_t *m, const magma_int_t *n,
                     const cuDoubleComplex *alpha,
                     const cuDoubleComplex *A, const magma_int_t *lda,
                     const cuDoubleComplex *x, const magma_int_t *incx,
                     const cuDoubleComplex *beta,
                           cuDoubleComplex *y, const magma_int_t *incy );

void blasf77_zgerc(  const magma_int_t *m, const magma_int_t *n,
                     const cuDoubleComplex *alpha,
                     const cuDoubleComplex *x, const magma_int_t *incx,
                     const cuDoubleComplex *y, const magma_int_t *incy,
                           cuDoubleComplex *A, const magma_int_t *lda );

#if defined(PRECISION_z) || defined(PRECISION_c)
void blasf77_zgeru(  const magma_int_t *m, const magma_int_t *n,
                     const cuDoubleComplex *alpha,
                     const cuDoubleComplex *x, const magma_int_t *incx,
                     const cuDoubleComplex *y, const magma_int_t *incy,
                           cuDoubleComplex *A, const magma_int_t *lda );
#endif

void blasf77_zhemm(  const char *side, const char *uplo,
                     const magma_int_t *m, const magma_int_t *n,
                     const cuDoubleComplex *alpha,
                     const cuDoubleComplex *A, const magma_int_t *lda,
                     const cuDoubleComplex *B, const magma_int_t *ldb,
                     const cuDoubleComplex *beta,
                           cuDoubleComplex *C, const magma_int_t *ldc );

void blasf77_zhemv(  const char *uplo,
                     const magma_int_t *n,
                     const cuDoubleComplex *alpha,
                     const cuDoubleComplex *A, const magma_int_t *lda,
                     const cuDoubleComplex *x, const magma_int_t *incx,
                     const cuDoubleComplex *beta,
                           cuDoubleComplex *y, const magma_int_t *incy );

void blasf77_zher2(  const char *uplo,
                     const magma_int_t *n,
                     const cuDoubleComplex *alpha,
                     const cuDoubleComplex *x, const magma_int_t *incx,
                     const cuDoubleComplex *y, const magma_int_t *incy,
                           cuDoubleComplex *A, const magma_int_t *lda );

void blasf77_zher2k(  const char *uplo, const char *trans,
                     const magma_int_t *n, const magma_int_t *k,
                     const cuDoubleComplex *alpha,
                     const cuDoubleComplex *A, const magma_int_t *lda,
                     const cuDoubleComplex *B, const magma_int_t *ldb,
                     const double *beta,
                           cuDoubleComplex *C, const magma_int_t *ldc );

void blasf77_zherk(  const char *uplo, const char *trans,
                     const magma_int_t *n, const magma_int_t *k,
                     const double *alpha,
                     const cuDoubleComplex *A, const magma_int_t *lda,
                     const double *beta,
                           cuDoubleComplex *C, const magma_int_t *ldc );

void blasf77_zscal(  const magma_int_t *n,
                     const cuDoubleComplex *alpha,
                           cuDoubleComplex *x, const magma_int_t *incx );

#if defined(PRECISION_z) || defined(PRECISION_c)
void blasf77_zdscal( const magma_int_t *n,
                     const double *alpha,
                           cuDoubleComplex *x, const magma_int_t *incx );
#endif

void blasf77_zswap(  const magma_int_t *n,
                     cuDoubleComplex *x, const magma_int_t *incx,
                     cuDoubleComplex *y, const magma_int_t *incy );

void blasf77_zsymm(  const char *side, const char *uplo,
                     const magma_int_t *m, const magma_int_t *n,
                     const cuDoubleComplex *alpha,
                     const cuDoubleComplex *A, const magma_int_t *lda,
                     const cuDoubleComplex *B, const magma_int_t *ldb,
                     const cuDoubleComplex *beta,
                           cuDoubleComplex *C, const magma_int_t *ldc );

void blasf77_zsyr2k( const char *uplo, const char *trans,
                     const magma_int_t *n, const magma_int_t *k,
                     const cuDoubleComplex *alpha,
                     const cuDoubleComplex *A, const magma_int_t *lda,
                     const cuDoubleComplex *B, const magma_int_t *ldb,
                     const cuDoubleComplex *beta,
                           cuDoubleComplex *C, const magma_int_t *ldc );

void blasf77_zsyrk(  const char *uplo, const char *trans,
                     const magma_int_t *n, const magma_int_t *k,
                     const cuDoubleComplex *alpha,
                     const cuDoubleComplex *A, const magma_int_t *lda,
                     const cuDoubleComplex *beta,
                           cuDoubleComplex *C, const magma_int_t *ldc );

void blasf77_ztrmm(  const char *side, const char *uplo, const char *transa, const char *diag,
                     const magma_int_t *m, const magma_int_t *n,
                     const cuDoubleComplex *alpha,
                     const cuDoubleComplex *A, const magma_int_t *lda,
                           cuDoubleComplex *B, const magma_int_t *ldb );

void blasf77_ztrmv(  const char *uplo, const char *transa, const char *diag,
                     const magma_int_t *n,
                     const cuDoubleComplex *A, const magma_int_t *lda,
                           cuDoubleComplex *x, const magma_int_t *incx );

void blasf77_ztrsm(  const char *side, const char *uplo, const char *transa, const char *diag,
                     const magma_int_t *m, const magma_int_t *n,
                     const cuDoubleComplex *alpha,
                     const cuDoubleComplex *A, const magma_int_t *lda,
                           cuDoubleComplex *B, const magma_int_t *ldb );

void blasf77_ztrsv(  const char *uplo, const char *transa, const char *diag,
                     const magma_int_t *n,
                     const cuDoubleComplex *A, const magma_int_t *lda,
                           cuDoubleComplex *x, const magma_int_t *incx );

/*
 * LAPACK functions (alphabetical order)
 */
void   lapackf77_zbdsqr( const char *uplo,
                         const magma_int_t *n, const magma_int_t *ncvt, const magma_int_t *nru,  const magma_int_t *ncc,
                         double *d, double *e,
                         cuDoubleComplex *Vt, const magma_int_t *ldvt,
                         cuDoubleComplex *U, const magma_int_t *ldu,
                         cuDoubleComplex *C, const magma_int_t *ldc,
                         double *work,
                         magma_int_t *info );

void   lapackf77_zgebak( const char *job, const char *side,
                         const magma_int_t *n,
                         const magma_int_t *ilo, const magma_int_t *ihi,
                         const double *scale, const magma_int_t *m,
                         cuDoubleComplex *V, const magma_int_t *ldv,
                         magma_int_t *info );

void   lapackf77_zgebal( const char *job,
                         const magma_int_t *n,
                         cuDoubleComplex *A, const magma_int_t *lda,
                         magma_int_t *ilo, magma_int_t *ihi,
                         double *scale,
                         magma_int_t *info );

void   lapackf77_zgebd2( const magma_int_t *m, const magma_int_t *n,
                         cuDoubleComplex *A, const magma_int_t *lda,
                         double *d, double *e,
                         cuDoubleComplex *tauq,
                         cuDoubleComplex *taup,
                         cuDoubleComplex *work,
                         magma_int_t *info );

void   lapackf77_zgebrd( const magma_int_t *m, const magma_int_t *n,
                         cuDoubleComplex *A, const magma_int_t *lda,
                         double *d, double *e,
                         cuDoubleComplex *tauq,
                         cuDoubleComplex *taup,
                         cuDoubleComplex *work, const magma_int_t *lwork,
                         magma_int_t *info );

void   lapackf77_zgeev(  const char *jobvl, const char *jobvr,
                         const magma_int_t *n,
                         cuDoubleComplex *A,    const magma_int_t *lda,
                         WSPLIT,
                         cuDoubleComplex *Vl,   const magma_int_t *ldvl,
                         cuDoubleComplex *Vr,   const magma_int_t *ldvr,
                         cuDoubleComplex *work, const magma_int_t *lwork,
                         DWORKFORZ
                         magma_int_t *info );

void   lapackf77_zgehd2( const magma_int_t *n,
                         const magma_int_t *ilo, const magma_int_t *ihi,
                         cuDoubleComplex *A, const magma_int_t *lda,
                         cuDoubleComplex *tau,
                         cuDoubleComplex *work,
                         magma_int_t *info );

void   lapackf77_zgehrd( const magma_int_t *n,
                         const magma_int_t *ilo, const magma_int_t *ihi,
                         cuDoubleComplex *A, const magma_int_t *lda,
                         cuDoubleComplex *tau,
                         cuDoubleComplex *work, const magma_int_t *lwork,
                         magma_int_t *info );

void   lapackf77_zgelqf( const magma_int_t *m, const magma_int_t *n,
                         cuDoubleComplex *A, const magma_int_t *lda,
                         cuDoubleComplex *tau,
                         cuDoubleComplex *work, const magma_int_t *lwork,
                         magma_int_t *info );

void   lapackf77_zgels(  const char *trans,
                         const magma_int_t *m, const magma_int_t *n, const magma_int_t *nrhs,
                         cuDoubleComplex *A, const magma_int_t *lda,
                         cuDoubleComplex *B, const magma_int_t *ldb,
                         cuDoubleComplex *work, const magma_int_t *lwork,
                         magma_int_t *info );

void   lapackf77_zgeqlf( const magma_int_t *m, const magma_int_t *n,
                         cuDoubleComplex *A, const magma_int_t *lda,
                         cuDoubleComplex *tau,
                         cuDoubleComplex *work, const magma_int_t *lwork,
                         magma_int_t *info );

void   lapackf77_zgeqp3( const magma_int_t *m, const magma_int_t *n,
                         cuDoubleComplex *A, const magma_int_t *lda,
                         magma_int_t *jpvt,
                         cuDoubleComplex *tau,
                         cuDoubleComplex *work, const magma_int_t *lwork,
                         DWORKFORZ
                         magma_int_t *info );

void   lapackf77_zgeqrf( const magma_int_t *m, const magma_int_t *n,
                         cuDoubleComplex *A, const magma_int_t *lda,
                         cuDoubleComplex *tau,
                         cuDoubleComplex *work, const magma_int_t *lwork,
                         magma_int_t *info );

void   lapackf77_zgesvd( const char *jobu, const char *jobvt,
                         const magma_int_t *m, const magma_int_t *n,
                         cuDoubleComplex *A, const magma_int_t *lda,
                         double *s,
                         cuDoubleComplex *U,  const magma_int_t *ldu,
                         cuDoubleComplex *Vt, const magma_int_t *ldvt,
                         cuDoubleComplex *work, const magma_int_t *lwork,
                         DWORKFORZ
                         magma_int_t *info );

void   lapackf77_zgetrf( const magma_int_t *m, const magma_int_t *n,
                         cuDoubleComplex *A, const magma_int_t *lda,
                         magma_int_t *ipiv,
                         magma_int_t *info );

void   lapackf77_zgetri( const magma_int_t *n,
                         cuDoubleComplex *A, const magma_int_t *lda,
                         const magma_int_t *ipiv,
                         cuDoubleComplex *work, const magma_int_t *lwork,
                         magma_int_t *info );

void   lapackf77_zgetrs( const char* trans,
                         const magma_int_t *n, const magma_int_t *nrhs,
                         const cuDoubleComplex *A, const magma_int_t *lda,
                         const magma_int_t *ipiv,
                         cuDoubleComplex *B, const magma_int_t *ldb,
                         magma_int_t *info );

void   lapackf77_zhbtrd( const char *vect, const char *uplo,
                         const magma_int_t *n, const magma_int_t *kd,
                         cuDoubleComplex *Ab, const magma_int_t *ldab,
                         double *d, double *e,
                         cuDoubleComplex *Q, const magma_int_t *ldq,
                         cuDoubleComplex *work,
                         magma_int_t *info );

void   lapackf77_zheev(  const char *jobz, const char *uplo,
                         const magma_int_t *n,
                         cuDoubleComplex *A, const magma_int_t *lda,
                         double *w,
                         cuDoubleComplex *work, const magma_int_t *lwork,
                         DWORKFORZ
                         magma_int_t *info );

void   lapackf77_zheevd( const char *jobz, const char *uplo,
                         const magma_int_t *n,
                         cuDoubleComplex *A, const magma_int_t *lda,
                         double *w,
                         cuDoubleComplex *work, const magma_int_t *lwork,
                         DWORKFORZ_AND_LD
                         magma_int_t *iwork, const magma_int_t *liwork,
                         magma_int_t *info );

void   lapackf77_zhegs2( const magma_int_t *itype, const char *uplo,
                         const magma_int_t *n,
                         cuDoubleComplex *A, const magma_int_t *lda,
                         cuDoubleComplex *B, const magma_int_t *ldb,
                         magma_int_t *info );

void   lapackf77_zhegvd( const magma_int_t *itype, const char *jobz, const char *uplo,
                         const magma_int_t *n,
                         cuDoubleComplex *A, const magma_int_t *lda,
                         cuDoubleComplex *B, const magma_int_t *ldb,
                         double *w,
                         cuDoubleComplex *work, const magma_int_t *lwork,
                         DWORKFORZ_AND_LD
                         magma_int_t *iwork, const magma_int_t *liwork,
                         magma_int_t *info );

void   lapackf77_zhetd2( const char *uplo,
                         const magma_int_t *n,
                         cuDoubleComplex *A, const magma_int_t *lda,
                         double *d, double *e,
                         cuDoubleComplex *tau,
                         magma_int_t *info );

void   lapackf77_zhetrd( const char *uplo,
                         const magma_int_t *n,
                         cuDoubleComplex *A, const magma_int_t *lda,
                         double *d, double *e,
                         cuDoubleComplex *tau,
                         cuDoubleComplex *work, const magma_int_t *lwork,
                         magma_int_t *info );

void   lapackf77_zhetrf( const char *uplo,
                         const magma_int_t *n,
                         cuDoubleComplex *A, const magma_int_t *lda,
                         magma_int_t *ipiv,
                         cuDoubleComplex *work, const magma_int_t *lwork,
                         magma_int_t *info );

void   lapackf77_zhseqr( const char *job, const char *compz,
                         const magma_int_t *n,
                         const magma_int_t *ilo, const magma_int_t *ihi,
                         cuDoubleComplex *H, const magma_int_t *ldh,
                         WSPLIT,
                         cuDoubleComplex *Z, const magma_int_t *ldz,
                         cuDoubleComplex *work, const magma_int_t *lwork,
                         magma_int_t *info );

void   lapackf77_zlabrd( const magma_int_t *m, const magma_int_t *n, const magma_int_t *nb,
                         cuDoubleComplex *A, const magma_int_t *lda,
                         double *d, double *e,
                         cuDoubleComplex *tauq,
                         cuDoubleComplex *taup,
                         cuDoubleComplex *X, const magma_int_t *ldx,
                         cuDoubleComplex *Y, const magma_int_t *ldy );

void   lapackf77_zladiv( cuDoubleComplex *ret_val, cuDoubleComplex *x, 
                         cuDoubleComplex *y );

void   lapackf77_zlacgv( const magma_int_t *n,
                         cuDoubleComplex *x, const magma_int_t *incx );

void   lapackf77_zlacpy( const char *uplo,
                         const magma_int_t *m, const magma_int_t *n,
                         const cuDoubleComplex *A, const magma_int_t *lda,
                         cuDoubleComplex *B, const magma_int_t *ldb );

void   lapackf77_zlahef( const char *uplo,
                         const magma_int_t *n, const magma_int_t *kn,
                         magma_int_t *kb,
                         cuDoubleComplex *A, const magma_int_t lda,
                         magma_int_t *ipiv,
                         cuDoubleComplex *work, const magma_int_t *ldwork,
                         magma_int_t *info );

double lapackf77_zlange( const char *norm,
                         const magma_int_t *m, const magma_int_t *n,
                         const cuDoubleComplex *A, const magma_int_t *lda,
                         double *work );

double lapackf77_zlanhe( const char *norm, const char *uplo,
                         const magma_int_t *n,
                         const cuDoubleComplex *A, const magma_int_t *lda,
                         double * work );

double lapackf77_zlanht( const char* norm, const magma_int_t* n,
                         const double* d, const cuDoubleComplex* e );

double lapackf77_zlansy( const char *norm, const char *uplo,
                         const magma_int_t *n,
                         const cuDoubleComplex *A, const magma_int_t *lda,
                         double * work );

void lapackf77_zlaqp2 (  magma_int_t *m, magma_int_t *n, magma_int_t *offset,
                         cuDoubleComplex *a, magma_int_t *lda, magma_int_t *jpvt, 
                         cuDoubleComplex *tau,
                         double *vn1, double *vn2, cuDoubleComplex *work);

void lapackf77_zlarf  (  char *, magma_int_t *, magma_int_t *,
                         cuDoubleComplex *, magma_int_t *, cuDoubleComplex *, cuDoubleComplex *,
                         magma_int_t *, cuDoubleComplex *);

void   lapackf77_zlarfb( const char *side, const char *trans, const char *direct, const char *storev,
                         const magma_int_t *m, const magma_int_t *n, const magma_int_t *k,
                         const cuDoubleComplex *V, const magma_int_t *ldv,
                         const cuDoubleComplex *T, const magma_int_t *ldt,
                         cuDoubleComplex *C, const magma_int_t *ldc,
                         cuDoubleComplex *work, const magma_int_t *ldwork );

void   lapackf77_zlarfg( const magma_int_t *n,
                         cuDoubleComplex *alpha,
                         cuDoubleComplex *x, const magma_int_t *incx,
                         cuDoubleComplex *tau );

void   lapackf77_zlarft( const char *direct, const char *storev,
                         const magma_int_t *n, const magma_int_t *k,
                         cuDoubleComplex *V, const magma_int_t *ldv,
                         const cuDoubleComplex *tau,
                         cuDoubleComplex *T, const magma_int_t *ldt );

void   lapackf77_zlarnv( const magma_int_t *idist, magma_int_t *iseed, const magma_int_t *n,
                         cuDoubleComplex *x );

void   lapackf77_zlartg( cuDoubleComplex *F,
                         cuDoubleComplex *G,
                         double *cs,
                         cuDoubleComplex *SN,
                         cuDoubleComplex *R );

void   lapackf77_zlascl( const char *type,
                         const magma_int_t *kl, const magma_int_t *ku,
                         double *cfrom,
                         double *cto,
                         const magma_int_t *m, const magma_int_t *n,
                         cuDoubleComplex *A, const magma_int_t *lda,
                         magma_int_t *info );

void   lapackf77_zlaset( const char *uplo,
                         const magma_int_t *m, const magma_int_t *n,
                         const cuDoubleComplex *alpha,
                         const cuDoubleComplex *beta,
                         cuDoubleComplex *A, const magma_int_t *lda );

void   lapackf77_zlaswp( const magma_int_t *n,
                         cuDoubleComplex *A, const magma_int_t *lda,
                         const magma_int_t *k1, const magma_int_t *k2,
                         magma_int_t *ipiv,
                         const magma_int_t *incx );

void   lapackf77_zlatrd( const char *uplo,
                         const magma_int_t *n, const magma_int_t *nb,
                         cuDoubleComplex *A, const magma_int_t *lda,
                         double *e,
                         cuDoubleComplex *tau,
                         cuDoubleComplex *work, const magma_int_t *ldwork );

void   lapackf77_zlauum( const char *uplo,
                         const magma_int_t *n,
                         cuDoubleComplex *A, const magma_int_t *lda,
                         magma_int_t *info );

void   lapackf77_zlavhe( const char *uplo, const char *trans, const char *diag,
                         magma_int_t *n, magma_int_t *nrhs,
                         cuDoubleComplex *A, magma_int_t *lda,
                         magma_int_t *ipiv,
                         cuDoubleComplex *B, magma_int_t *ldb,
                         magma_int_t *info );

void   lapackf77_zpotrf( const char *uplo,
                         const magma_int_t *n,
                         cuDoubleComplex *A, const magma_int_t *lda,
                         magma_int_t *info );

void   lapackf77_zpotri( const char *uplo,
                         const magma_int_t *n,
                         cuDoubleComplex *A, const magma_int_t *lda,
                         magma_int_t *info );

void   lapackf77_zpotrs( const char *uplo,
                         const magma_int_t *n, const magma_int_t *nrhs,
                         const cuDoubleComplex *A, const magma_int_t *lda,
                         cuDoubleComplex *B, const magma_int_t *ldb,
                         magma_int_t *info );

void   lapackf77_zstedc( const char *compz,
                         const magma_int_t *n,
                         double *d, double *e,
                         cuDoubleComplex *Z, const magma_int_t *ldz,
                         cuDoubleComplex *work, const magma_int_t *lwork,
                         DWORKFORZ_AND_LD
                         magma_int_t *iwork, const magma_int_t *liwork,
                         magma_int_t *info );

void   lapackf77_zstein( const magma_int_t *n,
                         const double *d, const double *e,
                         const magma_int_t *m,
                         const double *w,
                         const magma_int_t *iblock,
                         const magma_int_t *isplit,
                         cuDoubleComplex *Z, const magma_int_t *ldz,
                         double *work, magma_int_t *iwork, magma_int_t *ifailv,
                         magma_int_t *info );

void   lapackf77_zstemr( const char *jobz, const char *range,
                         const magma_int_t *n,
                         double *d, double *e,
                         const double *vl, const double *vu,
                         const magma_int_t *il, const magma_int_t *iu,
                         magma_int_t *m,
                         double *w,
                         cuDoubleComplex *Z, const magma_int_t *ldz,
                         const magma_int_t *nzc, magma_int_t *isuppz, magma_int_t *tryrac,
                         double *work, const magma_int_t *lwork,
                         magma_int_t *iwork, const magma_int_t *liwork,
                         magma_int_t *info );

void   lapackf77_zsteqr( const char *compz,
                         const magma_int_t *n,
                         double *d, double *e,
                         cuDoubleComplex *Z, const magma_int_t *ldz,
                         double *work,
                         magma_int_t *info );

#if defined(PRECISION_z) || defined(PRECISION_c)
void   lapackf77_zsymv(  const char *uplo,
                         const magma_int_t *n,
                         const cuDoubleComplex *alpha,
                         const cuDoubleComplex *A, const magma_int_t *lda,
                         const cuDoubleComplex *x, const magma_int_t *incx,
                         const cuDoubleComplex *beta,
                               cuDoubleComplex *y, const magma_int_t *incy );
#endif

void   lapackf77_ztrevc( const char *side, const char *howmny,
                         magma_int_t *select, const magma_int_t *n,
                         cuDoubleComplex *T,  const magma_int_t *ldt,
                         cuDoubleComplex *Vl, const magma_int_t *ldvl,
                         cuDoubleComplex *Vr, const magma_int_t *ldvr,
                         const magma_int_t *mm, magma_int_t *m,
                         cuDoubleComplex *work,
                         DWORKFORZ
                         magma_int_t *info );

void   lapackf77_ztrtri( const char *uplo, const char *diag,
                         const magma_int_t *n,
                         cuDoubleComplex *A, const magma_int_t *lda,
                         magma_int_t *info );

void   lapackf77_zung2r( const magma_int_t *m, const magma_int_t *n, const magma_int_t *k,
                         cuDoubleComplex *A, const magma_int_t *lda,
                         const cuDoubleComplex *tau,
                         cuDoubleComplex *work,
                         magma_int_t *info );

void   lapackf77_zungbr( const char *vect,
                         const magma_int_t *m, const magma_int_t *n, const magma_int_t *k,
                         cuDoubleComplex *A, const magma_int_t *lda,
                         const cuDoubleComplex *tau,
                         cuDoubleComplex *work, const magma_int_t *lwork,
                         magma_int_t *info );

void   lapackf77_zunghr( const magma_int_t *n,
                         const magma_int_t *ilo, const magma_int_t *ihi,
                         cuDoubleComplex *A, const magma_int_t *lda,
                         const cuDoubleComplex *tau,
                         cuDoubleComplex *work, const magma_int_t *lwork,
                         magma_int_t *info );

void   lapackf77_zunglq( const magma_int_t *m, const magma_int_t *n, const magma_int_t *k,
                         cuDoubleComplex *A, const magma_int_t *lda,
                         const cuDoubleComplex *tau,
                         cuDoubleComplex *work, const magma_int_t *lwork,
                         magma_int_t *info );

void   lapackf77_zungql( const magma_int_t *m, const magma_int_t *n, const magma_int_t *k,
                         cuDoubleComplex *A, const magma_int_t *lda,
                         const cuDoubleComplex *tau,
                         cuDoubleComplex *work, const magma_int_t *lwork,
                         magma_int_t *info );

void   lapackf77_zungqr( const magma_int_t *m, const magma_int_t *n, const magma_int_t *k,
                         cuDoubleComplex *A, const magma_int_t *lda,
                         const cuDoubleComplex *tau,
                         cuDoubleComplex *work, const magma_int_t *lwork,
                         magma_int_t *info );

void   lapackf77_zungtr( const char *uplo,
                         const magma_int_t *n,
                         cuDoubleComplex *A, const magma_int_t *lda,
                         const cuDoubleComplex *tau,
                         cuDoubleComplex *work, const magma_int_t *lwork,
                         magma_int_t *info );

void   lapackf77_zunm2r( const char *side, const char *trans,
                         const magma_int_t *m, const magma_int_t *n, const magma_int_t *k,
                         const cuDoubleComplex *A, const magma_int_t *lda,
                         const cuDoubleComplex *tau,
                         cuDoubleComplex *C, const magma_int_t *ldc,
                         cuDoubleComplex *work,
                         magma_int_t *info );

void   lapackf77_zunmbr( const char *vect, const char *side, const char *trans,
                         const magma_int_t *m, const magma_int_t *n, const magma_int_t *k,
                         const cuDoubleComplex *A, const magma_int_t *lda,
                         const cuDoubleComplex *tau,
                         cuDoubleComplex *C, const magma_int_t *ldc,
                         cuDoubleComplex *work, const magma_int_t *lwork,
                         magma_int_t *info );

void   lapackf77_zunmlq( const char *side, const char *trans,
                         const magma_int_t *m, const magma_int_t *n, const magma_int_t *k,
                         const cuDoubleComplex *A, const magma_int_t *lda,
                         const cuDoubleComplex *tau,
                         cuDoubleComplex *C, const magma_int_t *ldc,
                         cuDoubleComplex *work, const magma_int_t *lwork,
                         magma_int_t *info );

void   lapackf77_zunmql( const char *side, const char *trans,
                         const magma_int_t *m, const magma_int_t *n, const magma_int_t *k,
                         const cuDoubleComplex *A, const magma_int_t *lda,
                         const cuDoubleComplex *tau,
                         cuDoubleComplex *C, const magma_int_t *ldc,
                         cuDoubleComplex *work, const magma_int_t *lwork,
                         magma_int_t *info );

void   lapackf77_zunmqr( const char *side, const char *trans,
                         const magma_int_t *m, const magma_int_t *n, const magma_int_t *k,
                         const cuDoubleComplex *A, const magma_int_t *lda,
                         const cuDoubleComplex *tau,
                         cuDoubleComplex *C, const magma_int_t *ldc,
                         cuDoubleComplex *work, const magma_int_t *lwork,
                         magma_int_t *info );

void   lapackf77_zunmtr( const char *side, const char *uplo, const char *trans,
                         const magma_int_t *m, const magma_int_t *n,
                         const cuDoubleComplex *A, const magma_int_t *lda,
                         const cuDoubleComplex *tau,
                         cuDoubleComplex *C, const magma_int_t *ldc,
                         cuDoubleComplex *work, const magma_int_t *lwork,
                         magma_int_t *info );

/*
 * Real precision extras
 */
void   lapackf77_dstebz( const char *range, const char *order,
                         const magma_int_t *n,
                         double *vl, double *vu,
                         magma_int_t *il, magma_int_t *iu,
                         double *abstol,
                         double *d, double *e,
                         const magma_int_t *m, const magma_int_t *nsplit,
                         double *w,
                         magma_int_t *iblock, magma_int_t *isplit,
                         double *work,
                         magma_int_t *iwork,
                         magma_int_t *info );

double lapackf77_dlamc3( double* a, double* b );

void   lapackf77_dlamrg( magma_int_t* n1, magma_int_t* n2,
                         double* a,
                         magma_int_t* dtrd1, magma_int_t* dtrd2, magma_int_t* index );

double lapackf77_dlapy3(double *, double *, double *);

void   lapackf77_dlaed4( magma_int_t* n, magma_int_t* i,
                         double* d,
                         double* z,
                         double* delta,
                         double* rho,
                         double* dlam, magma_int_t* info );

/*
 * Testing functions
 */
#if defined(PRECISION_z) || defined(PRECISION_c)
void   lapackf77_zbdt01( const magma_int_t *m, const magma_int_t *n, const magma_int_t *kd,
                         cuDoubleComplex *A, const magma_int_t *lda,
                         cuDoubleComplex *Q, const magma_int_t *ldq,
                         double *d, double *e,
                         cuDoubleComplex *Pt, const magma_int_t *ldpt,
                         cuDoubleComplex *work,
                         double *rwork,
                         double *resid );

void   lapackf77_zget22( const char *transa, const char *transe, const char *transw, const magma_int_t *n,
                         cuDoubleComplex *A, const magma_int_t *lda,
                         cuDoubleComplex *E, const magma_int_t *lde,
                         cuDoubleComplex *w,
                         cuDoubleComplex *work,
                         double *rwork,
                         double *result );

void   lapackf77_zhet21( const magma_int_t *itype, const char *uplo,
                         const magma_int_t *n, const magma_int_t *kband,
                         cuDoubleComplex *A, const magma_int_t *lda,
                         double *d, double *e,
                         cuDoubleComplex *U, const magma_int_t *ldu,
                         cuDoubleComplex *V, const magma_int_t *ldv,
                         cuDoubleComplex *tau,
                         cuDoubleComplex *work,
                         double *rwork,
                         double *result );

void   lapackf77_zhst01( const magma_int_t *n, const magma_int_t *ilo, const magma_int_t *ihi,
                         cuDoubleComplex *A, const magma_int_t *lda,
                         cuDoubleComplex *H, const magma_int_t *ldh,
                         cuDoubleComplex *Q, const magma_int_t *ldq,
                         cuDoubleComplex *work, const magma_int_t *lwork,
                         double *rwork,
                         double *result );

void   lapackf77_zstt21( const magma_int_t *n, const magma_int_t *kband,
                         double *AD,
                         double *AE,
                         double *SD,
                         double *SE,
                         cuDoubleComplex *U, const magma_int_t *ldu,
                         cuDoubleComplex *work,
                         double *rwork,
                         double *result );

void   lapackf77_zunt01( const char *rowcol, const magma_int_t *m, const magma_int_t *n,
                         cuDoubleComplex *U, const magma_int_t *ldu,
                         cuDoubleComplex *work, const magma_int_t *lwork,
                         double *rwork,
                         double *resid );
#else
void   lapackf77_zbdt01( const magma_int_t *m, const magma_int_t *n, const magma_int_t *kd,
                         cuDoubleComplex *A, const magma_int_t *lda,
                         cuDoubleComplex *Q, const magma_int_t *ldq,
                         double *d, double *e,
                         cuDoubleComplex *Pt, const magma_int_t *ldpt,
                         cuDoubleComplex *work,
                         double *resid );

void   lapackf77_zget22( const char *transa, const char *transe, const char *transw, const magma_int_t *n,
                         cuDoubleComplex *A, const magma_int_t *lda,
                         cuDoubleComplex *E, const magma_int_t *lde,
                         cuDoubleComplex *wr,
                         cuDoubleComplex *wi,
                         double *work,
                         double *result );

void   lapackf77_zhet21( magma_int_t *itype, const char *uplo, const magma_int_t *n, const magma_int_t *kband,
                         cuDoubleComplex *A, const magma_int_t *lda,
                         double *d, double *e,
                         cuDoubleComplex *U, const magma_int_t *ldu,
                         cuDoubleComplex *V, const magma_int_t *ldv,
                         cuDoubleComplex *tau,
                         cuDoubleComplex *work,
                         double *result );

void   lapackf77_zhst01( const magma_int_t *n, const magma_int_t *ilo, const magma_int_t *ihi,
                         cuDoubleComplex *A, const magma_int_t *lda,
                         cuDoubleComplex *H, const magma_int_t *ldh,
                         cuDoubleComplex *Q, const magma_int_t *ldq,
                         cuDoubleComplex *work, const magma_int_t *lwork,
                         double *result );

void   lapackf77_zstt21( const magma_int_t *n, const magma_int_t *kband,
                         double *AD,
                         double *AE,
                         double *SD,
                         double *SE,
                         cuDoubleComplex *U, const magma_int_t *ldu,
                         cuDoubleComplex *work,
                         double *result );

void   lapackf77_zunt01( const char *rowcol, const magma_int_t *m, const magma_int_t *n,
                         cuDoubleComplex *U, const magma_int_t *ldu,
                         cuDoubleComplex *work, const magma_int_t *lwork,
                         double *resid );
#endif

void   lapackf77_zlarfy( const char *uplo, const magma_int_t *n,
                         cuDoubleComplex *V, const magma_int_t *incv,
                         cuDoubleComplex *tau,
                         cuDoubleComplex *C, const magma_int_t *ldc,
                         cuDoubleComplex *work );

void   lapackf77_zlarfx( const char *side, const magma_int_t *m, const magma_int_t *n,
                         cuDoubleComplex *V,
                         cuDoubleComplex *tau,
                         cuDoubleComplex *C, const magma_int_t *ldc,
                         cuDoubleComplex *work );

double lapackf77_zqpt01( const magma_int_t *m, const magma_int_t *n, const magma_int_t *k,
                         cuDoubleComplex *A,
                         cuDoubleComplex *Af, const magma_int_t *lda,
                         cuDoubleComplex *tau, magma_int_t *jpvt,
                         cuDoubleComplex *work, const magma_int_t *lwork );

void   lapackf77_zqrt02( const magma_int_t *m, const magma_int_t *n, const magma_int_t *k,
                         cuDoubleComplex *A,
                         cuDoubleComplex *AF,
                         cuDoubleComplex *Q,
                         cuDoubleComplex *R, const magma_int_t *lda,
                         cuDoubleComplex *tau,
                         cuDoubleComplex *work, const magma_int_t *lwork,
                         double *rwork,
                         double *result );

#ifdef __cplusplus
}
#endif

#undef DWORKFORZ
#undef DWORKFORZ_AND_LD
#undef WSPLIT
#undef PRECISION_z

#endif /* MAGMA_ZLAPACK_H */
