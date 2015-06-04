/*
    -- MAGMA (version 1.4.1-beta2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       December 2013

       @generated d Mon Dec  9 16:19:08 2013

*/
#include "common_magma.h"

extern "C" magma_int_t
magma_dposv    ( char uplo, magma_int_t n, magma_int_t nrhs,
                 double *A, magma_int_t lda,
                 double *B, magma_int_t ldb, magma_int_t *info )
{
/*  -- MAGMA (version 1.4.1-beta2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       December 2013

    Purpose
    =======
    DPOSV computes the solution to a real system of linear equations
       A * X = B,
    where A is an N-by-N symmetric positive definite matrix and X and B
    are N-by-NRHS matrices.
    The Cholesky decomposition is used to factor A as
       A = U**T * U,  if UPLO = 'U', or
       A = L * L**T,  if UPLO = 'L',
    where U is an upper triangular matrix and  L is a lower triangular
    matrix.  The factored form of A is then used to solve the system of
    equations A * X = B.

    Arguments
    =========
    UPLO    (input) CHARACTER*1
            = 'U':  Upper triangle of A is stored;
            = 'L':  Lower triangle of A is stored.

    N       (input) INTEGER
            The order of the matrix A.  N >= 0.

    NRHS    (input) INTEGER
            The number of right hand sides, i.e., the number of columns
            of the matrix B.  NRHS >= 0.

    A       (input/output) DOUBLE_PRECISION array, dimension (LDA,N)
            On entry, the symmetric matrix A.  If UPLO = 'U', the leading
            N-by-N upper triangular part of A contains the upper
            triangular part of the matrix A, and the strictly lower
            triangular part of A is not referenced.  If UPLO = 'L', the
            leading N-by-N lower triangular part of A contains the lower
            triangular part of the matrix A, and the strictly upper
            triangular part of A is not referenced.

            On exit, if INFO = 0, the factor U or L from the Cholesky
            factorization A = U**T*U or A = L*L**T.

    LDA     (input) INTEGER
            The leading dimension of the array A.  LDA >= max(1,N).

    B       (input/output) DOUBLE_PRECISION array, dimension (LDB,NRHS)
            On entry, the right hand side matrix B.
            On exit, the solution matrix X.

    LDB     (input) INTEGER
            The leading dimension of the array B.  LDB >= max(1,N).

    INFO    (output) INTEGER
            = 0:  successful exit
            < 0:  if INFO = -i, the i-th argument had an illegal value
    =====================================================================   */

    magma_int_t num_gpus, ldda, lddb;

    *info = 0 ;
    if( (uplo != 'U') && (uplo != 'u') && (uplo != 'L') && (uplo != 'l') )
        *info = -1;
    if( n < 0 )
        *info = -2;
    if( nrhs < 0)
        *info = -3;
    if ( lda < max(1, n) )
        *info = -5;
    if ( ldb < max(1, n) )
        *info = -7;
    if (*info != 0) {
        magma_xerbla( __func__, -(*info) );
        return *info;
    }

    /* Quick return if possible */
    if ( (n==0) || (nrhs == 0) ) {
        return *info;
    }

    /* If single-GPU and allocation suceeds, use GPU interface. */
    num_gpus = magma_num_gpus();
    double *dA, *dB;
    if ( num_gpus > 1 ) {
        goto CPU_INTERFACE;
    }
    ldda = ((n+31)/32)*32;
    lddb = ldda;
    if ( MAGMA_SUCCESS != magma_dmalloc( &dA, ldda*n )) {
        goto CPU_INTERFACE;
    }
    if ( MAGMA_SUCCESS != magma_dmalloc( &dB, lddb*nrhs )) {
        magma_free( dA );
        goto CPU_INTERFACE;
    }
    magma_dsetmatrix( n, n, A, lda, dA, ldda );
    magma_dpotrf_gpu( uplo, n, dA, ldda, info );
    if ( *info == MAGMA_ERR_DEVICE_ALLOC ) {
        magma_free( dA );
        magma_free( dB );
        goto CPU_INTERFACE;
    }
    magma_dgetmatrix( n, n, dA, ldda, A, lda );
    if ( *info == 0 ) {
        magma_dsetmatrix( n, nrhs, B, ldb, dB, lddb );
        magma_dpotrs_gpu( uplo, n, nrhs, dA, ldda, dB, lddb, info );
        magma_dgetmatrix( n, nrhs, dB, lddb, B, ldb );
    }
    magma_free( dA );
    magma_free( dB );
    return *info;

CPU_INTERFACE:
    /* If multi-GPU or allocation failed, use CPU interface and LAPACK.
     * Faster to use LAPACK for potrs than to copy A to GPU. */
    magma_dpotrf( uplo, n, A, lda, info );
    if ( *info == 0 ) {
        lapackf77_dpotrs( &uplo, &n, &nrhs, A, &lda, B, &ldb, info );
    }

    return *info;
}
