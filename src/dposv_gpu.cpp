/*
    -- MAGMA (version 1.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2011

       @generated d Sun Nov 13 20:48:14 2011

*/
#include "common_magma.h"

// === Define what BLAS to use ============================================
#define PRECISION_d
#if (defined(PRECISION_s) || defined(PRECISION_d)) 
  #define cublasDtrsm magmablas_dtrsm
#endif
// === End defining what BLAS to use =======================================

extern "C" magma_int_t
magma_dposv_gpu( char uplo, magma_int_t n, magma_int_t nrhs, 
                 double *dA, magma_int_t ldda, 
                 double *dB, magma_int_t lddb, magma_int_t *info )
{
/*  -- magma (version 1.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2011
 
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

    dA      (input/output) DOUBLE_PRECISION array on the GPU, dimension (LDDA,N)   
            On entry, the symmetric matrix dA.  If UPLO = 'U', the leading   
            N-by-N upper triangular part of dA contains the upper   
            triangular part of the matrix dA, and the strictly lower   
            triangular part of dA is not referenced.  If UPLO = 'L', the   
            leading N-by-N lower triangular part of dA contains the lower   
            triangular part of the matrix dA, and the strictly upper   
            triangular part of dA is not referenced.   

            On exit, if INFO = 0, the factor U or L from the Cholesky   
            factorization dA = U**T*U or dA = L*L**T.   

    LDDA    (input) INTEGER
            The leading dimension of the array A.  LDA >= max(1,N).

    dB      (input/output) DOUBLE_PRECISION array on the GPU, dimension (LDB,NRHS)
            On entry, the right hand side matrix B.
            On exit, the solution matrix X.

    LDDB    (input) INTEGER
            The leading dimension of the array B.  LDB >= max(1,N).

    INFO    (output) INTEGER
            = 0:  successful exit
            < 0:  if INFO = -i, the i-th argument had an illegal value
    =====================================================================   */

    magma_int_t ret;
    
    *info = 0 ; 
    if( (uplo != 'U') && (uplo != 'u') && (uplo != 'L') && (uplo != 'l') )
        *info = -1; 
    if( n < 0 )
        *info = -2; 
    if( nrhs < 0) 
        *info = -3; 
    if ( ldda < max(1, n) )
        *info = -5; 
    if ( lddb < max(1, n) )
        *info = -7;
    if (*info != 0) {
        magma_xerbla( __func__, -(*info) );
        return MAGMA_ERR_ILLEGAL_VALUE;
    }

    /* Quick return if possible */
    if ( (n == 0) || (nrhs == 0) ) {
        return MAGMA_SUCCESS;
    }

    ret = magma_dpotrf_gpu(uplo, n, dA, ldda, info);
    if ( (ret != MAGMA_SUCCESS) || ( *info != 0 ) ) {
        return ret;
    }

    ret = magma_dpotrs_gpu(uplo, n, nrhs, dA, ldda, dB, lddb, info);
    if ( (ret != MAGMA_SUCCESS) || ( *info != 0 ) ) {
        return ret;
    }

    return MAGMA_SUCCESS;
}
