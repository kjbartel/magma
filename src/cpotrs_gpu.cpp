/*
    -- MAGMA (version 1.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2011

       @generated c Sun Nov 13 20:48:14 2011

*/
#include "common_magma.h"

// === Define what BLAS to use ============================================
#define PRECISION_c
#if (defined(PRECISION_s) || defined(PRECISION_d)) 
  #define cublasCtrsm magmablas_ctrsm
#endif
// === End defining what BLAS to use =======================================

extern "C" magma_int_t
magma_cpotrs_gpu(char uplo, magma_int_t n, magma_int_t nrhs, 
                 cuFloatComplex *dA, magma_int_t ldda, 
                 cuFloatComplex *dB, magma_int_t lddb, magma_int_t *info)
{
/*  -- magma (version 1.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2011
 
    Purpose
    =======

    CPOTRS solves a system of linear equations A*X = B with a Hermitian
    positive definite matrix A using the Cholesky factorization
    A = U**H*U or A = L*L**H computed by CPOTRF.

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

    dA      (input) COMPLEX array on the GPU, dimension (LDDA,N)
            The triangular factor U or L from the Cholesky factorization
            A = U**H*U or A = L*L**H, as computed by CPOTRF.

    LDDA    (input) INTEGER
            The leading dimension of the array A.  LDDA >= max(1,N).

    dB      (input/output) COMPLEX array on the GPU, dimension (LDDB,NRHS)
            On entry, the right hand side matrix B.
            On exit, the solution matrix X.

    LDDB    (input) INTEGER
            The leading dimension of the array B.  LDDB >= max(1,N).

    INFO    (output) INTEGER
            = 0:  successful exit
            < 0:  if INFO = -i, the i-th argument had an illegal value
    =====================================================================   */

    cuFloatComplex c_one = MAGMA_C_ONE;
    
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

    if( (uplo=='U') || (uplo=='u') ){
        if ( nrhs == 1) {
            cublasCtrsv(MagmaUpper, MagmaConjTrans, MagmaNonUnit, n, dA, ldda, dB, 1 );
            cublasCtrsv(MagmaUpper, MagmaNoTrans,   MagmaNonUnit, n, dA, ldda, dB, 1 );
        } else {
            cublasCtrsm(MagmaLeft, MagmaUpper, MagmaConjTrans, MagmaNonUnit, n, nrhs, c_one, dA, ldda, dB, lddb);
            cublasCtrsm(MagmaLeft, MagmaUpper, MagmaNoTrans,   MagmaNonUnit, n, nrhs, c_one, dA, ldda, dB, lddb);
        }
    }
    else{
        if ( nrhs == 1) {
            cublasCtrsv(MagmaLower, MagmaNoTrans,   MagmaNonUnit, n, dA, ldda, dB, 1 );
            cublasCtrsv(MagmaLower, MagmaConjTrans, MagmaNonUnit, n, dA, ldda, dB, 1 );
        } else {
            cublasCtrsm(MagmaLeft, MagmaLower, MagmaNoTrans,   MagmaNonUnit, n, nrhs, c_one, dA, ldda, dB, lddb);
            cublasCtrsm(MagmaLeft, MagmaLower, MagmaConjTrans, MagmaNonUnit, n, nrhs, c_one, dA, ldda, dB, lddb);
        }
    }

    return MAGMA_SUCCESS;
}
