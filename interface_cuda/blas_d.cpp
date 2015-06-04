/*
 *   -- MAGMA (version 1.2.1) --
 *      Univ. of Tennessee, Knoxville
 *      Univ. of California, Berkeley
 *      Univ. of Colorado, Denver
 *      June 2012
 *
 * @author Mark Gates
 * @generated d Thu Jun 28 12:31:07 2012
 */

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "magma.h"
#include "error.h"

#ifdef HAVE_CUBLAS

// For now, magma constants are the same as cublas v1 constants (character).
// This will change in the future.
#define cublas_side_const(  x )  (x)
#define cublas_uplo_const(  x )  (x)
#define cublas_trans_const( x )  (x)
#define cublas_diag_const(  x )  (x)

// ========================================
// Level 1 BLAS

// --------------------
extern "C"
void magma_dswap(
    magma_int_t n,
    double *dx, magma_int_t incx,
    double *dy, magma_int_t incy )
{
    cublasDswap( n, dx, incx, dy, incy );
}

// --------------------
extern "C"
magma_int_t magma_idamax(
    magma_int_t n,
    double *dx, magma_int_t incx )
{
    return cublasIdamax( n, dx, incx );
}

// ========================================
// Level 2 BLAS

// --------------------
extern "C"
void magma_dgemv(
    magma_trans_t transA,
    magma_int_t m, magma_int_t n,
    double alpha, double const* dA, magma_int_t lda,
                           double const* dx, magma_int_t incx,
    double beta,  double*       dy, magma_int_t incy )
{
    cublasDgemv(
        cublas_trans_const( transA ),
        m, n,
        alpha, dA, lda,
               dx, incx,
        beta,  dy, incy );
}

// --------------------
extern "C"
void magma_dsymv(
    magma_uplo_t uplo,
    magma_int_t n,
    double alpha, double const* dA, magma_int_t lda,
                           double const* dx, magma_int_t incx,
    double beta,  double*       dy, magma_int_t incy )
{
    cublasDsymv(
        cublas_uplo_const( uplo ),
        n,
        alpha, dA, lda,
               dx, incx,
        beta,  dy, incy );
}

// --------------------
extern "C"
void magma_dtrsv(
    magma_uplo_t uplo, magma_trans_t trans, magma_diag_t diag, 
    magma_int_t n, 
    double const *dA, magma_int_t lda, 
    double       *dx, magma_int_t incx )
{
    cublasDtrsv(
        cublas_uplo_const( uplo ),
        cublas_trans_const( trans ),
        cublas_diag_const( diag ),
        n,
        dA, lda,
        dx, incx );
}

// ========================================
// Level 3 BLAS

// --------------------
extern "C"
void magma_dgemm(
    magma_trans_t transA, magma_trans_t transB,
    magma_int_t m, magma_int_t n, magma_int_t k,
    double alpha, double const* dA, magma_int_t lda,
                           double const* dB, magma_int_t ldb,
    double beta,  double*       dC, magma_int_t ldc )
{
    cublasDgemm(
        cublas_trans_const( transA ),
        cublas_trans_const( transB ),
        m, n, k,
        alpha, dA, lda,
               dB, ldb,
        beta,  dC, ldc );
}

// --------------------
extern "C"
void magma_dsymm(
    magma_side_t side, magma_uplo_t uplo,
    magma_int_t m, magma_int_t n,
    double alpha, double const* dA, magma_int_t lda,
                           double const* dB, magma_int_t ldb,
    double beta,  double*       dC, magma_int_t ldc )
{
    cublasDsymm(
        cublas_side_const( side ),
        cublas_uplo_const( uplo ),
        m, n,
        alpha, dA, lda,
               dB, ldb,
        beta,  dC, ldc );
}

// --------------------
extern "C"
void magma_dsyrk(
    magma_uplo_t uplo, magma_trans_t trans,
    magma_int_t n, magma_int_t k,
    double alpha, double const* dA, magma_int_t lda,
    double beta,  double*       dC, magma_int_t ldc )
{
    cublasDsyrk(
        cublas_uplo_const( uplo ),
        cublas_trans_const( trans ),
        n, k,
        alpha, dA, lda,
        beta,  dC, ldc );
}

// --------------------
extern "C"
void magma_dsyr2k(
    magma_uplo_t uplo, magma_trans_t trans,
    magma_int_t n, magma_int_t k,
    double alpha, double const *dA, magma_int_t lda,
                           double const *dB, magma_int_t ldb,
    double beta,           double       *dC, magma_int_t ldc )
{
    cublasDsyr2k(
        cublas_uplo_const( uplo ),
        cublas_trans_const( trans ),
        n, k,
        alpha, dA, lda,
               dB, ldb,
        beta,  dC, ldc );
}

// --------------------
extern "C"
void magma_dtrmm(
    magma_side_t side, magma_uplo_t uplo, magma_trans_t trans, magma_diag_t diag,
    magma_int_t m, magma_int_t n,
    double alpha, double const *dA, magma_int_t lda,
                           double       *dB, magma_int_t ldb )
{
    cublasDtrmm(
        cublas_side_const( side ),
        cublas_uplo_const( uplo ),
        cublas_trans_const( trans ),
        cublas_diag_const( diag ),
        m, n,
        alpha, dA, lda,
               dB, ldb );
}

// --------------------
extern "C"
void magma_dtrsm(
    magma_side_t side, magma_uplo_t uplo, magma_trans_t trans, magma_diag_t diag,
    magma_int_t m, magma_int_t n,
    double alpha, double const* dA, magma_int_t lda,
                           double*       dB, magma_int_t ldb )
{
    cublasDtrsm(
        cublas_side_const( side ),
        cublas_uplo_const( uplo ),
        cublas_trans_const( trans ),
        cublas_diag_const( diag ),
        m, n,
        alpha, dA, lda,
               dB, ldb );
}

#endif // HAVE_CUBLAS
