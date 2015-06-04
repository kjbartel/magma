/*
 *   -- MAGMA (version 1.3.0) --
 *      Univ. of Tennessee, Knoxville
 *      Univ. of California, Berkeley
 *      Univ. of Colorado, Denver
 *      November 2012
 *
 * @author Mark Gates
 * @generated s Wed Nov 14 22:53:36 2012
 */

#include <stdlib.h>
#include <stdio.h>

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
void magma_sswap(
    magma_int_t n,
    float *dx, magma_int_t incx,
    float *dy, magma_int_t incy )
{
    cublasSswap( n, dx, incx, dy, incy );
}

// --------------------
extern "C"
magma_int_t magma_isamax(
    magma_int_t n,
    float *dx, magma_int_t incx )
{
    return cublasIsamax( n, dx, incx );
}

// ========================================
// Level 2 BLAS

// --------------------
extern "C"
void magma_sgemv(
    magma_trans_t transA,
    magma_int_t m, magma_int_t n,
    float alpha, float const* dA, magma_int_t lda,
                           float const* dx, magma_int_t incx,
    float beta,  float*       dy, magma_int_t incy )
{
    cublasSgemv(
        cublas_trans_const( transA ),
        m, n,
        alpha, dA, lda,
               dx, incx,
        beta,  dy, incy );
}

// --------------------
extern "C"
void magma_ssymv(
    magma_uplo_t uplo,
    magma_int_t n,
    float alpha, float const* dA, magma_int_t lda,
                           float const* dx, magma_int_t incx,
    float beta,  float*       dy, magma_int_t incy )
{
    cublasSsymv(
        cublas_uplo_const( uplo ),
        n,
        alpha, dA, lda,
               dx, incx,
        beta,  dy, incy );
}

// --------------------
extern "C"
void magma_strsv(
    magma_uplo_t uplo, magma_trans_t trans, magma_diag_t diag, 
    magma_int_t n, 
    float const *dA, magma_int_t lda, 
    float       *dx, magma_int_t incx )
{
    cublasStrsv(
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
void magma_sgemm(
    magma_trans_t transA, magma_trans_t transB,
    magma_int_t m, magma_int_t n, magma_int_t k,
    float alpha, float const* dA, magma_int_t lda,
                           float const* dB, magma_int_t ldb,
    float beta,  float*       dC, magma_int_t ldc )
{
    cublasSgemm(
        cublas_trans_const( transA ),
        cublas_trans_const( transB ),
        m, n, k,
        alpha, dA, lda,
               dB, ldb,
        beta,  dC, ldc );
}

// --------------------
extern "C"
void magma_ssymm(
    magma_side_t side, magma_uplo_t uplo,
    magma_int_t m, magma_int_t n,
    float alpha, float const* dA, magma_int_t lda,
                           float const* dB, magma_int_t ldb,
    float beta,  float*       dC, magma_int_t ldc )
{
    cublasSsymm(
        cublas_side_const( side ),
        cublas_uplo_const( uplo ),
        m, n,
        alpha, dA, lda,
               dB, ldb,
        beta,  dC, ldc );
}

// --------------------
extern "C"
void magma_ssyrk(
    magma_uplo_t uplo, magma_trans_t trans,
    magma_int_t n, magma_int_t k,
    float alpha, float const* dA, magma_int_t lda,
    float beta,  float*       dC, magma_int_t ldc )
{
    cublasSsyrk(
        cublas_uplo_const( uplo ),
        cublas_trans_const( trans ),
        n, k,
        alpha, dA, lda,
        beta,  dC, ldc );
}

// --------------------
extern "C"
void magma_ssyr2k(
    magma_uplo_t uplo, magma_trans_t trans,
    magma_int_t n, magma_int_t k,
    float alpha, float const *dA, magma_int_t lda,
                           float const *dB, magma_int_t ldb,
    float beta,           float       *dC, magma_int_t ldc )
{
    cublasSsyr2k(
        cublas_uplo_const( uplo ),
        cublas_trans_const( trans ),
        n, k,
        alpha, dA, lda,
               dB, ldb,
        beta,  dC, ldc );
}

// --------------------
extern "C"
void magma_strmm(
    magma_side_t side, magma_uplo_t uplo, magma_trans_t trans, magma_diag_t diag,
    magma_int_t m, magma_int_t n,
    float alpha, float const *dA, magma_int_t lda,
                           float       *dB, magma_int_t ldb )
{
    cublasStrmm(
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
void magma_strsm(
    magma_side_t side, magma_uplo_t uplo, magma_trans_t trans, magma_diag_t diag,
    magma_int_t m, magma_int_t n,
    float alpha, float const* dA, magma_int_t lda,
                           float*       dB, magma_int_t ldb )
{
    cublasStrsm(
        cublas_side_const( side ),
        cublas_uplo_const( uplo ),
        cublas_trans_const( trans ),
        cublas_diag_const( diag ),
        m, n,
        alpha, dA, lda,
               dB, ldb );
}

#endif // HAVE_CUBLAS
