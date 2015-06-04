/*
 *   -- MAGMA (version 1.2.0) --
 *      Univ. of Tennessee, Knoxville
 *      Univ. of California, Berkeley
 *      Univ. of Colorado, Denver
 *      April 2012
 *
 * @author Mark Gates
 * @generated c Tue May 15 18:17:56 2012
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
void magma_cswap(
    magma_int_t n,
    cuFloatComplex *dx, magma_int_t incx,
    cuFloatComplex *dy, magma_int_t incy )
{
    cublasCswap( n, dx, incx, dy, incy );
}

// --------------------
extern "C"
magma_int_t magma_icamax(
    magma_int_t n,
    cuFloatComplex *dx, magma_int_t incx )
{
    return cublasIcamax( n, dx, incx );
}

// ========================================
// Level 2 BLAS

// --------------------
extern "C"
void magma_cgemv(
    magma_trans_t transA,
    magma_int_t m, magma_int_t n,
    cuFloatComplex alpha, cuFloatComplex const* dA, magma_int_t lda,
                           cuFloatComplex const* dx, magma_int_t incx,
    cuFloatComplex beta,  cuFloatComplex*       dy, magma_int_t incy )
{
    cublasCgemv(
        cublas_trans_const( transA ),
        m, n,
        alpha, dA, lda,
               dx, incx,
        beta,  dy, incy );
}

// --------------------
extern "C"
void magma_chemv(
    magma_uplo_t uplo,
    magma_int_t n,
    cuFloatComplex alpha, cuFloatComplex const* dA, magma_int_t lda,
                           cuFloatComplex const* dx, magma_int_t incx,
    cuFloatComplex beta,  cuFloatComplex*       dy, magma_int_t incy )
{
    cublasChemv(
        cublas_uplo_const( uplo ),
        n,
        alpha, dA, lda,
               dx, incx,
        beta,  dy, incy );
}

// --------------------
extern "C"
void magma_ctrsv(
    magma_uplo_t uplo, magma_trans_t trans, magma_diag_t diag, 
    magma_int_t n, 
    cuFloatComplex const *dA, magma_int_t lda, 
    cuFloatComplex       *dx, magma_int_t incx )
{
    cublasCtrsv(
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
void magma_cgemm(
    magma_trans_t transA, magma_trans_t transB,
    magma_int_t m, magma_int_t n, magma_int_t k,
    cuFloatComplex alpha, cuFloatComplex const* dA, magma_int_t lda,
                           cuFloatComplex const* dB, magma_int_t ldb,
    cuFloatComplex beta,  cuFloatComplex*       dC, magma_int_t ldc )
{
    cublasCgemm(
        cublas_trans_const( transA ),
        cublas_trans_const( transB ),
        m, n, k,
        alpha, dA, lda,
               dB, ldb,
        beta,  dC, ldc );
}

// --------------------
extern "C"
void magma_chemm(
    magma_side_t side, magma_uplo_t uplo,
    magma_int_t m, magma_int_t n,
    cuFloatComplex alpha, cuFloatComplex const* dA, magma_int_t lda,
                           cuFloatComplex const* dB, magma_int_t ldb,
    cuFloatComplex beta,  cuFloatComplex*       dC, magma_int_t ldc )
{
    cublasChemm(
        cublas_side_const( side ),
        cublas_uplo_const( uplo ),
        m, n,
        alpha, dA, lda,
               dB, ldb,
        beta,  dC, ldc );
}

// --------------------
extern "C"
void magma_cherk(
    magma_uplo_t uplo, magma_trans_t trans,
    magma_int_t n, magma_int_t k,
    float alpha, cuFloatComplex const* dA, magma_int_t lda,
    float beta,  cuFloatComplex*       dC, magma_int_t ldc )
{
    cublasCherk(
        cublas_uplo_const( uplo ),
        cublas_trans_const( trans ),
        n, k,
        alpha, dA, lda,
        beta,  dC, ldc );
}

// --------------------
extern "C"
void magma_cher2k(
    magma_uplo_t uplo, magma_trans_t trans,
    magma_int_t n, magma_int_t k,
    cuFloatComplex alpha, cuFloatComplex const *dA, magma_int_t lda,
                           cuFloatComplex const *dB, magma_int_t ldb,
    float beta,           cuFloatComplex       *dC, magma_int_t ldc )
{
    cublasCher2k(
        cublas_uplo_const( uplo ),
        cublas_trans_const( trans ),
        n, k,
        alpha, dA, lda,
               dB, ldb,
        beta,  dC, ldc );
}

// --------------------
extern "C"
void magma_ctrmm(
    magma_side_t side, magma_uplo_t uplo, magma_trans_t trans, magma_diag_t diag,
    magma_int_t m, magma_int_t n,
    cuFloatComplex alpha, cuFloatComplex const *dA, magma_int_t lda,
                           cuFloatComplex       *dB, magma_int_t ldb )
{
    cublasCtrmm(
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
void magma_ctrsm(
    magma_side_t side, magma_uplo_t uplo, magma_trans_t trans, magma_diag_t diag,
    magma_int_t m, magma_int_t n,
    cuFloatComplex alpha, cuFloatComplex const* dA, magma_int_t lda,
                           cuFloatComplex*       dB, magma_int_t ldb )
{
    cublasCtrsm(
        cublas_side_const( side ),
        cublas_uplo_const( uplo ),
        cublas_trans_const( trans ),
        cublas_diag_const( diag ),
        m, n,
        alpha, dA, lda,
               dB, ldb );
}

#endif // HAVE_CUBLAS
