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

// ========================================
// copying vectors
extern "C"
void magma_csetvector(
    magma_int_t n,
    cuFloatComplex const* hx_src, magma_int_t incx,
    cuFloatComplex*       dy_dst, magma_int_t incy )
{
    cublasStatus_t status;
    status = cublasSetVector(
        n, sizeof(cuFloatComplex),
        hx_src, incx,
        dy_dst, incy );
    check_error( status );
}

// --------------------
extern "C"
void magma_cgetvector(
    magma_int_t n,
    cuFloatComplex const* dx_src, magma_int_t incx,
    cuFloatComplex*       hy_dst, magma_int_t incy )
{
    cublasStatus_t status;
    status = cublasGetVector(
        n, sizeof(cuFloatComplex),
        dx_src, incx,
        hy_dst, incy );
    check_error( status );
}

// --------------------
extern "C"
void magma_csetvector_async(
    magma_int_t n,
    cuFloatComplex const* hx_src, magma_int_t incx,
    cuFloatComplex*       dy_dst, magma_int_t incy,
    cudaStream_t stream )
{
    cublasStatus_t status;
    status = cublasSetVectorAsync(
        n, sizeof(cuFloatComplex),
        hx_src, incx,
        dy_dst, incy, stream );
    check_error( status );
}

// --------------------
extern "C"
void magma_cgetvector_async(
    magma_int_t n,
    cuFloatComplex const* dx_src, magma_int_t incx,
    cuFloatComplex*       hy_dst, magma_int_t incy,
    cudaStream_t stream )
{
    cublasStatus_t status;
    status = cublasGetVectorAsync(
        n, sizeof(cuFloatComplex),
        dx_src, incx,
        hy_dst, incy, stream );
    check_error( status );
}


// ========================================
// copying sub-matrices (contiguous columns)
extern "C"
void magma_csetmatrix(
    magma_int_t m, magma_int_t n,
    cuFloatComplex const* hA_src, magma_int_t lda,
    cuFloatComplex*       dB_dst, magma_int_t ldb )
{
    cublasStatus_t status;
    status = cublasSetMatrix(
        m, n, sizeof(cuFloatComplex),
        hA_src, lda,
        dB_dst, ldb );
    check_error( status );
}

// --------------------
extern "C"
void magma_cgetmatrix(
    magma_int_t m, magma_int_t n,
    cuFloatComplex const* dA_src, magma_int_t lda,
    cuFloatComplex*       hB_dst, magma_int_t ldb )
{
    cublasStatus_t status;
    status = cublasGetMatrix(
        m, n, sizeof(cuFloatComplex),
        dA_src, lda,
        hB_dst, ldb );
    check_error( status );
}

// --------------------
extern "C"
void magma_csetmatrix_async(
    magma_int_t m, magma_int_t n,
    cuFloatComplex const* hA_src, magma_int_t lda,
    cuFloatComplex*       dB_dst, magma_int_t ldb,
    cudaStream_t stream )
{
    cublasStatus_t status;
    status = cublasSetMatrixAsync(
        m, n, sizeof(cuFloatComplex),
        hA_src, lda,
        dB_dst, ldb, stream );
    check_error( status );
}

// --------------------
extern "C"
void magma_cgetmatrix_async(
    magma_int_t m, magma_int_t n,
    cuFloatComplex const* dA_src, magma_int_t lda,
    cuFloatComplex*       hB_dst, magma_int_t ldb,
    cudaStream_t stream )
{
    cublasStatus_t status;
    status = cublasGetMatrixAsync(
        m, n, sizeof(cuFloatComplex),
        dA_src, lda,
        hB_dst, ldb, stream );
    check_error( status );
}

// --------------------
extern "C"
void magma_ccopymatrix(
    magma_int_t m, magma_int_t n,
    cuFloatComplex const* dA_src, magma_int_t lda,
    cuFloatComplex*       dB_dst, magma_int_t ldb )
{
    cudaError_t status;
    status = cudaMemcpy2D(
        dB_dst, ldb*sizeof(cuFloatComplex),
        dA_src, lda*sizeof(cuFloatComplex),
        m*sizeof(cuFloatComplex), n, cudaMemcpyDeviceToDevice );
    assert( status == cudaSuccess );
}

// --------------------
extern "C"
void magma_ccopymatrix_async(
    magma_int_t m, magma_int_t n,
    cuFloatComplex const* dA_src, magma_int_t lda,
    cuFloatComplex*       dB_dst, magma_int_t ldb,
    cudaStream_t stream )
{
    cudaError_t status;
    status = cudaMemcpy2DAsync(
        dB_dst, ldb*sizeof(cuFloatComplex),
        dA_src, lda*sizeof(cuFloatComplex),
        m*sizeof(cuFloatComplex), n, cudaMemcpyDeviceToDevice, stream );
    assert( status == cudaSuccess );
}

#endif // HAVE_CUBLAS
