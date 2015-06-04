/*
 *   -- MAGMA (version 1.2.0) --
 *      Univ. of Tennessee, Knoxville
 *      Univ. of California, Berkeley
 *      Univ. of Colorado, Denver
 *      April 2012
 *
 * @author Mark Gates
 * @precisions normal z -> s d c
 */

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "magma.h"

#ifdef HAVE_CUBLAS

// generic, type-independent routines to copy data.
// type-safe versions which avoid the user needing sizeof(...) are in copy_[sdcz].cpp

// ========================================
// copying vectors
extern "C"
void magma_setvector(
    magma_int_t n, size_t elemSize,
    void const* hx_src, magma_int_t incx,
    void*       dy_dst, magma_int_t incy )
{
    cublasStatus_t status;
    status = cublasSetVector(
        n, elemSize,
        hx_src, incx,
        dy_dst, incy );
    assert( status == CUBLAS_STATUS_SUCCESS );
}

// --------------------
extern "C"
void magma_getvector(
    magma_int_t n, size_t elemSize,
    void const* dx_src, magma_int_t incx,
    void*       hy_dst, magma_int_t incy )
{
    cublasStatus_t status;
    status = cublasGetVector(
        n, elemSize,
        dx_src, incx,
        hy_dst, incy );
    assert( status == CUBLAS_STATUS_SUCCESS );
}

// --------------------
extern "C"
void magma_setvector_async(
    magma_int_t n, size_t elemSize,
    void const* hx_src, magma_int_t incx,
    void*       dy_dst, magma_int_t incy,
    cudaStream_t stream )
{
    cublasStatus_t status;
    status = cublasSetVectorAsync(
        n, elemSize,
        hx_src, incx,
        dy_dst, incy, stream );
    assert( status == CUBLAS_STATUS_SUCCESS );
}

// --------------------
extern "C"
void magma_getvector_async(
    magma_int_t n, size_t elemSize,
    void const* dx_src, magma_int_t incx,
    void*       hy_dst, magma_int_t incy,
    cudaStream_t stream )
{
    cublasStatus_t status;
    status = cublasGetVectorAsync(
        n, elemSize,
        dx_src, incx,
        hy_dst, incy, stream );
    assert( status == CUBLAS_STATUS_SUCCESS );
}


// ========================================
// copying sub-matrices (contiguous columns)
extern "C"
void magma_setmatrix(
    magma_int_t m, magma_int_t n, size_t elemSize,
    void const* hA_src, magma_int_t lda,
    void*       dB_dst, magma_int_t ldb )
{
    cublasStatus_t status;
    status = cublasSetMatrix(
        m, n, elemSize,
        hA_src, lda,
        dB_dst, ldb );
    assert( status == CUBLAS_STATUS_SUCCESS );
}

// --------------------
extern "C"
void magma_getmatrix(
    magma_int_t m, magma_int_t n, size_t elemSize,
    void const* dA_src, magma_int_t lda,
    void*       hB_dst, magma_int_t ldb )
{
    cublasStatus_t status;
    status = cublasGetMatrix(
        m, n, elemSize,
        dA_src, lda,
        hB_dst, ldb );
    assert( status == CUBLAS_STATUS_SUCCESS );
}

// --------------------
extern "C"
void magma_setmatrix_async(
    magma_int_t m, magma_int_t n, size_t elemSize,
    void const* hA_src, magma_int_t lda,
    void*       dB_dst, magma_int_t ldb,
    cudaStream_t stream )
{
    cublasStatus_t status;
    status = cublasSetMatrixAsync(
        m, n, elemSize,
        hA_src, lda,
        dB_dst, ldb, stream );
    assert( status == CUBLAS_STATUS_SUCCESS );
}

// --------------------
extern "C"
void magma_getmatrix_async(
    magma_int_t m, magma_int_t n, size_t elemSize,
    void const* dA_src, magma_int_t lda,
    void*       hB_dst, magma_int_t ldb,
    cudaStream_t stream )
{
    cublasStatus_t status;
    status = cublasGetMatrixAsync(
        m, n, elemSize,
        dA_src, lda,
        hB_dst, ldb, stream );
    assert( status == CUBLAS_STATUS_SUCCESS );
}

// --------------------
extern "C"
void magma_copymatrix(
    magma_int_t m, magma_int_t n, size_t elemSize,
    void const* dA_src, magma_int_t lda,
    void*       dB_dst, magma_int_t ldb )
{
    cudaError_t status;
    status = cudaMemcpy2D(
        dB_dst, ldb*elemSize,
        dA_src, lda*elemSize,
        m*elemSize, n, cudaMemcpyDeviceToDevice );
    assert( status == cudaSuccess );
}

// --------------------
extern "C"
void magma_copymatrix_async(
    magma_int_t m, magma_int_t n, size_t elemSize,
    void const* dA_src, magma_int_t lda,
    void*       dB_dst, magma_int_t ldb,
    cudaStream_t stream )
{
    cudaError_t status;
    status = cudaMemcpy2DAsync(
        dB_dst, ldb*elemSize,
        dA_src, lda*elemSize,
        m*elemSize, n, cudaMemcpyDeviceToDevice, stream );
    assert( status == cudaSuccess );
}

#endif // HAVE_CUBLAS
