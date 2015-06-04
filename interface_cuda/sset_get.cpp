/*
 *   -- MAGMA (version 1.2.0) --
 *      Univ. of Tennessee, Knoxville
 *      Univ. of California, Berkeley
 *      Univ. of Colorado, Denver
 *      April 2012
 *
 * @author Mark Gates
 * @generated s Tue May 15 18:17:56 2012
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
void magma_ssetvector(
    magma_int_t n,
    float const* hx_src, magma_int_t incx,
    float*       dy_dst, magma_int_t incy )
{
    cublasStatus_t status;
    status = cublasSetVector(
        n, sizeof(float),
        hx_src, incx,
        dy_dst, incy );
    check_error( status );
}

// --------------------
extern "C"
void magma_sgetvector(
    magma_int_t n,
    float const* dx_src, magma_int_t incx,
    float*       hy_dst, magma_int_t incy )
{
    cublasStatus_t status;
    status = cublasGetVector(
        n, sizeof(float),
        dx_src, incx,
        hy_dst, incy );
    check_error( status );
}

// --------------------
extern "C"
void magma_ssetvector_async(
    magma_int_t n,
    float const* hx_src, magma_int_t incx,
    float*       dy_dst, magma_int_t incy,
    cudaStream_t stream )
{
    cublasStatus_t status;
    status = cublasSetVectorAsync(
        n, sizeof(float),
        hx_src, incx,
        dy_dst, incy, stream );
    check_error( status );
}

// --------------------
extern "C"
void magma_sgetvector_async(
    magma_int_t n,
    float const* dx_src, magma_int_t incx,
    float*       hy_dst, magma_int_t incy,
    cudaStream_t stream )
{
    cublasStatus_t status;
    status = cublasGetVectorAsync(
        n, sizeof(float),
        dx_src, incx,
        hy_dst, incy, stream );
    check_error( status );
}


// ========================================
// copying sub-matrices (contiguous columns)
extern "C"
void magma_ssetmatrix(
    magma_int_t m, magma_int_t n,
    float const* hA_src, magma_int_t lda,
    float*       dB_dst, magma_int_t ldb )
{
    cublasStatus_t status;
    status = cublasSetMatrix(
        m, n, sizeof(float),
        hA_src, lda,
        dB_dst, ldb );
    check_error( status );
}

// --------------------
extern "C"
void magma_sgetmatrix(
    magma_int_t m, magma_int_t n,
    float const* dA_src, magma_int_t lda,
    float*       hB_dst, magma_int_t ldb )
{
    cublasStatus_t status;
    status = cublasGetMatrix(
        m, n, sizeof(float),
        dA_src, lda,
        hB_dst, ldb );
    check_error( status );
}

// --------------------
extern "C"
void magma_ssetmatrix_async(
    magma_int_t m, magma_int_t n,
    float const* hA_src, magma_int_t lda,
    float*       dB_dst, magma_int_t ldb,
    cudaStream_t stream )
{
    cublasStatus_t status;
    status = cublasSetMatrixAsync(
        m, n, sizeof(float),
        hA_src, lda,
        dB_dst, ldb, stream );
    check_error( status );
}

// --------------------
extern "C"
void magma_sgetmatrix_async(
    magma_int_t m, magma_int_t n,
    float const* dA_src, magma_int_t lda,
    float*       hB_dst, magma_int_t ldb,
    cudaStream_t stream )
{
    cublasStatus_t status;
    status = cublasGetMatrixAsync(
        m, n, sizeof(float),
        dA_src, lda,
        hB_dst, ldb, stream );
    check_error( status );
}

// --------------------
extern "C"
void magma_scopymatrix(
    magma_int_t m, magma_int_t n,
    float const* dA_src, magma_int_t lda,
    float*       dB_dst, magma_int_t ldb )
{
    cudaError_t status;
    status = cudaMemcpy2D(
        dB_dst, ldb*sizeof(float),
        dA_src, lda*sizeof(float),
        m*sizeof(float), n, cudaMemcpyDeviceToDevice );
    assert( status == cudaSuccess );
}

// --------------------
extern "C"
void magma_scopymatrix_async(
    magma_int_t m, magma_int_t n,
    float const* dA_src, magma_int_t lda,
    float*       dB_dst, magma_int_t ldb,
    cudaStream_t stream )
{
    cudaError_t status;
    status = cudaMemcpy2DAsync(
        dB_dst, ldb*sizeof(float),
        dA_src, lda*sizeof(float),
        m*sizeof(float), n, cudaMemcpyDeviceToDevice, stream );
    assert( status == cudaSuccess );
}

#endif // HAVE_CUBLAS
