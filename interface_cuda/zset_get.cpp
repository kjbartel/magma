/*
 *   -- MAGMA (version 1.2.1) --
 *      Univ. of Tennessee, Knoxville
 *      Univ. of California, Berkeley
 *      Univ. of Colorado, Denver
 *      June 2012
 *
 * @author Mark Gates
 * @precisions normal z -> s d c
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
void magma_zsetvector(
    magma_int_t n,
    cuDoubleComplex const* hx_src, magma_int_t incx,
    cuDoubleComplex*       dy_dst, magma_int_t incy )
{
    cublasStatus_t status;
    status = cublasSetVector(
        n, sizeof(cuDoubleComplex),
        hx_src, incx,
        dy_dst, incy );
    check_error( status );
}

// --------------------
extern "C"
void magma_zgetvector(
    magma_int_t n,
    cuDoubleComplex const* dx_src, magma_int_t incx,
    cuDoubleComplex*       hy_dst, magma_int_t incy )
{
    cublasStatus_t status;
    status = cublasGetVector(
        n, sizeof(cuDoubleComplex),
        dx_src, incx,
        hy_dst, incy );
    check_error( status );
}

// --------------------
extern "C"
void magma_zsetvector_async(
    magma_int_t n,
    cuDoubleComplex const* hx_src, magma_int_t incx,
    cuDoubleComplex*       dy_dst, magma_int_t incy,
    cudaStream_t stream )
{
    cublasStatus_t status;
    status = cublasSetVectorAsync(
        n, sizeof(cuDoubleComplex),
        hx_src, incx,
        dy_dst, incy, stream );
    check_error( status );
}

// --------------------
extern "C"
void magma_zgetvector_async(
    magma_int_t n,
    cuDoubleComplex const* dx_src, magma_int_t incx,
    cuDoubleComplex*       hy_dst, magma_int_t incy,
    cudaStream_t stream )
{
    cublasStatus_t status;
    status = cublasGetVectorAsync(
        n, sizeof(cuDoubleComplex),
        dx_src, incx,
        hy_dst, incy, stream );
    check_error( status );
}


// ========================================
// copying sub-matrices (contiguous columns)
extern "C"
void magma_zsetmatrix(
    magma_int_t m, magma_int_t n,
    cuDoubleComplex const* hA_src, magma_int_t lda,
    cuDoubleComplex*       dB_dst, magma_int_t ldb )
{
    cublasStatus_t status;
    status = cublasSetMatrix(
        m, n, sizeof(cuDoubleComplex),
        hA_src, lda,
        dB_dst, ldb );
    check_error( status );
}

// --------------------
extern "C"
void magma_zgetmatrix(
    magma_int_t m, magma_int_t n,
    cuDoubleComplex const* dA_src, magma_int_t lda,
    cuDoubleComplex*       hB_dst, magma_int_t ldb )
{
    cublasStatus_t status;
    status = cublasGetMatrix(
        m, n, sizeof(cuDoubleComplex),
        dA_src, lda,
        hB_dst, ldb );
    check_error( status );
}

// --------------------
extern "C"
void magma_zsetmatrix_async(
    magma_int_t m, magma_int_t n,
    cuDoubleComplex const* hA_src, magma_int_t lda,
    cuDoubleComplex*       dB_dst, magma_int_t ldb,
    cudaStream_t stream )
{
    cublasStatus_t status;
    status = cublasSetMatrixAsync(
        m, n, sizeof(cuDoubleComplex),
        hA_src, lda,
        dB_dst, ldb, stream );
    check_error( status );
}

// --------------------
extern "C"
void magma_zgetmatrix_async(
    magma_int_t m, magma_int_t n,
    cuDoubleComplex const* dA_src, magma_int_t lda,
    cuDoubleComplex*       hB_dst, magma_int_t ldb,
    cudaStream_t stream )
{
    cublasStatus_t status;
    status = cublasGetMatrixAsync(
        m, n, sizeof(cuDoubleComplex),
        dA_src, lda,
        hB_dst, ldb, stream );
    check_error( status );
}

// --------------------
extern "C"
void magma_zcopymatrix(
    magma_int_t m, magma_int_t n,
    cuDoubleComplex const* dA_src, magma_int_t lda,
    cuDoubleComplex*       dB_dst, magma_int_t ldb )
{
    cudaError_t status;
    status = cudaMemcpy2D(
        dB_dst, ldb*sizeof(cuDoubleComplex),
        dA_src, lda*sizeof(cuDoubleComplex),
        m*sizeof(cuDoubleComplex), n, cudaMemcpyDeviceToDevice );
    assert( status == cudaSuccess );
}

// --------------------
extern "C"
void magma_zcopymatrix_async(
    magma_int_t m, magma_int_t n,
    cuDoubleComplex const* dA_src, magma_int_t lda,
    cuDoubleComplex*       dB_dst, magma_int_t ldb,
    cudaStream_t stream )
{
    cudaError_t status;
    status = cudaMemcpy2DAsync(
        dB_dst, ldb*sizeof(cuDoubleComplex),
        dA_src, lda*sizeof(cuDoubleComplex),
        m*sizeof(cuDoubleComplex), n, cudaMemcpyDeviceToDevice, stream );
    assert( status == cudaSuccess );
}

#endif // HAVE_CUBLAS
