/*
 *   -- MAGMA (version 1.2.1) --
 *      Univ. of Tennessee, Knoxville
 *      Univ. of California, Berkeley
 *      Univ. of Colorado, Denver
 *      June 2012
 *
 * @author Mark Gates
 * @generated d Thu Jun 28 12:31:08 2012
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
void magma_dsetvector(
    magma_int_t n,
    double const* hx_src, magma_int_t incx,
    double*       dy_dst, magma_int_t incy )
{
    cublasStatus_t status;
    status = cublasSetVector(
        n, sizeof(double),
        hx_src, incx,
        dy_dst, incy );
    check_error( status );
}

// --------------------
extern "C"
void magma_dgetvector(
    magma_int_t n,
    double const* dx_src, magma_int_t incx,
    double*       hy_dst, magma_int_t incy )
{
    cublasStatus_t status;
    status = cublasGetVector(
        n, sizeof(double),
        dx_src, incx,
        hy_dst, incy );
    check_error( status );
}

// --------------------
extern "C"
void magma_dsetvector_async(
    magma_int_t n,
    double const* hx_src, magma_int_t incx,
    double*       dy_dst, magma_int_t incy,
    cudaStream_t stream )
{
    cublasStatus_t status;
    status = cublasSetVectorAsync(
        n, sizeof(double),
        hx_src, incx,
        dy_dst, incy, stream );
    check_error( status );
}

// --------------------
extern "C"
void magma_dgetvector_async(
    magma_int_t n,
    double const* dx_src, magma_int_t incx,
    double*       hy_dst, magma_int_t incy,
    cudaStream_t stream )
{
    cublasStatus_t status;
    status = cublasGetVectorAsync(
        n, sizeof(double),
        dx_src, incx,
        hy_dst, incy, stream );
    check_error( status );
}


// ========================================
// copying sub-matrices (contiguous columns)
extern "C"
void magma_dsetmatrix(
    magma_int_t m, magma_int_t n,
    double const* hA_src, magma_int_t lda,
    double*       dB_dst, magma_int_t ldb )
{
    cublasStatus_t status;
    status = cublasSetMatrix(
        m, n, sizeof(double),
        hA_src, lda,
        dB_dst, ldb );
    check_error( status );
}

// --------------------
extern "C"
void magma_dgetmatrix(
    magma_int_t m, magma_int_t n,
    double const* dA_src, magma_int_t lda,
    double*       hB_dst, magma_int_t ldb )
{
    cublasStatus_t status;
    status = cublasGetMatrix(
        m, n, sizeof(double),
        dA_src, lda,
        hB_dst, ldb );
    check_error( status );
}

// --------------------
extern "C"
void magma_dsetmatrix_async(
    magma_int_t m, magma_int_t n,
    double const* hA_src, magma_int_t lda,
    double*       dB_dst, magma_int_t ldb,
    cudaStream_t stream )
{
    cublasStatus_t status;
    status = cublasSetMatrixAsync(
        m, n, sizeof(double),
        hA_src, lda,
        dB_dst, ldb, stream );
    check_error( status );
}

// --------------------
extern "C"
void magma_dgetmatrix_async(
    magma_int_t m, magma_int_t n,
    double const* dA_src, magma_int_t lda,
    double*       hB_dst, magma_int_t ldb,
    cudaStream_t stream )
{
    cublasStatus_t status;
    status = cublasGetMatrixAsync(
        m, n, sizeof(double),
        dA_src, lda,
        hB_dst, ldb, stream );
    check_error( status );
}

// --------------------
extern "C"
void magma_dcopymatrix(
    magma_int_t m, magma_int_t n,
    double const* dA_src, magma_int_t lda,
    double*       dB_dst, magma_int_t ldb )
{
    cudaError_t status;
    status = cudaMemcpy2D(
        dB_dst, ldb*sizeof(double),
        dA_src, lda*sizeof(double),
        m*sizeof(double), n, cudaMemcpyDeviceToDevice );
    assert( status == cudaSuccess );
}

// --------------------
extern "C"
void magma_dcopymatrix_async(
    magma_int_t m, magma_int_t n,
    double const* dA_src, magma_int_t lda,
    double*       dB_dst, magma_int_t ldb,
    cudaStream_t stream )
{
    cudaError_t status;
    status = cudaMemcpy2DAsync(
        dB_dst, ldb*sizeof(double),
        dA_src, lda*sizeof(double),
        m*sizeof(double), n, cudaMemcpyDeviceToDevice, stream );
    assert( status == cudaSuccess );
}

#endif // HAVE_CUBLAS
