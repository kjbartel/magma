/*
    -- MAGMA (version 1.6.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date November 2014
       
       @author Mark Gates

       @generated from zswap.cu normal z -> s, Sat Nov 15 19:53:59 2014

*/
#include "common_magma.h"

#define NB 64


/* Vector is divided into ceil(n/nb) blocks.
   Each thread swaps one element, x[tid] <---> y[tid].
*/
__global__ void sswap_kernel(
    int n,
    float *x, int incx,
    float *y, int incy )
{
    float tmp;
    int ind = threadIdx.x + blockDim.x*blockIdx.x;
    if ( ind < n ) {
        x += ind*incx;
        y += ind*incy;
        tmp = *x;
        *x  = *y;
        *y  = tmp;
    }
}


/**
    Purpose:
    =============
    Swap vector x and y; \f$ x <-> y \f$.

    @param[in]
    n       Number of elements in vector x and y. n >= 0.

    @param[in,out]
    dx      REAL array on GPU device.
            The n element vector x of dimension (1 + (n-1)*incx).

    @param[in]
    incx    Stride between consecutive elements of dx. incx != 0.

    @param[in,out]
    dy      REAL array on GPU device.
            The n element vector y of dimension (1 + (n-1)*incy).

    @param[in]
    incy    Stride between consecutive elements of dy. incy != 0.

    @ingroup magma_sblas1
    ********************************************************************/
extern "C" void 
magmablas_sswap_q(
    magma_int_t n,
    magmaFloat_ptr dx, magma_int_t incx, 
    magmaFloat_ptr dy, magma_int_t incy,
    magma_queue_t queue )
{
    dim3 blocks( (n+NB-1) / NB );
    dim3 threads( NB );
    sswap_kernel<<< blocks, threads, 0, queue >>>( n, dx, incx, dy, incy );
}


/**
    @see magmablas_sswap_q
    @ingroup magma_sblas1
    ********************************************************************/
extern "C" void 
magmablas_sswap(
    magma_int_t n,
    magmaFloat_ptr dx, magma_int_t incx, 
    magmaFloat_ptr dy, magma_int_t incy)
{
    magmablas_sswap_q( n, dx, incx, dy, incy, magma_stream );
}
