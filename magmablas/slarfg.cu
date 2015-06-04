/*
    -- MAGMA (version 1.6.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date November 2014

       @generated from zlarfg.cu normal z -> s, Sat Nov 15 19:53:57 2014
       
       @author Mark Gates
*/
#include "common_magma.h"
#include "magma_templates.h"

#define REAL

// 512 is maximum number of threads for CUDA capability 1.x
#define NB 512


// ----------------------------------------
// CUDA kernel for magma_slarfg.
// Uses one block of NB (currently 512) threads.
// Each thread sums dx[ tx + k*NB ]^2 for k = 0, 1, ...,
// then does parallel sum reduction to get norm-squared.
// 
// Currently setup to use NB threads, no matter how small dx is.
// This was slightly faster (5%) than passing n to magma_sum_reduce.
// To use number of threads = min( NB, max( 1, n-1 )), pass n as
// argument to magma_sum_reduce, rather than as template parameter.
__global__ void
slarfg_kernel(
    int n,
    float* dalpha, float* dx, int incx,
    float* dtau )
{
    const int tx = threadIdx.x;
    __shared__ float swork[ NB ];
    // TODO is it faster for each thread to have its own scale (register)?
    // if so, communicate it via swork[0]
    __shared__ float sscale;
    __shared__ float sscale2;
    float tmp;
    
    // find max of [dalpha, dx], to use as scaling to avoid unnecesary under- and overflow
    if ( tx == 0 ) {
        tmp = *dalpha;
        #ifdef COMPLEX
        swork[tx] = max( fabs(real(tmp)), fabs(imag(tmp)) );
        #else
        swork[tx] = fabs(tmp);
        #endif
    }
    else {
        swork[tx] = 0;
    }
    for( int j = tx; j < n-1; j += NB ) {
        tmp = dx[j*incx];
        #ifdef COMPLEX
        swork[tx] = max( swork[tx], max( fabs(real(tmp)), fabs(imag(tmp)) ));
        #else
        swork[tx] = max( swork[tx], fabs(tmp) );
        #endif
    }
    magma_max_reduce< NB >( tx, swork );
    if ( tx == 0 )
        sscale = swork[0];
    __syncthreads();
    
    // sum norm^2 of dx/sscale
    // dx has length n-1
    swork[tx] = 0;
    if ( sscale > 0 ) {
        for( int j = tx; j < n-1; j += NB ) {
            tmp = dx[j*incx] / sscale;
            swork[tx] += real(tmp)*real(tmp) + imag(tmp)*imag(tmp);
        }
        magma_sum_reduce< NB >( tx, swork );
        //magma_sum_reduce( blockDim.x, tx, swork );
    }
    
    if ( tx == 0 ) {
        float alpha = *dalpha;
        if ( swork[0] == 0 && imag(alpha) == 0 ) {
            // H = I
            *dtau = MAGMA_S_ZERO;
        }
        else {
            // beta = norm( [dalpha, dx] )
            float beta;
            tmp  = alpha / sscale;
            beta = sscale * sqrt( real(tmp)*real(tmp) + imag(tmp)*imag(tmp) + swork[0] );
            beta = -copysign( beta, real(alpha) );
            // todo: deal with badly scaled vectors (see lapack's larfg)
            *dtau   = MAGMA_S_MAKE( (beta - real(alpha)) / beta, -imag(alpha) / beta );
            *dalpha = MAGMA_S_MAKE( beta, 0 );
            sscale2 = 1 / (alpha - beta);
        }
    }
    
    // scale x (if norm was not 0)
    __syncthreads();
    if ( swork[0] != 0 ) {
        for( int j = tx; j < n-1; j += NB ) {
            dx[j*incx] *= sscale2;
        }
    }
}


/**
    Purpose
    -------
    SLARFG generates a real elementary reflector (Householder matrix)
    H of order n, such that

         H * ( alpha ) = ( beta ),   H**H * H = I.
             (   x   )   (   0  )

    where alpha and beta are scalars, with beta real and beta = ±norm([alpha, x]),
    and x is an (n-1)-element real vector. H is represented in the form

         H = I - tau * ( 1 ) * ( 1 v**H ),
                       ( v )

    where tau is a real scalar and v is a real (n-1)-element vector.
    Note that H is not symmetric.

    If the elements of x are all zero and dalpha is real, then tau = 0
    and H is taken to be the unit matrix.

    Otherwise  1 <= real(tau) <= 2  and  abs(tau-1) <= 1.

    Arguments
    ---------
    @param[in]
    n       INTEGER
            The order of the elementary reflector.

    @param[in,out]
    dalpha  REAL* on the GPU.
            On entry, pointer to the value alpha, i.e., the first entry of the vector.
            On exit, it is overwritten with the value beta.

    @param[in,out]
    dx      REAL array, dimension (1+(N-2)*abs(INCX)), on the GPU
            On entry, the (n-1)-element vector x.
            On exit, it is overwritten with the vector v.

    @param[in]
    incx    INTEGER
            The increment between elements of X. INCX > 0.

    @param[out]
    dtau    REAL* on the GPU.
            Pointer to the value tau.

    @param[in]
    queue   magma_queue_t
            Queue to execute in.

    @ingroup magma_saux1
    ********************************************************************/
extern "C"
void magmablas_slarfg_q(
    magma_int_t n,
    magmaFloat_ptr dalpha,
    magmaFloat_ptr dx, magma_int_t incx,
    magmaFloat_ptr dtau,
    magma_queue_t queue )
{
    dim3 blocks( 1 );
    dim3 threads( NB );
    //dim3 threads( min( NB, max( n-1, 1 )));
    slarfg_kernel<<< blocks, threads, 0, queue >>>( n, dalpha, dx, incx, dtau );
}


/**
    @see magmablas_slarfg_q
    @ingroup magma_saux1
    ********************************************************************/
extern "C"
void magmablas_slarfg(
    magma_int_t n,
    magmaFloat_ptr dalpha,
    magmaFloat_ptr dx, magma_int_t incx,
    magmaFloat_ptr dtau )
{
    dim3 blocks( 1 );
    dim3 threads( NB );
    //dim3 threads( min( NB, max( n-1, 1 )));
    slarfg_kernel<<< blocks, threads >>>( n, dalpha, dx, incx, dtau );
}
