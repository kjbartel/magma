/*
    -- MAGMA (version 1.3.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2012

       @generated s Wed Nov 14 22:53:47 2012

*/
#include "common_magma.h"

// 512 is maximum number of threads for CUDA capability 1.x
#define BLOCK_SIZE 512

#define PRECISION_s

__global__
void magma_slarfg_gpu_kernel( int n, float* dx0, float* dx, 
                              float *dtau, float *dxnorm )
{
    const int i = threadIdx.x;
    const int j = i + BLOCK_SIZE * blockIdx.x;
    __shared__ float scale;
    __shared__ float xnorm;    
  
    float dxi;

    if ( j < n-1)
        dxi = dx[j];
  
    if ( i == 0 ) {
        xnorm = *dxnorm;
        if ( xnorm == 0 ) {
            *dtau = MAGMA_S_ZERO;
        }
        else {

#if (defined(PRECISION_s) || defined(PRECISION_d))
            float alpha = *dx0;
            float beta  = xnorm; // sqrt( alpha*alpha + xnorm*xnorm );
            beta  = -copysign( beta, alpha );

            // todo: deal with badly scaled vectors (see lapack's larfg)
            *dtau = (beta - alpha) / beta;
            *dx0  = beta;

            scale = 1 / (alpha - beta);
#else
            float alpha = *dx0;
            float alphar =  MAGMA_S_REAL(alpha), alphai = MAGMA_S_IMAG(alpha);
            float beta  = xnorm; // sqrt( alphar*alphar + alphai*alphai + xnorm*xnorm );
            beta  = -copysign( beta, alphar );

            // todo: deal with badly scaled vectors (see lapack's larfg)
            *dtau = MAGMA_S_MAKE((beta - alphar)/beta, -alphai/beta);
            *dx0  = MAGMA_S_MAKE(beta, 0.);
            
            alpha = MAGMA_S_MAKE( MAGMA_S_REAL(alpha) - beta, MAGMA_S_IMAG(alpha));
            scale = MAGMA_S_DIV( MAGMA_S_ONE, alpha);
#endif
        }
    }

    // scale x
    __syncthreads();
    if ( xnorm != 0 && j < n-1)
        dx[j] = MAGMA_S_MUL(dxi, scale);
}

/*
   Generates Householder elementary reflector H = I - tau v v^T to reduce
     H [ dx0 ] = [ beta ]
       [ dx  ]   [ 0    ]
   with beta = ±norm( [dx0, dx] ).
   Stores v over dx; first element of v is 1 and is not stored.
   Stores beta over dx0.
   Stores tau.  
*/
extern "C" void
magma_slarfg_gpu(int n, float *dx0, float *dx, 
                 float *dtau, float *dxnorm)
{
    dim3 blocks((n+BLOCK_SIZE-1) / BLOCK_SIZE);
    dim3 threads( BLOCK_SIZE );

    magma_slarfg_gpu_kernel<<< blocks, threads >>>( n, dx0, dx, dtau, dxnorm );
}
