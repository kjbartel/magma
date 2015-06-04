/*
    -- MAGMA (version 1.3.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2012

       @generated c Wed Nov 14 22:53:47 2012

*/
#include "common_magma.h"

// 512 is maximum number of threads for CUDA capability 1.x
#define BLOCK_SIZE 512

#define PRECISION_c

__global__
void magma_clarfg_gpu_kernel( int n, cuFloatComplex* dx0, cuFloatComplex* dx, 
                              cuFloatComplex *dtau, float *dxnorm )
{
    const int i = threadIdx.x;
    const int j = i + BLOCK_SIZE * blockIdx.x;
    __shared__ cuFloatComplex scale;
    __shared__ float xnorm;    
  
    cuFloatComplex dxi;

    if ( j < n-1)
        dxi = dx[j];
  
    if ( i == 0 ) {
        xnorm = *dxnorm;
        if ( xnorm == 0 ) {
            *dtau = MAGMA_C_ZERO;
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
            cuFloatComplex alpha = *dx0;
            float alphar =  MAGMA_C_REAL(alpha), alphai = MAGMA_C_IMAG(alpha);
            float beta  = xnorm; // sqrt( alphar*alphar + alphai*alphai + xnorm*xnorm );
            beta  = -copysign( beta, alphar );

            // todo: deal with badly scaled vectors (see lapack's larfg)
            *dtau = MAGMA_C_MAKE((beta - alphar)/beta, -alphai/beta);
            *dx0  = MAGMA_C_MAKE(beta, 0.);
            
            alpha = MAGMA_C_MAKE( MAGMA_C_REAL(alpha) - beta, MAGMA_C_IMAG(alpha));
            scale = MAGMA_C_DIV( MAGMA_C_ONE, alpha);
#endif
        }
    }

    // scale x
    __syncthreads();
    if ( xnorm != 0 && j < n-1)
        dx[j] = MAGMA_C_MUL(dxi, scale);
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
magma_clarfg_gpu(int n, cuFloatComplex *dx0, cuFloatComplex *dx, 
                 cuFloatComplex *dtau, float *dxnorm)
{
    dim3 blocks((n+BLOCK_SIZE-1) / BLOCK_SIZE);
    dim3 threads( BLOCK_SIZE );

    magma_clarfg_gpu_kernel<<< blocks, threads >>>( n, dx0, dx, dtau, dxnorm );
}
