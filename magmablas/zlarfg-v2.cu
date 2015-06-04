/*
    -- MAGMA (version 1.3.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2012

       @precisions normal z -> s d c

*/
#include "common_magma.h"

// 512 is maximum number of threads for CUDA capability 1.x
#define BLOCK_SIZE 512

#define PRECISION_z

__global__
void magma_zlarfg_gpu_kernel( int n, cuDoubleComplex* dx0, cuDoubleComplex* dx, 
                              cuDoubleComplex *dtau, double *dxnorm )
{
    const int i = threadIdx.x;
    const int j = i + BLOCK_SIZE * blockIdx.x;
    __shared__ cuDoubleComplex scale;
    __shared__ double xnorm;    
  
    cuDoubleComplex dxi;

    if ( j < n-1)
        dxi = dx[j];
  
    if ( i == 0 ) {
        xnorm = *dxnorm;
        if ( xnorm == 0 ) {
            *dtau = MAGMA_Z_ZERO;
        }
        else {

#if (defined(PRECISION_s) || defined(PRECISION_d))
            double alpha = *dx0;
            double beta  = xnorm; // sqrt( alpha*alpha + xnorm*xnorm );
            beta  = -copysign( beta, alpha );

            // todo: deal with badly scaled vectors (see lapack's larfg)
            *dtau = (beta - alpha) / beta;
            *dx0  = beta;

            scale = 1 / (alpha - beta);
#else
            cuDoubleComplex alpha = *dx0;
            double alphar =  MAGMA_Z_REAL(alpha), alphai = MAGMA_Z_IMAG(alpha);
            double beta  = xnorm; // sqrt( alphar*alphar + alphai*alphai + xnorm*xnorm );
            beta  = -copysign( beta, alphar );

            // todo: deal with badly scaled vectors (see lapack's larfg)
            *dtau = MAGMA_Z_MAKE((beta - alphar)/beta, -alphai/beta);
            *dx0  = MAGMA_Z_MAKE(beta, 0.);
            
            alpha = MAGMA_Z_MAKE( MAGMA_Z_REAL(alpha) - beta, MAGMA_Z_IMAG(alpha));
            scale = MAGMA_Z_DIV( MAGMA_Z_ONE, alpha);
#endif
        }
    }

    // scale x
    __syncthreads();
    if ( xnorm != 0 && j < n-1)
        dx[j] = MAGMA_Z_MUL(dxi, scale);
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
magma_zlarfg_gpu(int n, cuDoubleComplex *dx0, cuDoubleComplex *dx, 
                 cuDoubleComplex *dtau, double *dxnorm)
{
    dim3 blocks((n+BLOCK_SIZE-1) / BLOCK_SIZE);
    dim3 threads( BLOCK_SIZE );

    magma_zlarfg_gpu_kernel<<< blocks, threads >>>( n, dx0, dx, dtau, dxnorm );
}
