/*
    -- MAGMA (version 1.3.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2012

       @generated d Wed Nov 14 22:53:47 2012

*/
#include "common_magma.h"

// 512 is maximum number of threads for CUDA capability 1.x
#define BLOCK_SIZE 512

#define PRECISION_d

__global__
void magma_dlarfg_gpu_kernel( int n, double* dx0, double* dx, 
                              double *dtau, double *dxnorm )
{
    const int i = threadIdx.x;
    const int j = i + BLOCK_SIZE * blockIdx.x;
    __shared__ double scale;
    __shared__ double xnorm;    
  
    double dxi;

    if ( j < n-1)
        dxi = dx[j];
  
    if ( i == 0 ) {
        xnorm = *dxnorm;
        if ( xnorm == 0 ) {
            *dtau = MAGMA_D_ZERO;
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
            double alpha = *dx0;
            double alphar =  MAGMA_D_REAL(alpha), alphai = MAGMA_D_IMAG(alpha);
            double beta  = xnorm; // sqrt( alphar*alphar + alphai*alphai + xnorm*xnorm );
            beta  = -copysign( beta, alphar );

            // todo: deal with badly scaled vectors (see lapack's larfg)
            *dtau = MAGMA_D_MAKE((beta - alphar)/beta, -alphai/beta);
            *dx0  = MAGMA_D_MAKE(beta, 0.);
            
            alpha = MAGMA_D_MAKE( MAGMA_D_REAL(alpha) - beta, MAGMA_D_IMAG(alpha));
            scale = MAGMA_D_DIV( MAGMA_D_ONE, alpha);
#endif
        }
    }

    // scale x
    __syncthreads();
    if ( xnorm != 0 && j < n-1)
        dx[j] = MAGMA_D_MUL(dxi, scale);
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
magma_dlarfg_gpu(int n, double *dx0, double *dx, 
                 double *dtau, double *dxnorm)
{
    dim3 blocks((n+BLOCK_SIZE-1) / BLOCK_SIZE);
    dim3 threads( BLOCK_SIZE );

    magma_dlarfg_gpu_kernel<<< blocks, threads >>>( n, dx0, dx, dtau, dxnorm );
}
