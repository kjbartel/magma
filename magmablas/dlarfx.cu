/*
    -- MAGMA (version 1.5.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date September 2014

       @generated from zlarfx.cu normal z -> d, Wed Sep 17 15:08:23 2014

*/
#include "common_magma.h"
#include "commonblas_d.h"
#include "magma_templates.h"

// 512 is maximum number of threads for CUDA capability 1.x
#define BLOCK_SIZE 512

#define BLOCK_SIZEx  32
#define BLOCK_SIZEy  16


//==============================================================================

__global__
void magma_dlarfx_kernel( int m, double *v, double *tau,
                         double *c, int ldc, double *xnorm,
                         double *T, int it )
{
    if ( !MAGMA_D_EQUAL(*tau, MAGMA_D_ZERO) ) {
        const int tx = threadIdx.x;
        //double *dc = c + (blockIdx.x-it-1) * ldc;
        double *dc = c + (blockIdx.x) * ldc;

        __shared__ double sum[ BLOCK_SIZE ];
        double lsum;

        /* NOTE HERE C is the C at position C(i, 0) 
         * if blockIdx.x<it it performs the V(i:n,i)' * V(i:n,1:i-1)' used for computing T
         * if blockIdx.x>it it perform  w := v' * C  */
        lsum = MAGMA_D_ZERO;
        for( int j = tx; j < m; j += BLOCK_SIZE ){
            if (j==0){
               lsum += MAGMA_D_MUL( MAGMA_D_ONE, dc[j] );
               v[j] = MAGMA_D_ONE;
            }
            else
               lsum += MAGMA_D_MUL( MAGMA_D_CNJG( v[j] ), dc[j] );
        }
        sum[tx] = lsum;
        magma_sum_reduce< BLOCK_SIZE >( tx, sum );

        /*  C := C - v * w  */
        __syncthreads();
        double z__1 = - MAGMA_D_CNJG(*tau) * sum[0];
        if (blockIdx.x>it){
           for( int j = m-tx-1; j>=0 ; j -= BLOCK_SIZE )
                 dc[j] += z__1 * v[j];
           __syncthreads();

           /* Adjust the rest of the column norms */
           /*
           if (tx==0){
             double temp = MAGMA_D_ABS( dc[0] ) / xnorm[blockIdx.x-it-1];
             temp = (temp + 1.) * (1. - temp);
             xnorm[blockIdx.x-it-1] = xnorm[blockIdx.x-it-1] * sqrt(temp); 
           }
           */
        }
        else
        {
           if (blockIdx.x==it)
              *(T+it) = *tau;
           else
              *(T+blockIdx.x) = MAGMA_D_CNJG(z__1);
        }
    }
    else if (blockIdx.x<=it)// in case tau is zero put the corresponding column of T to zero
    {
        *(T+blockIdx.x) = MAGMA_D_ZERO;
    }

}

//==============================================================================
extern "C"
__global__
void magma_dtrmv_kernel(const double *T, int ldt, double *t)
{
   const int tx = threadIdx.x;
   T += tx;

   __shared__ double tlocal[ BLOCK_SIZE ];
   double res = MAGMA_D_MAKE(0., 0.);

   tlocal[tx] = t[tx];
   __syncthreads();

   #pragma unroll
   for(int j=0; j<blockDim.x; j++)
      res +=  T[j*ldt]*tlocal[j];

   t[tx] = res;
}

extern "C"
__global__
void magma_dtrmv_kernel2(const double *T, int ldt, double *t, 
                         double *y, double *tau)
{
   const int tx = threadIdx.x;
   T += blockIdx.x;

   __shared__ double sum[ 128 ];

   sum[tx] = T[tx*ldt]*t[tx];
   magma_sum_reduce_n(blockDim.x, tx, sum);

   __syncthreads();

   if (tx==0){
      y[blockIdx.x] = sum[0];
      if (blockIdx.x==0)
         y[gridDim.x] = tau[0];
   }
}

//==============================================================================
extern "C"
__global__
void magma_dtrmv_tkernel(double *T, int ldt, double *t, double *y)
{
   const int tx = threadIdx.x;
   T += blockIdx.x*ldt;

   __shared__ double sum[ 128 ];

   sum[tx] = MAGMA_D_CNJG(T[tx])*t[tx];
   magma_sum_reduce_n(blockDim.x, tx, sum);

   __syncthreads();

   if (tx==0)
      y[blockIdx.x] = sum[0];
}

//==============================================================================

/*
    Apply a real elementary reflector H to a real M-by-N
    matrix C from the left. H is represented in the form
          H = I - tau * v * v'
    where tau is a real scalar and v is a real vector.
    If tau = 0, then H is taken to be the unit matrix.

    To apply H' (the conjugate transpose of H), supply conjg(tau) 
    instead tau.

    The norms of v(:, 1:n) are given as input in xnorm(1:n). On exit, the norms
    are adjusted to hold the norms of v(2:m,2:n). This is a difference with the 
    LAPACK's dlarf routine. 
 */
extern "C" void
magma_dlarfx_gpu(magma_int_t m, magma_int_t n, double *v, double *tau,
                double *c, magma_int_t ldc, double *xnorm, 
                double *T, magma_int_t i, double *work )
{
    magma_int_t N = n + i + 1;

    if (i==0)
        magma_dlarfx_kernel<<< N, BLOCK_SIZE, 0, magma_stream >>>( m, v, tau, c, ldc, xnorm, T+i*N, i);
    else
        magma_dlarfx_kernel<<< N, BLOCK_SIZE, 0, magma_stream >>>( m, v, tau, c, ldc, xnorm, work, i);

    if (i > 0){
        //magma_dtrmv_kernel<<< 1, i, 0, magma_stream >>>( T, N, T+i*N);
        magma_dtrmv_kernel2<<< i, i, 0, magma_stream  >>>( T, N, work, T+i*N, tau);
    }
}

//==============================================================================
