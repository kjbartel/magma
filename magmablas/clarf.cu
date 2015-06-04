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

// ----------------------------------------
// Does sum reduction of array x, leaving total in x[0].
// Contents of x are destroyed in the process.
// With k threads, can reduce array up to 2*k in size.
// Assumes number of threads <= 1024 (which is max number of threads up to CUDA capability 3.0)
// Having n as template parameter allows compiler to evaluate some conditions at compile time.
template< int n >
__device__ void sum_reduce( /*int n,*/ int i, cuFloatComplex* x )
{
    __syncthreads();
    if ( n > 1024 ) { if ( i < 1024 && i + 1024 < n ) { x[i] += x[i+1024]; }  __syncthreads(); }
    if ( n >  512 ) { if ( i <  512 && i +  512 < n ) { x[i] += x[i+ 512]; }  __syncthreads(); }
    if ( n >  256 ) { if ( i <  256 && i +  256 < n ) { x[i] += x[i+ 256]; }  __syncthreads(); }
    if ( n >  128 ) { if ( i <  128 && i +  128 < n ) { x[i] += x[i+ 128]; }  __syncthreads(); }
    if ( n >   64 ) { if ( i <   64 && i +   64 < n ) { x[i] += x[i+  64]; }  __syncthreads(); }
    if ( n >   32 ) { if ( i <   32 && i +   32 < n ) { x[i] += x[i+  32]; }  __syncthreads(); }
    // probably don't need __syncthreads for < 16 threads
    // because of implicit warp level synchronization.
    if ( n >   16 ) { if ( i <   16 && i +   16 < n ) { x[i] += x[i+  16]; }  __syncthreads(); }
    if ( n >    8 ) { if ( i <    8 && i +    8 < n ) { x[i] += x[i+   8]; }  __syncthreads(); }
    if ( n >    4 ) { if ( i <    4 && i +    4 < n ) { x[i] += x[i+   4]; }  __syncthreads(); }
    if ( n >    2 ) { if ( i <    2 && i +    2 < n ) { x[i] += x[i+   2]; }  __syncthreads(); }
    if ( n >    1 ) { if ( i <    1 && i +    1 < n ) { x[i] += x[i+   1]; }  __syncthreads(); }
}
// end sum_reduce


__global__
void magma_clarf_kernel( int m, cuFloatComplex *v, cuFloatComplex *tau,
                         cuFloatComplex *c, int ldc, float *xnorm )
{
    if ( !MAGMA_C_EQUAL(*tau, MAGMA_C_ZERO) ) {
        const int i = threadIdx.x;
        cuFloatComplex *dc = c + blockIdx.x * ldc;

        __shared__ cuFloatComplex sum[ BLOCK_SIZE ];

        /*  w := v' * C  */
        sum[i] = MAGMA_C_ZERO;
        for( int j = i; j < m; j += BLOCK_SIZE ){
            if (j==0)
               sum[i] += MAGMA_C_MUL( MAGMA_C_ONE, dc[j] );
            else
               sum[i] += MAGMA_C_MUL( MAGMA_C_CNJG( v[j] ), dc[j] );
        }
        sum_reduce< BLOCK_SIZE >( i, sum );

        /*  C := C - v * w  */
        __syncthreads();
        cuFloatComplex z__1 = - MAGMA_C_CNJG(*tau) * sum[0];
        for( int j = m-i-1; j>=0 ; j -= BLOCK_SIZE ) {
             if (j==0)
                dc[j] += z__1;
             else
                dc[j] += z__1 * v[j];
        }
        __syncthreads();

        /* Adjust the rest of the column norms */
        if (i==0){
            float temp = MAGMA_C_ABS( dc[0] ) / xnorm[blockIdx.x];
            temp = (temp + 1.) * (1. - temp);
            xnorm[blockIdx.x] = xnorm[blockIdx.x] * sqrt(temp); 
        }
    }
}

/*
    Apply a complex elementary reflector H to a complex M-by-N
    matrix C from the left. H is represented in the form
          H = I - tau * v * v'
    where tau is a complex scalar and v is a complex vector.
    If tau = 0, then H is taken to be the unit matrix.

    To apply H' (the conjugate transpose of H), supply conjg(tau) 
    instead tau.
 */
extern "C" void
magma_clarf_gpu(int m, int n, cuFloatComplex *v, cuFloatComplex *tau,
                cuFloatComplex *c, int ldc, float *xnorm)
{
    dim3  blocks( n );
    dim3 threads( BLOCK_SIZE );

    magma_clarf_kernel<<< blocks, threads >>>( m, v, tau, c, ldc, xnorm);
}
