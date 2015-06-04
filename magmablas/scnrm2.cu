/*
    -- MAGMA (version 1.3.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2012

       @generated c Wed Nov 14 22:53:53 2012

*/
#include "common_magma.h"

// 512 is maximum number of threads for CUDA capability 1.x
#define BLOCK_SIZE 512

#define PRECISION_c

// ----------------------------------------
// Does sum reduction of array x, leaving total in x[0].
// Contents of x are destroyed in the process.
// With k threads, can reduce array up to 2*k in size.
// Assumes number of threads <= 1024 (which is max number of threads up to CUDA capability 3.0)
// Having n as template parameter allows compiler to evaluate some conditions at compile time.
template< int n >
__device__ void sum_reduce( /*int n,*/ int i, float* x )
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


__global__ void
magmablas_scnrm2_kernel( int m, cuFloatComplex *da, int ldda, float *dxnorm )
{
    const int i = threadIdx.x;
    cuFloatComplex *dx = da + blockIdx.x * ldda;

    __shared__ float sum[ BLOCK_SIZE ];
    float re;

    // get norm of dx
    sum[i] = 0;
    for( int j = i; j < m; j += BLOCK_SIZE ) {

#if (defined(PRECISION_s) || defined(PRECISION_d))
        re = dx[j];
        sum[i] += re*re;
#else
        re = MAGMA_C_REAL( dx[j] );
        float im = MAGMA_C_IMAG( dx[j] );
        sum[i] += re*re + im*im;
#endif

    }
    sum_reduce< BLOCK_SIZE >( i, sum );
    
    if (i==0)
       dxnorm[blockIdx.x] = sqrt(sum[0]);
}

/*
   Compute the scnrm2 of da, da+ldda, ..., da +(num-1)*ldda where the vectors are
   of size m. The resulting norms are written in the dxnorm array. 
*/
extern "C" void
magmablas_scnrm2(int m, int num, cuFloatComplex *da, magma_int_t ldda, 
                 float *dxnorm) 
{
    dim3  blocks( num );
    dim3 threads( BLOCK_SIZE );
    
    magmablas_scnrm2_kernel<<< blocks, threads >>>( m, da, ldda, dxnorm );
}

