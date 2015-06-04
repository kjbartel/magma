/*
    -- MAGMA (version 1.5.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date September 2014

       @generated from dznrm2.cu normal z -> c, Wed Sep 17 15:08:23 2014

*/
#include "common_magma.h"
#include "commonblas_c.h"
#include "magma_templates.h"

// 512 is maximum number of threads for CUDA capability 1.x
#define BLOCK_SIZE  512
#define BLOCK_SIZEx  32
#define BLOCK_SIZEy  16

#define PRECISION_c


//==============================================================================

__global__ void
magmablas_scnrm2_kernel( int m, magmaFloatComplex *da, int ldda, float *dxnorm )
{
    const int tx = threadIdx.x;
    magmaFloatComplex *dx = da + blockIdx.x * ldda;

    __shared__ float sum[ BLOCK_SIZE ];
    float re, lsum;

    // get norm of dx
    lsum = 0;
    for( int j = tx; j < m; j += BLOCK_SIZE ) {
        #if (defined(PRECISION_s) || defined(PRECISION_d))
        re = dx[j];
        lsum += re*re;
        #else
        re = MAGMA_C_REAL( dx[j] );
        float im = MAGMA_C_IMAG( dx[j] );
        lsum += re*re + im*im;
        #endif
    }
    sum[tx] = lsum;
    magma_sum_reduce< BLOCK_SIZE >( tx, sum );
    
    if (tx==0)
        dxnorm[blockIdx.x] = sqrt(sum[0]);
}


//==============================================================================
__global__ void
magmablas_scnrm2_check_kernel( int m, magmaFloatComplex *da, int ldda, float *dxnorm, 
                               float *lsticc )
{
    const int tx = threadIdx.x;
    magmaFloatComplex *dx = da + blockIdx.x * ldda;

    __shared__ float sum[ BLOCK_SIZE ];
    float re, lsum;

    // get norm of dx only if lsticc[blockIdx+1] != 0
    if ( lsticc[blockIdx.x + 1] == 0 )
        return;

    lsum = 0;
    for( int j = tx; j < m; j += BLOCK_SIZE ) {
        #if (defined(PRECISION_s) || defined(PRECISION_d))
        re = dx[j];
        lsum += re*re;
        #else
        re = MAGMA_C_REAL( dx[j] );
        float im = MAGMA_C_IMAG( dx[j] );
        lsum += re*re + im*im;
        #endif
    }
    sum[tx] = lsum;
    magma_sum_reduce< BLOCK_SIZE >( tx, sum );
    
    if (tx==0)
        dxnorm[blockIdx.x] = sqrt(sum[0]);
}

extern "C" void
magmablas_scnrm2_check(
    magma_int_t m, magma_int_t n, magmaFloatComplex *da, magma_int_t ldda, 
    float *dxnorm, float *lsticc) 
{
    dim3  blocks( n );
    dim3 threads( BLOCK_SIZE );
    
    magmablas_scnrm2_check_kernel<<< blocks, threads >>>( m, da, ldda, dxnorm, lsticc );
}


//==============================================================================
__global__ void
magmablas_scnrm2_smkernel( int m, int n, magmaFloatComplex *da, int ldda,
                           float *dxnorm )
{
    const int tx = threadIdx.x;
    const int ty = threadIdx.y;
    __shared__ float sum[ BLOCK_SIZEx ][ BLOCK_SIZEy + 1];
    float re, lsum;

    for( int k = ty; k < n; k += BLOCK_SIZEy ) {
        magmaFloatComplex *dx = da + k * ldda;

        // get norm of dx
        lsum = 0;
        for( int j = tx; j < m; j += BLOCK_SIZEx ) {
            #if (defined(PRECISION_s) || defined(PRECISION_d))
            re = dx[j];
            lsum += re*re;
            #else
            re = MAGMA_C_REAL( dx[j] );
            float im = MAGMA_C_IMAG( dx[j] );
            lsum += re*re + im*im;
            #endif
        }
        sum[tx][ty] = lsum;
        magma_sum_reduce_2d< BLOCK_SIZEx, BLOCK_SIZEy+1 >( tx, ty, sum );

        if (tx == 0)
            dxnorm[k] = sqrt(sum[0][ty]);
        __syncthreads();
    }
}


//==============================================================================
/*
    Compute the scnrm2 of each column of m-by-n matrix dA.
    The resulting norms are written in the dxnorm array.
    This routine uses only one SM (block).
*/
extern "C" void
magmablas_scnrm2_sm(
    magma_int_t m, magma_int_t n, magmaFloatComplex *da, magma_int_t ldda,
    float *dxnorm)
{
    dim3  blocks( 1 );
    dim3 threads( BLOCK_SIZEx, BLOCK_SIZEy );

    magmablas_scnrm2_smkernel<<< blocks, threads, 0, magma_stream >>>( m, n, da, ldda, dxnorm );
}

//==============================================================================
extern "C"
__global__ void
magma_scnrm2_adjust_kernel(float *xnorm, magmaFloatComplex *c)
{
    const int tx = threadIdx.x;

    __shared__ float sum[ BLOCK_SIZE ];
    float temp;

    temp = MAGMA_C_ABS( c[tx] ) / xnorm[0];
    sum[tx] = -temp * temp;
    magma_sum_reduce_n( blockDim.x, tx, sum );

    __syncthreads();
    if (tx == 0)
        xnorm[0] = xnorm[0] * sqrt(1+sum[0]);
}


/*
    Adjust the norm of c to give the norm of c[k+1:], assumin that
    c was changed with orthogonal transformations.
*/
extern "C" void
magmablas_scnrm2_adjust(magma_int_t k, float *xnorm, magmaFloatComplex *c)
{
    magma_scnrm2_adjust_kernel<<< 1, k, 0, magma_stream >>> (xnorm, c);
}

//==============================================================================

#define BS 256

__global__ void
magma_scnrm2_row_check_adjust_kernel(int n, float tol, float *xnorm, float *xnorm2, 
                                     magmaFloatComplex *c, int ldc, float *lsticc)
{
    const int tx = threadIdx.x + blockIdx.x*BS;
    lsticc[tx+1] = 0;

    if (tx < n) {
        float temp = MAGMA_C_ABS( c[tx*ldc] ) / xnorm[tx];
        temp = max( 0.0, ((1.0 + temp) * (1.0 - temp)) );
        
        
        float temp2 = xnorm[tx] / xnorm2[tx];
        temp2 = temp * (temp2 * temp2);
        
        if (temp2 <= tol) {
            lsticc[tx+1] = 1;
        } else {
            xnorm[tx] *= sqrt(temp);
        }
    }
    if (tx == 0)
        lsticc[0] = 0;
    magma_sum_reduce_n( blockDim.x, tx, lsticc );
}

/*
    Adjust the norm of c[,1:k] to give the norm of c[k+1:,1:k], assuming that
    c was changed with orthogonal transformations.
    It also do checks for QP3
*/
extern "C" void
magmablas_scnrm2_row_check_adjust(
    magma_int_t k, float tol, float *xnorm, float *xnorm2, 
    magmaFloatComplex *c, magma_int_t ldc, float *lsticc)
{
    int nblocks = (k+BS-1)/BS;
    magma_scnrm2_row_check_adjust_kernel<<< nblocks, BS >>> (k, tol, xnorm, xnorm2, c, ldc, lsticc);
}

//==============================================================================

/*
    Compute the scnrm2 of each column of m-by-n matrix dA.
    The resulting norms are written in the dxnorm array. 
    The computation can be done using n blocks (default) or on one SM (commented).
*/
extern "C" void
magmablas_scnrm2_cols(
    magma_int_t m, magma_int_t n,
    magmaFloatComplex *da, magma_int_t ldda, 
    float *dxnorm) 
{
    dim3  blocks( n );
    dim3 threads( BLOCK_SIZE );
    
    magmablas_scnrm2_kernel<<< blocks, threads, 0, magma_stream >>>( m, da, ldda, dxnorm );

    // The following would do the computation on one SM
    // magmablas_scnrm2_sm(m, n, da, ldda, dxnorm);
}

//==============================================================================
