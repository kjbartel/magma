/*
    -- MAGMA (version 1.2.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       June 2012

       @generated d Thu Jun 28 12:31:19 2012

*/
#include "common_magma.h"

#define BLOCK_SIZE 64

/*********************************************************/
/*
*  Blocked version: swap several pair of line
 */
typedef struct {
    double *A1;
    double *A2;
    int n, lda1, lda2, npivots;
    short ipiv[BLOCK_SIZE];
} magmagpu_dswapblk_params_t;

__global__ void magmagpu_dswapblkrm( magmagpu_dswapblk_params_t params )
{
    unsigned int y = threadIdx.x + blockDim.x*blockIdx.x;
    if( y < params.n )
    {
        double *A1 = params.A1 + y - params.lda1;
        double *A2 = params.A2 + y;
      
        for( int i = 0; i < params.npivots; i++ )
        {
            A1 += params.lda1;
            if ( params.ipiv[i] == -1 )
                continue;
            double tmp1  = *A1;
            double *tmp2 = A2 + params.ipiv[i]*params.lda2;
            *A1   = *tmp2;
            *tmp2 = tmp1;
        }
    }
}

__global__ void magmagpu_dswapblkcm( magmagpu_dswapblk_params_t params )
{
    unsigned int y = threadIdx.x + blockDim.x*blockIdx.x;
    unsigned int offset1 = __mul24( y, params.lda1);
    unsigned int offset2 = __mul24( y, params.lda2);
    if( y < params.n )
    {
        double *A1 = params.A1 + offset1 - 1;
        double *A2 = params.A2 + offset2;
      
        for( int i = 0; i < params.npivots; i++ )
        {
            A1++;
            if ( params.ipiv[i] == -1 )
                continue;
            double tmp1  = *A1;
            double *tmp2 = A2 + params.ipiv[i];
            *A1   = *tmp2;
            *tmp2 = tmp1;
        }
    }
    __syncthreads();
}

extern "C" void 
magmablas_dswapblk( char storev, magma_int_t n, 
                    double *dA1T, magma_int_t lda1,
                    double *dA2T, magma_int_t lda2,
                    magma_int_t i1, magma_int_t i2,
                    magma_int_t *ipiv, magma_int_t inci, magma_int_t offset )
{
    int  blocksize = 64;
    dim3 blocks( (n+blocksize-1) / blocksize, 1, 1);
    int  k, im;

    if ( (storev == 'C') || (storev == 'c') ) {
        for( k=(i1-1); k<i2; k+=BLOCK_SIZE )
        {
            int sb = min(BLOCK_SIZE, i2-k);
            magmagpu_dswapblk_params_t params = { dA1T+k, dA2T, n, lda1, lda2, sb };
            for( int j = 0; j < sb; j++ )
            {
                im = ipiv[(k+j)*inci] - 1;
                if ( (k+j) == im)
                    params.ipiv[j] = -1;
                else
                    params.ipiv[j] = im - offset;
            }
            magmagpu_dswapblkcm<<< blocks, blocksize, 0, magma_stream >>>( params );
        }
    }else {
        for( k=(i1-1); k<i2; k+=BLOCK_SIZE )
        {
            int sb = min(BLOCK_SIZE, i2-k);
            magmagpu_dswapblk_params_t params = { dA1T+k*lda1, dA2T, n, lda1, lda2, sb };
            for( int j = 0; j < sb; j++ )
            {
                im = ipiv[(k+j)*inci] - 1;
                if ( (k+j) == im)
                    params.ipiv[j] = -1;
                else
                    params.ipiv[j] = im - offset;
            }
            magmagpu_dswapblkrm<<< blocks, blocksize, 0, magma_stream >>>( params );
        }
    }
}

