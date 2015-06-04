/*
    -- MAGMA (version 1.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2011

       @generated d Sun Nov 13 20:48:37 2011

*/
#include "common_magma.h"

#define BLOCK_SIZE 64

/*********************************************************
 *
 * SWAP BLAS: permute to set of N elements
 *
 ********************************************************/
/*
 *  First version: line per line
 */
typedef struct {
    double *A1;
    double *A2;
    int n, lda1, lda2;
} magmagpu_dswap_params_t;

__global__ void magmagpu_dswap( magmagpu_dswap_params_t params )
{
    unsigned int x = threadIdx.x + __mul24(blockDim.x, blockIdx.x);
    unsigned int offset1 = __mul24( x, params.lda1);
    unsigned int offset2 = __mul24( x, params.lda2);
    if( x < params.n )
    {
        double *A1  = params.A1 + offset1;
        double *A2  = params.A2 + offset2;
        double temp = *A1;
        *A1 = *A2;
        *A2 = temp;
    }
}

extern "C" void 
magmablas_dswap( int n, double *dA1T, int lda1, 
                 double *dA2T, int lda2)
{
    int  blocksize = 64;
    dim3 blocks( (n+blocksize-1) / blocksize, 1, 1);
    magmagpu_dswap_params_t params = { dA1T, dA2T, n, lda1, lda2 };
    magmagpu_dswap<<< blocks, blocksize, 0, magma_stream >>>( params );
}

