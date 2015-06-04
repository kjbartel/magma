/*
    -- MAGMA (version 1.2.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       June 2012

       @generated d Thu Jun 28 12:31:18 2012

*/
#include "common_magma.h"

#define BLOCK_SIZE 64

typedef struct {
        double *A;
        int n, lda, j0;
        short ipiv[BLOCK_SIZE];
} dlaswp_params_t;

__global__ void mydlaswp_( dlaswp_params_t params )
{
        unsigned int tid = threadIdx.x + __mul24(blockDim.x, blockIdx.x);
        if( tid < params.n )
        {
                int lda = params.lda;
                double *A = params.A + tid + lda * params.j0;

                for( int i = 0; i < BLOCK_SIZE; i++ )
                {
                         int j = params.ipiv[i];
                        double *p1 = A + i*lda;
                        double *p2 = A + j*lda;
                        double temp = *p1;
                        *p1 = *p2;
                        *p2 = temp;
                }
        }
}

extern "C" void dlaswp2( dlaswp_params_t &params )
{
         int blocksize = 64;
        dim3 blocks = (params.n+blocksize-1) / blocksize;
        mydlaswp_<<< blocks, blocksize, 0, magma_stream >>>( params );
}


extern "C" void 
magmablas_dpermute_long( double *dAT, magma_int_t lda,
                         magma_int_t *ipiv, magma_int_t nb, magma_int_t ind )
{
        // assert( (nb % BLOCK_SIZE) == 0 );
        for( int k = 0; k < nb; k += BLOCK_SIZE )
        {
                dlaswp_params_t params = { dAT, lda, lda, ind + k };
                for( int j = 0; j < BLOCK_SIZE; j++ )
                {
                        params.ipiv[j] = ipiv[ind + k + j] - k - 1;
                        ipiv[ind + k + j] += ind;
                }
                dlaswp2( params );
        }
}

#undef BLOCK_SIZE
