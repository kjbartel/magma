/*
    -- MAGMA (version 1.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2011

       @generated s Sun Nov 13 20:48:37 2011

*/
#include "common_magma.h"

#define BLOCK_SIZE 64

typedef struct {
        float *A;
        int n, lda, j0;
        short ipiv[BLOCK_SIZE];
} slaswp_params_t;

typedef struct {
        float *A;
        int n, lda, j0, npivots;
        short ipiv[BLOCK_SIZE];
} slaswp_params_t2;

/*********************************************************
 *
 * LAPACK Swap: permute a set of lines following ipiv
 *
 ********************************************************/
typedef struct {
    float *A;
    int n, ldx, ldy, j0, npivots;
    short ipiv[BLOCK_SIZE];
} slaswpx_params_t;


extern "C" void slaswp3( slaswp_params_t2 &params );

extern "C" void 
magmablas_spermute_long3( float *dAT, int lda, int *ipiv, int nb, int ind )
{
        int k;
        for( k = 0; k < nb-BLOCK_SIZE; k += BLOCK_SIZE )
        {
                slaswp_params_t2 params = { dAT, lda, lda, ind + k, BLOCK_SIZE };
                for( int j = 0; j < BLOCK_SIZE; j++ )
                {
                        params.ipiv[j] = ipiv[ind + k + j] - k - 1 - ind;
                }
                    slaswp3( params );
        }

        int num_pivots = nb - k;
        slaswp_params_t2 params = { dAT, lda, lda, ind + k, num_pivots};
        for( int j = 0; j < num_pivots; j++ )
        {
            params.ipiv[j] = ipiv[ind + k + j] - k - 1 - ind;
        }
        slaswp3( params );
}

#undef BLOCK_SIZE
