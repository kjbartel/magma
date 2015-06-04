/*
    -- MAGMA (version 1.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2011

       @precisions normal z -> s d c

*/
#include "common_magma.h"

#define BLOCK_SIZE 64

typedef struct {
        cuDoubleComplex *A;
        int n, lda, j0;
        short ipiv[BLOCK_SIZE];
} zlaswp_params_t;

typedef struct {
        cuDoubleComplex *A;
        int n, lda, j0, npivots;
        short ipiv[BLOCK_SIZE];
} zlaswp_params_t2;

/*********************************************************
 *
 * LAPACK Swap: permute a set of lines following ipiv
 *
 ********************************************************/
typedef struct {
    cuDoubleComplex *A;
    int n, ldx, ldy, j0, npivots;
    short ipiv[BLOCK_SIZE];
} zlaswpx_params_t;


extern "C" void zlaswp3( zlaswp_params_t2 &params );

extern "C" void 
magmablas_zpermute_long3( cuDoubleComplex *dAT, int lda, int *ipiv, int nb, int ind )
{
        int k;
        for( k = 0; k < nb-BLOCK_SIZE; k += BLOCK_SIZE )
        {
                zlaswp_params_t2 params = { dAT, lda, lda, ind + k, BLOCK_SIZE };
                for( int j = 0; j < BLOCK_SIZE; j++ )
                {
                        params.ipiv[j] = ipiv[ind + k + j] - k - 1 - ind;
                }
                    zlaswp3( params );
        }

        int num_pivots = nb - k;
        zlaswp_params_t2 params = { dAT, lda, lda, ind + k, num_pivots};
        for( int j = 0; j < num_pivots; j++ )
        {
            params.ipiv[j] = ipiv[ind + k + j] - k - 1 - ind;
        }
        zlaswp3( params );
}

#undef BLOCK_SIZE
