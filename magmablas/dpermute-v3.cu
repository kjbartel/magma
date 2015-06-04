/*
    -- MAGMA (version 1.2.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       May 2012

       @generated d Tue May 15 18:18:01 2012

*/
#include "common_magma.h"

#define BLOCK_SIZE 64

typedef struct {
        double *A;
        int n, lda, j0;
        short ipiv[BLOCK_SIZE];
} dlaswp_params_t;

typedef struct {
        double *A;
        int n, lda, j0, npivots;
        short ipiv[BLOCK_SIZE];
} dlaswp_params_t2;

/*********************************************************
 *
 * LAPACK Swap: permute a set of lines following ipiv
 *
 ********************************************************/
typedef struct {
    double *A;
    int n, ldx, ldy, j0, npivots;
    short ipiv[BLOCK_SIZE];
} dlaswpx_params_t;


extern "C" void dlaswp3( dlaswp_params_t2 &params );

extern "C" void 
magmablas_dpermute_long3( double *dAT, int lda, int *ipiv, int nb, int ind )
{
        int k;
        for( k = 0; k < nb-BLOCK_SIZE; k += BLOCK_SIZE )
        {
                dlaswp_params_t2 params = { dAT, lda, lda, ind + k, BLOCK_SIZE };
                for( int j = 0; j < BLOCK_SIZE; j++ )
                {
                        params.ipiv[j] = ipiv[ind + k + j] - k - 1 - ind;
                }
                    dlaswp3( params );
        }

        int num_pivots = nb - k;
        dlaswp_params_t2 params = { dAT, lda, lda, ind + k, num_pivots};
        for( int j = 0; j < num_pivots; j++ )
        {
            params.ipiv[j] = ipiv[ind + k + j] - k - 1 - ind;
        }
        dlaswp3( params );
}

#undef BLOCK_SIZE
