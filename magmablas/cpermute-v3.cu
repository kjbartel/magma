/*
    -- MAGMA (version 1.2.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       June 2012

       @generated c Thu Jun 28 12:31:18 2012

*/
#include "common_magma.h"

#define BLOCK_SIZE 64

typedef struct {
        cuFloatComplex *A;
        int n, lda, j0;
        short ipiv[BLOCK_SIZE];
} claswp_params_t;

typedef struct {
        cuFloatComplex *A;
        int n, lda, j0, npivots;
        short ipiv[BLOCK_SIZE];
} claswp_params_t2;

/*********************************************************
 *
 * LAPACK Swap: permute a set of lines following ipiv
 *
 ********************************************************/
typedef struct {
    cuFloatComplex *A;
    int n, ldx, ldy, j0, npivots;
    short ipiv[BLOCK_SIZE];
} claswpx_params_t;


extern "C" void claswp3( claswp_params_t2 &params );

extern "C" void 
magmablas_cpermute_long3( cuFloatComplex *dAT, magma_int_t lda,
                          magma_int_t *ipiv, magma_int_t nb, magma_int_t ind )
{
        int k;
        for( k = 0; k < nb-BLOCK_SIZE; k += BLOCK_SIZE )
        {
                claswp_params_t2 params = { dAT, lda, lda, ind + k, BLOCK_SIZE };
                for( int j = 0; j < BLOCK_SIZE; j++ )
                {
                        params.ipiv[j] = ipiv[ind + k + j] - k - 1 - ind;
                }
                    claswp3( params );
        }

        int num_pivots = nb - k;
        claswp_params_t2 params = { dAT, lda, lda, ind + k, num_pivots};
        for( int j = 0; j < num_pivots; j++ )
        {
            params.ipiv[j] = ipiv[ind + k + j] - k - 1 - ind;
        }
        claswp3( params );
}

#undef BLOCK_SIZE
