/**
 *
 *  @file pzplgsy.c
 *
 *  MAGMA compute
 *  MAGMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver,
 *  and INRIA Bordeaux Sud-Ouest
 *
 *  @version 2.3.1
 *  @author Mathieu Faverge
 *  @date 2011-06-01
 *  @precisions normal z -> c d s
 *
 **/
#include "common.h"

#define A(m, n) dA, m, n

/***************************************************************************//**
 *  Parallel tile Cholesky factorization - dynamic scheduling
 **/
void magma_pzplgsy( PLASMA_Complex64_t bump, magma_desc_t *dA, unsigned long long int seed,
                    magma_sequence_t *sequence, magma_request_t *request )
{
    magma_context_t *magma;
    MorseOption_t options;
    PLASMA_desc A = dA->desc;

    int m, n;
    int tempmm, tempnn;

    magma = magma_context_self();
    if (sequence->status != MAGMA_SUCCESS)
        return;

    morse_options_init( &options, magma, sequence, request );

    for (m = 0; m < A.mt; m++) {
        tempmm = m == A.mt-1 ? A.m-m*A.mb : A.mb;
        
        for (n = 0; n < A.nt; n++) {
            tempnn = n == A.nt-1 ? A.n-n*A.nb : A.nb;
            
            MORSE_zplgsy( 
                &options,
                bump, tempmm, tempnn, A(m, n),
                A.m, m*A.mb, n*A.nb, seed );
        }
    }
    
    morse_options_finalize( &options, magma );
}
