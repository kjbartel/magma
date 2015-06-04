/**
 *
 *  @file pzplghe.c
 *
 *  MAGMA compute
 *  MAGMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver,
 *  and INRIA Bordeaux Sud-Ouest
 *
 *  @version 2.3.1
 *  @author Mathieu Faverge
 *  @date 2011-06-01
 *  @precisions normal z -> c
 *
 **/
#include "common.h"

#define A(m, n) dA, m, n

/***************************************************************************//**
 *  Parallel tile Cholesky factorization - dynamic scheduling
 **/
void magma_pzplghe( double bump, magma_desc_t *dA, unsigned long long int seed,
                    magma_sequence_t *sequence, magma_request_t *request )
{
    magma_context_t *magma;
    MorseOption_t options;

    int m, n;
    int tempmm, tempnn;

    magma = magma_context_self();
    if (sequence->status != MAGMA_SUCCESS)
        return;

    morse_options_init( &options, magma, sequence, request );

    for (m = 0; m < dA->mt; m++) {
        tempmm = m == dA->mt-1 ? dA->m-m*dA->mb : dA->mb;

        for (n = 0; n < dA->nt; n++) {
            tempnn = n == dA->nt-1 ? dA->n-n*dA->nb : dA->nb;

            MORSE_zplghe( 
                &options,
                bump, tempmm, tempnn, A(m, n),
                dA->m, m*dA->mb, n*dA->nb, seed );
        
        }
    }
    morse_options_finalize( &options, magma );
}
