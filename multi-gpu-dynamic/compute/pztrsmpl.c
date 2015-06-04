/**
 *
 * @file pztrsmpl.c
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.3.1
 * @author Jakub Kurzak
 * @author Hatem Ltaief
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @precisions normal z -> s d c
 *
 **/
#include "common.h"

#define A(m, n) dA, m, n
#define B(m, n) dB, m, n
#define L(m, n) dL, m, n

#define IPIV(m,n) &(IPIV[(int64_t)dA->nb*((int64_t)(m)+(int64_t)dA->mt*(int64_t)(n))])

/***************************************************************************//**
 *  Parallel forward substitution for tile LU - dynamic scheduling
 **/
void magma_pztrsmpl(magma_desc_t *dA, magma_desc_t *dB, magma_desc_t *dL, int *IPIV,
                    magma_sequence_t *sequence, magma_request_t *request)
{
    magma_context_t *magma;
    MorseOption_t options;

    int k, m, n;
    int tempkm, tempnn, tempkmin, tempmm, tempkn;
    int ib;
    int maxk = min(dA->mt, dA->nt);

    magma = magma_context_self();
    if (sequence->status != MAGMA_SUCCESS)
        return;

    morse_options_init( &options, magma, sequence, request );

    ib = MAGMA_IB;

    for (k = 0; k < maxk; k++) {
        tempkm   = k == dA->mt-1 ? dA->m-k*dA->mb : dA->mb;
        tempkn   = k == dA->nt-1 ? dA->n-k*dA->nb : dA->nb;
        tempkmin = k == maxk-1 ? min(dA->m, dA->n)-k*dA->mb : dA->mb;
        for (n = 0; n < dB->nt; n++) {
            tempnn = n == dB->nt-1 ? dB->n-n*dB->nb : dB->nb;
            MORSE_zgessm(
                &options,
                tempkm, tempnn, tempkmin, ib,
                IPIV(k, k),
                L(k, k),
                A(k, k),
                B(k, n));
        }
        for (m = k+1; m < dA->mt; m++) {
            tempmm = m == dA->mt-1 ? dA->m-m*dA->mb : dA->mb;
            for (n = 0; n < dB->nt; n++) {
                tempnn  = n == dB->nt-1 ? dB->n-n*dB->nb : dB->nb;
                MORSE_zssssm(
                    &options,
                    dA->nb, tempnn, tempmm, tempnn, tempkn, ib,
                    B(k, n),
                    B(m, n),
                    L(m, k),
                    A(m, k),
                    IPIV(m, k));
            }
        }
    }
    
    morse_options_finalize( &options, magma );
}
