/**
 *
 *  @file pzpotrf.c
 *
 *  MAGMA compute
 *  MAGMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver,
 *  and INRIA Bordeaux Sud-Ouest
 *
 *  @version 2.3.1
 *  @author Mathieu Faverge
 *  @author Cedric Augonnet
 *  @author Jakub Kurzak
 *  @author Hatem Ltaief
 *  @date 2011-06-01
 *  @precisions normal z -> c d s
 *
 **/
#include "common.h"

#define A(m, n) dA, m, n

/***************************************************************************//**
 *  Parallel tile Cholesky factorization - dynamic scheduling
 **/

void magma_pzpotrf(PLASMA_enum uplo, magma_desc_t *dA, 
                   magma_sequence_t *sequence, 
                   magma_request_t *request)
{
    magma_context_t *magma;
    MorseOption_t options;
    int k, m, n;
    int tempkm, tempmm;

    PLASMA_Complex64_t zone  = (PLASMA_Complex64_t) 1.0;
    PLASMA_Complex64_t mzone = (PLASMA_Complex64_t)-1.0;

    magma = magma_context_self();
    if (sequence->status != MAGMA_SUCCESS)
        return;

    morse_options_init( &options, magma, sequence, request );

    /*
     *  PlasmaLower
     */
    if (uplo == PlasmaLower) {
        for (k = 0; k < dA->mt; k++) {
            tempkm = k == dA->mt-1 ? dA->m-k*dA->mb : dA->mb;

            MORSE_zpotrf(
                &options,
                PlasmaLower, tempkm,
                A(k, k), dA->nb*k);
            
            for (m = k+1; m < dA->mt; m++) {
                tempmm = m == dA->mt-1 ? dA->m-m*dA->mb : dA->mb;
                MORSE_ztrsm(
                    &options,
                    PlasmaRight, PlasmaLower, PlasmaConjTrans, PlasmaNonUnit,
                    tempmm, dA->mb,
                    zone, A(k, k),
                          A(m, k));
            }

            for (m = k+1; m < dA->mt; m++) {
                tempmm = m == dA->mt-1 ? dA->m-m*dA->mb : dA->mb;
                MORSE_zherk(
                    &options,
                    PlasmaLower, PlasmaNoTrans,
                    tempmm, dA->mb,
                    -1.0, A(m, k),
                     1.0, A(m, m));

                for (n = k+1; n < m; n++) {
                    MORSE_zgemm(
                        &options,
                        PlasmaNoTrans, PlasmaConjTrans,
                        tempmm, dA->mb, dA->mb,
                        mzone, A(m, k),
                               A(n, k),
                        zone,  A(m, n));
                }
            }
        }
    }
    /*
     *  PlasmaUpper
     */
    else {
        for (k = 0; k < dA->nt; k++) {
            tempkm = k == dA->nt-1 ? dA->n-k*dA->nb : dA->nb;

            MORSE_zpotrf(
                &options,
                PlasmaUpper, tempkm,
                A(k, k),
                dA->nb*k);

            for (m = k+1; m < dA->nt; m++) {
                tempmm = m == dA->nt-1 ? dA->n-m*dA->nb : dA->nb;
                MORSE_ztrsm(
                    &options,
                    PlasmaLeft, PlasmaUpper, PlasmaConjTrans, PlasmaNonUnit,
                    dA->nb, tempmm,
                    zone, A(k, k),
                          A(k, m));
            }

            for (m = k+1; m < dA->nt; m++) {
                tempmm = m == dA->nt-1 ? dA->n-m*dA->nb : dA->nb;
                MORSE_zherk(
                    &options,
                    PlasmaUpper, PlasmaConjTrans,
                    tempmm, dA->mb,
                    -1.0, A(k, m),
                     1.0, A(m, m));

                for (n = k+1; n < m; n++) {
                    MORSE_zgemm(
                        &options,
                        PlasmaConjTrans, PlasmaNoTrans,
                        dA->mb, tempmm, dA->mb,
                        mzone, A(k, n),
                               A(k, m),
                        zone,  A(n, m));
                }
            }
        }
    }

    morse_options_finalize( &options, magma );
}
