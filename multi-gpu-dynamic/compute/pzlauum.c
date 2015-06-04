/**
 *
 *  @file pzlauum.c
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
 *  Parallel UU' or L'L operation - dynamic scheduling
 **/
void magma_pzlauum(PLASMA_enum uplo, magma_desc_t *dA,
                   magma_sequence_t *sequence, 
                   magma_request_t *request)
{
    magma_context_t *magma;
    MorseOption_t options;
    int k, m, n;
    int tempkm, tempmm, tempnn;

    PLASMA_Complex64_t zone  = (PLASMA_Complex64_t) 1.0;

    magma = magma_context_self();
    if (sequence->status != MAGMA_SUCCESS)
        return;

    morse_options_init( &options, magma, sequence, request );

    /*
     *  PlasmaLower
     */
    if (uplo == PlasmaLower) {
        for (m = 0; m < dA->mt; m++) {
            tempmm = m == dA->mt-1 ? dA->m-m*dA->mb : dA->mb;
            for (n = 0; n < m; n++) {
                tempnn = n == dA->nt-1 ? dA->n-n*dA->nb : dA->nb;
                MORSE_zherk(
                    &options,
                    PlasmaLower, PlasmaConjTrans,
                    tempnn, tempmm,
                    1.0, A(m, n),
                    1.0, A(n, n));

                for (k = n+1; k < m; k++) {
                    tempkm = k == dA->mt-1 ? dA->m-k*dA->mb : dA->mb;
                    MORSE_zgemm(
                        &options,
                        PlasmaConjTrans, PlasmaNoTrans,
                        tempkm, tempnn, tempmm,
                        zone, A(m, k),
                              A(m, n),
                        zone, A(k, n) );
                }
            }
            for (n = 0; n < m; n++) {
                tempnn = n == dA->nt-1 ? dA->n-n*dA->nb : dA->nb;
                MORSE_ztrmm(
                        &options,
                        PlasmaLeft, PlasmaLower, PlasmaConjTrans, PlasmaNonUnit,
                        tempmm, tempnn,
                        zone, A(m, m),
                              A(m, n));
            }
            MORSE_zlauum(
                        &options,
                        PlasmaLower, tempmm, A(m, m));
        }
    }
    /*
     *  PlasmaUpper
     */
    else {
         for (m = 0; m < dA->mt; m++) {
            tempmm = m == dA->mt-1 ? dA->m-m*dA->mb : dA->mb;
            for (n = 0; n < m; n++) {
                tempnn = n == dA->nt-1 ? dA->n-n*dA->nb : dA->nb;
                MORSE_zherk(
                    &options,
                    PlasmaUpper, PlasmaNoTrans,
                    tempnn, tempmm,
                    1.0, A(n, m),
                    1.0, A(n, n));

                for (k = n+1; k < m; k++) {
                    tempkm = k == dA->mt-1 ? dA->m-k*dA->mb : dA->mb;
                    MORSE_zgemm(
                        &options,
                        PlasmaNoTrans, PlasmaConjTrans,
                        tempnn, tempkm, tempmm,
                        zone, A(n, m),
                              A(k, m),
                        zone, A(n, k));
                }
            }
            for (n = 0; n < m; n++) {
                tempnn = n == dA->nt-1 ? dA->n-n*dA->nb : dA->nb;
                MORSE_ztrmm(
                        &options,
                        PlasmaRight, PlasmaUpper, PlasmaConjTrans, PlasmaNonUnit,
                        tempnn, tempmm,
                        zone, A(m, m),
                              A(n, m));
            }
            MORSE_zlauum(
                        &options,
                        PlasmaUpper, tempmm, A(m, m));
        }
    }
}
