/**
 *
 *  @file pztrtri.c
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
 *  Parallel tile triangular matrix inverse - dynamic scheduling
 **/

void magma_pztrtri(PLASMA_enum uplo, PLASMA_enum diag, 
                   magma_desc_t *dA, 
                   magma_sequence_t *sequence, 
                   magma_request_t *request)
{
    magma_context_t *magma;
    MorseOption_t options;
    int k, m, n;
    int tempkn, tempmm, tempnn;

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
        for (n = 0; n < dA->nt; n++) {
            tempnn = n == dA->nt-1 ? dA->n-n*dA->nb : dA->nb;
            for (m = n+1; m < dA->mt; m++) {
                tempmm = m == dA->mt-1 ? dA->m-m*dA->mb : dA->mb;
                MORSE_ztrsm(
                    &options,
                    PlasmaRight, uplo, PlasmaNoTrans, diag,
                    tempmm, tempnn,
                    mzone, A(n, n),
                           A(m, n));
            }
            for (m = n+1; m < dA->mt; m++) {
                tempmm = m == dA->mt-1 ? dA->m-m*dA->mb : dA->mb;
                for (k = 0; k < n; k++) {
                    tempkn = k == dA->nt-1 ? dA->n-k*dA->nb : dA->nb;
                    MORSE_zgemm(
                        &options,
                        PlasmaNoTrans, PlasmaNoTrans,
                        tempmm, tempkn, tempnn,
                        zone,  A(m, n),
                               A(n, k),
                        zone,  A(m, k));
                }
            }
            for (m = 0; m < n; m++) {
                tempmm = m == dA->mt-1 ? dA->m-m*dA->mb : dA->mb;
                MORSE_ztrsm(
                        &options,
                        PlasmaLeft, uplo, PlasmaNoTrans, diag,
                        tempnn, tempmm,
                        zone, A(n, n),
                              A(n, m));
            }
            MORSE_ztrtri(
                &options,
                uplo, diag,
                tempnn,
                A(n, n),
                dA->nb*n);
        }
    }
    /*
     *  PlasmaUpper
     */
    else{
        for (m = 0; m < dA->mt; m++) {
            tempmm = m == dA->mt-1 ? dA->m-m*dA->mb : dA->mb;
            for (n = m+1; n < dA->nt; n++) {
                tempnn = n == dA->nt-1 ? dA->n-n*dA->nb : dA->nb;
                MORSE_ztrsm(
                    &options,
                    PlasmaLeft, uplo, PlasmaNoTrans, diag,
                    tempmm, tempnn,
                    mzone, A(m, m),
                           A(m, n));
            }
            for (n = 0; n < m; n++) {
                tempnn = n == dA->nt-1 ? dA->n-n*dA->nb : dA->nb;
                for (k = m+1; k < dA->nt; k++) {
                    tempkn = k == dA->nt-1 ? dA->n-k*dA->nb : dA->nb;
                    MORSE_zgemm(
                        &options,
                        PlasmaNoTrans, PlasmaNoTrans,
                        tempnn, tempkn, tempmm,
                        zone,  A(n, m),
                               A(m, k),
                        zone,  A(n, k));
                }
                MORSE_ztrsm(
                        &options,
                        PlasmaRight, uplo, PlasmaNoTrans, diag,
                        tempnn, tempmm,
                        zone, A(m, m),
                              A(n, m));
            }
            MORSE_ztrtri(
                &options,
                uplo, diag,
                tempmm,
                A(m, m),
                dA->mb*m);
        }
    }
}
