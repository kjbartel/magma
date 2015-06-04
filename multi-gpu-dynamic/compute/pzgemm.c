/**
 *
 * @file pzgemm.c
 *
 *  MAGMA auxiliary routines
 *  MAGMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.3.1
 * @author Emmanuel Agullo
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @precisions normal z -> s d c
 *
 **/
#include "common.h"

#define A(m, n) dA, m, n
#define B(m, n) dB, m, n
#define C(m, n) dC, m, n

/***************************************************************************//**
 *  Parallel tile matrix-matrix multiplication - dynamic scheduling
 **/
void magma_pzgemm(PLASMA_enum transA, PLASMA_enum transB,
                  PLASMA_Complex64_t alpha, magma_desc_t *dA, magma_desc_t *dB,
                  PLASMA_Complex64_t beta,  magma_desc_t *dC,
                  magma_sequence_t *sequence, magma_request_t *request)
{
    magma_context_t *magma;
    MorseOption_t options;
    int k, m, n;
    int tempmm, tempnn, tempkn, tempkm;

    PLASMA_Complex64_t zone = (PLASMA_Complex64_t) 1.0;
    PLASMA_Complex64_t zbeta;

    magma = magma_context_self();
    if (sequence->status != MAGMA_SUCCESS)
        return;

    morse_options_init( &options, magma, sequence, request );

    for (m = 0; m < dC->mt; m++) {
        tempmm = m == dC->mt-1 ? dC->m-m*dC->mb : dC->mb;
        for (n = 0; n < dC->nt; n++) {
            tempnn = n == dC->nt-1 ? dC->n-n*dC->nb : dC->nb;
            /*
             *  A: PlasmaNoTrans / B: PlasmaNoTrans
             */
            if (transA == PlasmaNoTrans) {
                if (transB == PlasmaNoTrans) {
                    for (k = 0; k < dA->nt; k++) {
                        tempkn = k == dA->nt-1 ? dA->n-k*dA->nb : dA->nb;
                        zbeta = k == 0 ? beta : zone;
                        MORSE_zgemm(
                            &options,
                            transA, transB,
                            tempmm, tempnn, tempkn, 
                            alpha, A(m, k), 
                                   B(k, n),
                            zbeta, C(m, n));
                    }
                }
                /*
                 *  A: PlasmaNoTrans / B: Plasma[Conj]Trans
                 */
                else {
                    for (k = 0; k < dA->nt; k++) {
                        tempkn = k == dA->nt-1 ? dA->n-k*dA->nb : dA->nb;
                        zbeta = k == 0 ? beta : zone;
                        MORSE_zgemm(
                            &options,
                            transA, transB,
                            tempmm, tempnn, tempkn, 
                            alpha, A(m, k),  /* lda * Z */
                                   B(n, k),  /* ldb * Z */
                            zbeta, C(m, n)); /* ldc * Y */
                    }
                }
            }
            /*
             *  A: Plasma[Conj]Trans / B: PlasmaNoTrans
             */
            else {
                if (transB == PlasmaNoTrans) {
                    for (k = 0; k < dA->mt; k++) {
                        tempkm = k == dA->mt-1 ? dA->m-k*dA->mb : dA->mb;
                        zbeta = k == 0 ? beta : zone;
                        MORSE_zgemm(
                            &options,
                            transA, transB,
                            tempmm, tempnn, tempkm, 
                            alpha, A(k, m),  /* lda * X */
                                   B(k, n),  /* ldb * Y */
                            zbeta, C(m, n)); /* ldc * Y */
                    }
                }
                /*
                 *  A: Plasma[Conj]Trans / B: Plasma[Conj]Trans
                 */
                else {
                    for (k = 0; k < dA->mt; k++) {
                        tempkm = k == dA->mt-1 ? dA->m-k*dA->mb : dA->mb;
                        zbeta = k == 0 ? beta : zone;
                        MORSE_zgemm(
                            &options,
                            transA, transB,
                            tempmm, tempnn, tempkm, 
                            alpha, A(k, m),  /* lda * X */
                                   B(n, k),  /* ldb * Z */
                            zbeta, C(m, n)); /* ldc * Y */
                    }
                }
            }
        }
    }
    morse_options_finalize( &options, magma );
}
