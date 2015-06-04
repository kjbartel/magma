/**
 *
 *  @file pztrsm.c
 *
 *  MAGMA compute
 *  MAGMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver,
 *  and INRIA Bordeaux Sud-Ouest
 *
 *  @version 2.3.1
 *  @author Mathieu Faverge
 *  @author Cedric Augonnet
 *  @date 2011-06-01
 *  @precisions normal z -> c d s
 *
 **/
#include "common.h"

#define A(m, n) dA, m, n
#define B(m, n) dB, m, n

/***************************************************************************//**
 *  Parallel tile triangular solve - dynamic scheduling
 **/
void magma_pztrsm(PLASMA_enum side, PLASMA_enum uplo, PLASMA_enum trans, PLASMA_enum diag,
                  PLASMA_Complex64_t alpha, magma_desc_t *dA, magma_desc_t *dB, 
                  magma_sequence_t *sequence, magma_request_t *request)
{
    magma_context_t *magma;
    MorseOption_t options;
    int k, m, n;
    int tempkm, tempkn, tempmm, tempnn;

    PLASMA_Complex64_t zone       = (PLASMA_Complex64_t) 1.0;
    PLASMA_Complex64_t mzone      = (PLASMA_Complex64_t)-1.0;
    PLASMA_Complex64_t minvalpha  = (PLASMA_Complex64_t)-1.0 / alpha;
    PLASMA_Complex64_t lalpha;

    magma = magma_context_self();
    if (sequence->status != MAGMA_SUCCESS)
        return;

    morse_options_init( &options, magma, sequence, request );

    /*
     *  PlasmaLeft / PlasmaUpper / PlasmaNoTrans
     */
    if (side == PlasmaLeft) {
        if (uplo == PlasmaUpper) {
            if (trans == PlasmaNoTrans) {
                for (k = 0; k < dB->mt; k++) {
                    tempkm = k == 0 ? dB->m-(dB->mt-1)*dB->mb : dB->mb;
                    lalpha = k == 0 ? alpha : zone;
                    for (n = 0; n < dB->nt; n++) {
                        tempnn = n == dB->nt-1 ? dB->n-n*dB->nb : dB->nb;
                        MORSE_ztrsm(
                            &options,
                            side, uplo, trans, diag,
                            tempkm, tempnn, 
                            lalpha, A(dB->mt-1-k, dB->mt-1-k),  /* lda * tempkm */
                                    B(dB->mt-1-k,        n)); /* ldb * tempnn */
                    }
                    for (m = k+1; m < dB->mt; m++) {
                        for (n = 0; n < dB->nt; n++) {
                            tempnn = n == dB->nt-1 ? dB->n-n*dB->nb : dB->nb;
                            MORSE_zgemm(
                                &options,
                                PlasmaNoTrans, PlasmaNoTrans,
                                dB->mb, tempnn, tempkm, 
                                mzone,  A(dB->mt-1-m, dB->mt-1-k),
                                        B(dB->mt-1-k, n       ),
                                lalpha, B(dB->mt-1-m, n       ));
                        }
                    }
                }
            }
            /*
             *  PlasmaLeft / PlasmaUpper / Plasma[Conj]Trans
             */
            else {
                for (k = 0; k < dB->mt; k++) {
                    tempkm = k == dB->mt-1 ? dB->m-k*dB->mb : dB->mb;
                    lalpha = k == 0 ? alpha : zone;
                    for (n = 0; n < dB->nt; n++) {
                        tempnn = n == dB->nt-1 ? dB->n-n*dB->nb : dB->nb;
                        MORSE_ztrsm(
                            &options,
                            side, uplo, trans, diag,
                            tempkm, tempnn, 
                            lalpha, A(k, k),
                                    B(k, n));
                    }
                    for (m = k+1; m < dB->mt; m++) {
                        tempmm = m == dB->mt-1 ? dB->m-m*dB->mb : dB->mb;
                        for (n = 0; n < dB->nt; n++) {
                            tempnn = n == dB->nt-1 ? dB->n-n*dB->nb : dB->nb;
                            MORSE_zgemm(
                                &options,
                                trans, PlasmaNoTrans,
                                tempmm, tempnn, dB->mb, 
                                mzone,  A(k, m),
                                        B(k, n),
                                lalpha, B(m, n));
                        }
                    }
                }
            }
        }
        /*
         *  PlasmaLeft / PlasmaLower / PlasmaNoTrans
         */
        else {
            if (trans == PlasmaNoTrans) {
                for (k = 0; k < dB->mt; k++) {
                    tempkm = k == dB->mt-1 ? dB->m-k*dB->mb : dB->mb;
                    lalpha = k == 0 ? alpha : zone;
                    for (n = 0; n < dB->nt; n++) {
                        tempnn = n == dB->nt-1 ? dB->n-n*dB->nb : dB->nb;
                        MORSE_ztrsm(
                            &options,
                            side, uplo, trans, diag,
                            tempkm, tempnn, 
                            lalpha, A(k, k),
                                    B(k, n));
                    }
                    for (m = k+1; m < dB->mt; m++) {
                        tempmm = m == dB->mt-1 ? dB->m-m*dB->mb : dB->mb;
                        for (n = 0; n < dB->nt; n++) {
                            tempnn = n == dB->nt-1 ? dB->n-n*dB->nb : dB->nb;
                            MORSE_zgemm(
                                &options,
                                PlasmaNoTrans, PlasmaNoTrans,
                                tempmm, tempnn, dB->mb, 
                                mzone,  A(m, k),
                                        B(k, n),
                                lalpha, B(m, n));
                        }
                    }
                }
            }
            /*
             *  PlasmaLeft / PlasmaLower / Plasma[Conj]Trans
             */
            else {
                for (k = 0; k < dB->mt; k++) {
                    tempkm = k == 0 ? dB->m-(dB->mt-1)*dB->mb : dB->mb;
                    lalpha = k == 0 ? alpha : zone;
                    for (n = 0; n < dB->nt; n++) {
                        tempnn = n == dB->nt-1 ? dB->n-n*dB->nb : dB->nb;
                        MORSE_ztrsm(
                            &options,
                            side, uplo, trans, diag,
                            tempkm, tempnn, 
                            lalpha, A(dB->mt-1-k, dB->mt-1-k),
                                    B(dB->mt-1-k,        n));
                    }
                    for (m = k+1; m < dB->mt; m++) {
                        tempmm = m == dB->mt-1 ? dB->m-m*dB->mb : dB->mb;
                        for (n = 0; n < dB->nt; n++) {
                            tempnn = n == dB->nt-1 ? dB->n-n*dB->nb : dB->nb;
                            MORSE_zgemm(
                                &options,
                                trans, PlasmaNoTrans,
                                dB->mb, tempnn, tempkm, 
                                mzone,  A(dB->mt-1-k, dB->mt-1-m),
                                        B(dB->mt-1-k, n       ),
                                lalpha, B(dB->mt-1-m, n       ));
                        }
                    }
                }
            }
        }
    }
    /*
     *  PlasmaRight / PlasmaUpper / PlasmaNoTrans
     */
    else {
        if (uplo == PlasmaUpper) {
            if (trans == PlasmaNoTrans) {
                for (k = 0; k < dB->nt; k++) {
                    tempkn = k == dB->nt-1 ? dB->n-k*dB->nb : dB->nb;
                    lalpha = k == 0 ? alpha : zone;
                    for (m = 0; m < dB->mt; m++) {
                        tempmm = m == dB->mt-1 ? dB->m-m*dB->mb : dB->mb;
                        MORSE_ztrsm(
                            &options,
                            side, uplo, trans, diag,
                            tempmm, tempkn, 
                            lalpha, A(k, k),  /* lda * tempkn */
                                    B(m, k)); /* ldb * tempkn */
                    }
                    for (m = 0; m < dB->mt; m++) {
                        tempmm = m == dB->mt-1 ? dB->m-m*dB->mb : dB->mb;
                        for (n = k+1; n < dB->nt; n++) {
                            tempnn = n == dB->nt-1 ? dB->n-n*dB->nb : dB->nb;
                            MORSE_zgemm(
                                &options,
                                PlasmaNoTrans, PlasmaNoTrans,
                                tempmm, tempnn, dB->mb, 
                                mzone,  B(m, k),  /* ldb * dB->mb   */
                                        A(k, n),  /* lda * tempnn */
                                lalpha, B(m, n)); /* ldb * tempnn */
                        }
                    }
                }
            }
            /*
             *  PlasmaRight / PlasmaUpper / Plasma[Conj]Trans
             */
            else {
                for (k = 0; k < dB->nt; k++) {
                    tempkn = k == 0 ? dB->n-(dB->nt-1)*dB->nb : dB->nb;
                    for (m = 0; m < dB->mt; m++) {
                        tempmm = m == dB->mt-1 ? dB->m-m*dB->mb : dB->mb;
                        MORSE_ztrsm(
                            &options,
                            side, uplo, trans, diag,
                            tempmm, tempkn, 
                            alpha, A(dB->nt-1-k, dB->nt-1-k),  /* lda * tempkn */
                                   B(       m, dB->nt-1-k)); /* ldb * tempkn */

                        for (n = k+1; n < dB->nt; n++) {
                            MORSE_zgemm(
                                &options,
                                PlasmaNoTrans, trans,
                                tempmm, dB->nb, tempkn, 
                                minvalpha, B(m,        dB->nt-1-k),  /* ldb  * tempkn */
                                           A(dB->nt-1-n, dB->nt-1-k), /* dA->mb * tempkn (Never last row) */
                                zone,      B(m,        dB->nt-1-n)); /* ldb  * dB->nb   */
                        }
                    }
                }
            }
        }
        /*
         *  PlasmaRight / PlasmaLower / PlasmaNoTrans
         */
        else {
            if (trans == PlasmaNoTrans) {
                for (k = 0; k < dB->nt; k++) {
                    tempkn = k == 0 ? dB->n-(dB->nt-1)*dB->nb : dB->nb;
                    lalpha = k == 0 ? alpha : zone;
                    for (m = 0; m < dB->mt; m++) {
                        tempmm = m == dB->mt-1 ? dB->m-m*dB->mb : dB->mb;
                        MORSE_ztrsm(
                            &options,
                            side, uplo, trans, diag,
                            tempmm, tempkn, 
                            lalpha, A(dB->nt-1-k, dB->nt-1-k),  /* lda * tempkn */
                                    B(       m, dB->nt-1-k)); /* ldb * tempkn */

                        for (n = k+1; n < dB->nt; n++) {
                            MORSE_zgemm(
                                &options,
                                PlasmaNoTrans, PlasmaNoTrans,
                                tempmm, dB->nb, tempkn, 
                                mzone,  B(m,        dB->nt-1-k),  /* ldb * tempkn */
                                        A(dB->nt-1-k, dB->nt-1-n),  /* lda * dB->nb   */
                                lalpha, B(m,        dB->nt-1-n)); /* ldb * dB->nb   */
                        }
                    }
                }
            }
            /*
             *  PlasmaRight / PlasmaLower / Plasma[Conj]Trans
             */
            else {
                for (k = 0; k < dB->nt; k++) {
                    tempkn = k == dB->nt-1 ? dB->n-k*dB->nb : dB->nb;
                    for (m = 0; m < dB->mt; m++) {
                        tempmm = m == dB->mt-1 ? dB->m-m*dB->mb : dB->mb;
                        MORSE_ztrsm(
                            &options,
                            side, uplo, trans, diag,
                            tempmm, tempkn, 
                            alpha, A(k, k),  /* lda * tempkn */
                                   B(m, k)); /* ldb * tempkn */

                        for (n = k+1; n < dB->nt; n++) {
                            tempnn = n == dB->nt-1 ? dB->n-n*dB->nb : dB->nb;
                            MORSE_zgemm(
                                &options,
                                PlasmaNoTrans, trans,
                                tempmm, tempnn, dB->mb, 
                                minvalpha, B(m, k),  /* ldb  * tempkn */
                                           A(n, k),  /* ldan * tempkn */
                                zone,      B(m, n)); /* ldb  * tempnn */
                        }
                    }
                }
            }
        }
    }

    morse_options_finalize( &options, magma );
}
