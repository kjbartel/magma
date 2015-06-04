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
    PLASMA_desc B = dB->desc;
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
                for (k = 0; k < B.mt; k++) {
                    tempkm = k == 0 ? B.m-(B.mt-1)*B.mb : B.mb;
                    lalpha = k == 0 ? alpha : zone;
                    for (n = 0; n < B.nt; n++) {
                        tempnn = n == B.nt-1 ? B.n-n*B.nb : B.nb;
                        MORSE_ztrsm(
                            &options,
                            side, uplo, trans, diag,
                            tempkm, tempnn, 
                            lalpha, A(B.mt-1-k, B.mt-1-k),  /* lda * tempkm */
                                    B(B.mt-1-k,        n)); /* ldb * tempnn */
                    }
                    for (m = k+1; m < B.mt; m++) {
                        for (n = 0; n < B.nt; n++) {
                            tempnn = n == B.nt-1 ? B.n-n*B.nb : B.nb;
                            MORSE_zgemm(
                                &options,
                                PlasmaNoTrans, PlasmaNoTrans,
                                B.mb, tempnn, tempkm, 
                                mzone,  A(B.mt-1-m, B.mt-1-k),
                                        B(B.mt-1-k, n       ),
                                lalpha, B(B.mt-1-m, n       ));
                        }
                    }
                }
            }
            /*
             *  PlasmaLeft / PlasmaUpper / Plasma[Conj]Trans
             */
            else {
                for (k = 0; k < B.mt; k++) {
                    tempkm = k == B.mt-1 ? B.m-k*B.mb : B.mb;
                    lalpha = k == 0 ? alpha : zone;
                    for (n = 0; n < B.nt; n++) {
                        tempnn = n == B.nt-1 ? B.n-n*B.nb : B.nb;
                        MORSE_ztrsm(
                            &options,
                            side, uplo, trans, diag,
                            tempkm, tempnn, 
                            lalpha, A(k, k),
                                    B(k, n));
                    }
                    for (m = k+1; m < B.mt; m++) {
                        tempmm = m == B.mt-1 ? B.m-m*B.mb : B.mb;
                        for (n = 0; n < B.nt; n++) {
                            tempnn = n == B.nt-1 ? B.n-n*B.nb : B.nb;
                            MORSE_zgemm(
                                &options,
                                trans, PlasmaNoTrans,
                                tempmm, tempnn, B.mb, 
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
                for (k = 0; k < B.mt; k++) {
                    tempkm = k == B.mt-1 ? B.m-k*B.mb : B.mb;
                    lalpha = k == 0 ? alpha : zone;
                    for (n = 0; n < B.nt; n++) {
                        tempnn = n == B.nt-1 ? B.n-n*B.nb : B.nb;
                        MORSE_ztrsm(
                            &options,
                            side, uplo, trans, diag,
                            tempkm, tempnn, 
                            lalpha, A(k, k),
                                    B(k, n));
                    }
                    for (m = k+1; m < B.mt; m++) {
                        tempmm = m == B.mt-1 ? B.m-m*B.mb : B.mb;
                        for (n = 0; n < B.nt; n++) {
                            tempnn = n == B.nt-1 ? B.n-n*B.nb : B.nb;
                            MORSE_zgemm(
                                &options,
                                PlasmaNoTrans, PlasmaNoTrans,
                                tempmm, tempnn, B.mb, 
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
                for (k = 0; k < B.mt; k++) {
                    tempkm = k == 0 ? B.m-(B.mt-1)*B.mb : B.mb;
                    lalpha = k == 0 ? alpha : zone;
                    for (n = 0; n < B.nt; n++) {
                        tempnn = n == B.nt-1 ? B.n-n*B.nb : B.nb;
                        MORSE_ztrsm(
                            &options,
                            side, uplo, trans, diag,
                            tempkm, tempnn, 
                            lalpha, A(B.mt-1-k, B.mt-1-k),
                                    B(B.mt-1-k,        n));
                    }
                    for (m = k+1; m < B.mt; m++) {
                        tempmm = m == B.mt-1 ? B.m-m*B.mb : B.mb;
                        for (n = 0; n < B.nt; n++) {
                            tempnn = n == B.nt-1 ? B.n-n*B.nb : B.nb;
                            MORSE_zgemm(
                                &options,
                                trans, PlasmaNoTrans,
                                B.mb, tempnn, tempkm, 
                                mzone,  A(B.mt-1-k, B.mt-1-m),
                                        B(B.mt-1-k, n       ),
                                lalpha, B(B.mt-1-m, n       ));
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
                for (k = 0; k < B.nt; k++) {
                    tempkn = k == B.nt-1 ? B.n-k*B.nb : B.nb;
                    lalpha = k == 0 ? alpha : zone;
                    for (m = 0; m < B.mt; m++) {
                        tempmm = m == B.mt-1 ? B.m-m*B.mb : B.mb;
                        MORSE_ztrsm(
                            &options,
                            side, uplo, trans, diag,
                            tempmm, tempkn, 
                            lalpha, A(k, k),  /* lda * tempkn */
                                    B(m, k)); /* ldb * tempkn */
                    }
                    for (m = 0; m < B.mt; m++) {
                        tempmm = m == B.mt-1 ? B.m-m*B.mb : B.mb;
                        for (n = k+1; n < B.nt; n++) {
                            tempnn = n == B.nt-1 ? B.n-n*B.nb : B.nb;
                            MORSE_zgemm(
                                &options,
                                PlasmaNoTrans, PlasmaNoTrans,
                                tempmm, tempnn, B.mb, 
                                mzone,  B(m, k),  /* ldb * B.mb   */
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
                for (k = 0; k < B.nt; k++) {
                    tempkn = k == 0 ? B.n-(B.nt-1)*B.nb : B.nb;
                    for (m = 0; m < B.mt; m++) {
                        tempmm = m == B.mt-1 ? B.m-m*B.mb : B.mb;
                        MORSE_ztrsm(
                            &options,
                            side, uplo, trans, diag,
                            tempmm, tempkn, 
                            alpha, A(B.nt-1-k, B.nt-1-k),  /* lda * tempkn */
                                   B(       m, B.nt-1-k)); /* ldb * tempkn */

                        for (n = k+1; n < B.nt; n++) {
                            MORSE_zgemm(
                                &options,
                                PlasmaNoTrans, trans,
                                tempmm, B.nb, tempkn, 
                                minvalpha, B(m,        B.nt-1-k),  /* ldb  * tempkn */
                                           A(B.nt-1-n, B.nt-1-k), /* A.mb * tempkn (Never last row) */
                                zone,      B(m,        B.nt-1-n)); /* ldb  * B.nb   */
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
                for (k = 0; k < B.nt; k++) {
                    tempkn = k == 0 ? B.n-(B.nt-1)*B.nb : B.nb;
                    lalpha = k == 0 ? alpha : zone;
                    for (m = 0; m < B.mt; m++) {
                        tempmm = m == B.mt-1 ? B.m-m*B.mb : B.mb;
                        MORSE_ztrsm(
                            &options,
                            side, uplo, trans, diag,
                            tempmm, tempkn, 
                            lalpha, A(B.nt-1-k, B.nt-1-k),  /* lda * tempkn */
                                    B(       m, B.nt-1-k)); /* ldb * tempkn */

                        for (n = k+1; n < B.nt; n++) {
                            MORSE_zgemm(
                                &options,
                                PlasmaNoTrans, PlasmaNoTrans,
                                tempmm, B.nb, tempkn, 
                                mzone,  B(m,        B.nt-1-k),  /* ldb * tempkn */
                                        A(B.nt-1-k, B.nt-1-n),  /* lda * B.nb   */
                                lalpha, B(m,        B.nt-1-n)); /* ldb * B.nb   */
                        }
                    }
                }
            }
            /*
             *  PlasmaRight / PlasmaLower / Plasma[Conj]Trans
             */
            else {
                for (k = 0; k < B.nt; k++) {
                    tempkn = k == B.nt-1 ? B.n-k*B.nb : B.nb;
                    for (m = 0; m < B.mt; m++) {
                        tempmm = m == B.mt-1 ? B.m-m*B.mb : B.mb;
                        MORSE_ztrsm(
                            &options,
                            side, uplo, trans, diag,
                            tempmm, tempkn, 
                            alpha, A(k, k),  /* lda * tempkn */
                                   B(m, k)); /* ldb * tempkn */

                        for (n = k+1; n < B.nt; n++) {
                            tempnn = n == B.nt-1 ? B.n-n*B.nb : B.nb;
                            MORSE_zgemm(
                                &options,
                                PlasmaNoTrans, trans,
                                tempmm, tempnn, B.mb, 
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
