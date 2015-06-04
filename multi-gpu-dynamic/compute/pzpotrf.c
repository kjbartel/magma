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
    PLASMA_desc A = dA->desc;
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
        for (k = 0; k < A.mt; k++) {
            tempkm = k == A.mt-1 ? A.m-k*A.mb : A.mb;

            MORSE_zpotrf(
                &options,
                PlasmaLower, tempkm,
                A(k, k), A.nb*k);
            
            for (m = k+1; m < A.mt; m++) {
                tempmm = m == A.mt-1 ? A.m-m*A.mb : A.mb;
                MORSE_ztrsm(
                    &options,
                    PlasmaRight, PlasmaLower, PlasmaConjTrans, PlasmaNonUnit,
                    tempmm, A.mb,
                    zone, A(k, k),
                          A(m, k));
            }

            for (m = k+1; m < A.mt; m++) {
                tempmm = m == A.mt-1 ? A.m-m*A.mb : A.mb;
                MORSE_zherk(
                    &options,
                    PlasmaLower, PlasmaNoTrans,
                    tempmm, A.mb,
                    -1.0, A(m, k),
                     1.0, A(m, m));

                for (n = k+1; n < m; n++) {
                    MORSE_zgemm(
                        &options,
                        PlasmaNoTrans, PlasmaConjTrans,
                        tempmm, A.mb, A.mb,
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
        for (k = 0; k < A.nt; k++) {
            tempkm = k == A.nt-1 ? A.n-k*A.nb : A.nb;

            MORSE_zpotrf(
                &options,
                PlasmaUpper, tempkm,
                A(k, k),
                A.nb*k);

            for (m = k+1; m < A.nt; m++) {
                tempmm = m == A.nt-1 ? A.n-m*A.nb : A.nb;
                MORSE_ztrsm(
                    &options,
                    PlasmaLeft, PlasmaUpper, PlasmaConjTrans, PlasmaNonUnit,
                    A.nb, tempmm,
                    zone, A(k, k),
                          A(k, m));
            }

            for (m = k+1; m < A.nt; m++) {
                tempmm = m == A.nt-1 ? A.n-m*A.nb : A.nb;
                MORSE_zherk(
                    &options,
                    PlasmaUpper, PlasmaConjTrans,
                    tempmm, A.mb,
                    -1.0, A(k, m),
                     1.0, A(m, m));

                for (n = k+1; n < m; n++) {
                    MORSE_zgemm(
                        &options,
                        PlasmaConjTrans, PlasmaNoTrans,
                        A.mb, tempmm, A.mb,
                        mzone, A(k, n),
                               A(k, m),
                        zone,  A(n, m));
                }
            }
        }
    }

    morse_options_finalize( &options, magma );
}
