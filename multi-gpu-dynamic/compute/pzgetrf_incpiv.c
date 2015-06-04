/**
 *
 * @file pzgetrf_incpiv.c
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

#if defined(MORSE_NON_EXPLICIT_COPY)
#define DIAG(m)   dA, m, m
#else
#define DIAG(m)   dD, m, 0
#endif
#define A(m, n)   dA, m, n
#define L(m, n)   dL, m, n
#define IPIV(m,n) &(IPIV[(int64_t)dA->mb*((int64_t)(m)+(int64_t)dA->mt*(int64_t)(n))])

/***************************************************************************//**
 *  Parallel tile LU factorization - dynamic scheduling
 **/
void magma_pzgetrf_incpiv(magma_desc_t *dA, magma_desc_t *dL, int *IPIV,
                          magma_sequence_t *sequence, magma_request_t *request)
{
    magma_context_t *magma;
    MorseOption_t options;
#if !defined(MORSE_NON_EXPLICIT_COPY)
    magma_desc_t descD;
    magma_desc_t *dD = &descD;
#endif

    int k, m, n;
    int tempkm, tempkn, tempmm, tempnn;
    int ib;
    int maxk = min(dA->mt, dA->nt);
    size_t h_work_size;
    size_t d_work_size;

    magma = magma_context_self();
    if (sequence->status != MAGMA_SUCCESS)
        return;

    morse_options_init( &options, magma, sequence, request );
    
#if !defined(MORSE_NON_EXPLICIT_COPY)
    magma_zdesc_alloc2( descD, dA->mb, dA->nb, maxk*dA->mb, dA->nb, 0, 0, maxk*dA->mb, dA->nb );
#endif

    ib = MAGMA_IB;

    h_work_size = sizeof(PLASMA_Complex64_t)*( 2*ib + 2*dL->nb )*2*dA->mb;
    d_work_size = sizeof(PLASMA_Complex64_t)*(   ib          )*2*dA->mb;

    morse_options_ws_alloc( &options, h_work_size, d_work_size );

    for (k = 0; k < maxk; k++) {
        tempkm = k == dA->mt-1 ? dA->m-k*dA->mb : dA->mb;
        tempkn = k == dA->nt-1 ? dA->n-k*dA->nb : dA->nb;
        MORSE_zgetrl(
            &options,
            tempkm, tempkn, ib,
            A(k, k), 
            L(k, k),
            IPIV(k, k),
            k == dA->mt-1, dA->nb*k);

#if !defined(MORSE_NON_EXPLICIT_COPY)
        if ( k < (dA->nt-1) ) {
            MORSE_zlacpy(
                &options,
                PlasmaUpperLower, tempkm, tempkn,
                A(k, k), 
                DIAG(k));
        }
#endif

        for (n = k+1; n < dA->nt; n++) {
            tempnn = n == dA->nt-1 ? dA->n-n*dA->nb : dA->nb;
            MORSE_zgessm(
                &options,
                tempkm, tempnn, tempkm, ib,
                IPIV(k, k),
                L(k, k),
                DIAG(k),
                A(k, n));
        }
        for (m = k+1; m < dA->mt; m++) {
            tempmm = m == dA->mt-1 ? dA->m-m*dA->mb : dA->mb;
            MORSE_ztstrf(
                &options,
                tempmm, tempkn, ib, dL->nb,
                A(k, k),
                A(m, k),
                L(m, k),
                IPIV(m, k),
                m == dA->mt-1, dA->nb*k);

            for (n = k+1; n < dA->nt; n++) {
                tempnn = n == dA->nt-1 ? dA->n-n*dA->nb : dA->nb;
                MORSE_zssssm(
                    &options,
                    dA->nb, tempnn, tempmm, tempnn, dA->nb, ib,
                    A(k, n),
                    A(m, n),
                    L(m, k),
                    A(m, k),
                    IPIV(m, k));
            }
        }
    }

    morse_options_ws_free( &options );
    morse_options_finalize( &options, magma );
    
#if !defined(MORSE_NON_EXPLICIT_COPY)
    magma_sequence_wait(magma, sequence);
    magma_desc_mat_free( dD );
#endif
}
