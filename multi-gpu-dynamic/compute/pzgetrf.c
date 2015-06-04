/**
 *
 * @file pzgetrf.c
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
#define IPIV(m,n) &(IPIV[(int64_t)A.mb*((int64_t)(m)+(int64_t)A.mt*(int64_t)(n))])

/***************************************************************************//**
 *  Parallel tile LU factorization - dynamic scheduling
 **/
void magma_pzgetrf(magma_desc_t *dA, magma_desc_t *dL, int *IPIV,
                   magma_sequence_t *sequence, magma_request_t *request)
{
    magma_context_t *magma;
    MorseOption_t options;
    PLASMA_desc A = dA->desc;
    PLASMA_desc L = dL->desc;
#if !defined(MORSE_NON_EXPLICIT_COPY)
    magma_desc_t descD;
    magma_desc_t *dD = &descD;
#endif

    int k, m, n;
    int tempkm, tempkn, tempmm, tempnn;
    int ib;
    int maxk = min(A.mt, A.nt);
    size_t h_work_size;
    size_t d_work_size;

    magma = magma_context_self();
    if (sequence->status != MAGMA_SUCCESS)
        return;

    morse_options_init( &options, magma, sequence, request );
    
#if !defined(MORSE_NON_EXPLICIT_COPY)
    magma_zdesc_alloc2( descD, A.mb, A.nb, maxk*A.mb, A.nb, 0, 0, maxk*A.mb, A.nb );
#endif

    ib = MAGMA_IB;

    h_work_size = sizeof(PLASMA_Complex64_t)*( 2*ib + 2*L.nb )*2*A.mb;
    d_work_size = sizeof(PLASMA_Complex64_t)*(   ib          )*2*A.mb;

    morse_options_ws_alloc( &options, h_work_size, d_work_size );

    for (k = 0; k < maxk; k++) {
        tempkm = k == A.mt-1 ? A.m-k*A.mb : A.mb;
        tempkn = k == A.nt-1 ? A.n-k*A.nb : A.nb;
        MORSE_zgetrl(
            &options,
            tempkm, tempkn, ib,
            A(k, k), 
            L(k, k),
            IPIV(k, k),
            k == A.mt-1, A.nb*k);

#if !defined(MORSE_NON_EXPLICIT_COPY)
        if ( k < (A.nt-1) ) {
            MORSE_zlacpy(
                &options,
                PlasmaUpperLower, tempkm, tempkn,
                A(k, k), 
                DIAG(k));
        }
#endif

        for (n = k+1; n < A.nt; n++) {
            tempnn = n == A.nt-1 ? A.n-n*A.nb : A.nb;
            MORSE_zgessm(
                &options,
                tempkm, tempnn, tempkm, ib,
                IPIV(k, k),
                L(k, k),
                DIAG(k),
                A(k, n));
        }
        for (m = k+1; m < A.mt; m++) {
            tempmm = m == A.mt-1 ? A.m-m*A.mb : A.mb;
            MORSE_ztstrf(
                &options,
                tempmm, tempkn, ib, L.nb,
                A(k, k),
                A(m, k),
                L(m, k),
                IPIV(m, k),
                m == A.mt-1, A.nb*k);

            for (n = k+1; n < A.nt; n++) {
                tempnn = n == A.nt-1 ? A.n-n*A.nb : A.nb;
                MORSE_zssssm(
                    &options,
                    A.nb, tempnn, tempmm, tempnn, A.nb, ib,
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
