/**
 *
 * @file pzgeqrf.c
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
#define T(m, n)   dT, m, n

/***************************************************************************//**
 *  Parallel tile QR factorization - dynamic scheduling
 **/
void magma_pzgeqrf(magma_desc_t *dA, magma_desc_t *dT,
                   magma_sequence_t *sequence, magma_request_t *request)
{
    magma_context_t *magma;
    MorseOption_t options;
    PLASMA_desc A = dA->desc;
    PLASMA_desc T = dT->desc;
#if !defined(MORSE_NON_EXPLICIT_COPY)
    magma_desc_t descD;
    magma_desc_t *dD = &descD;
#endif

    int k, m, n;
    int tempkm, tempkn, tempmm, tempnn;
    int ib;
    int minMN = min(A.mt, A.nt);
    size_t h_work_size, d_work_size;

    magma = magma_context_self();
    if (sequence->status != MAGMA_SUCCESS)
      return;

    /* Be sure the global NB is equal to the one used in those matrices */
    MAGMA_NB = T.nb;

    morse_options_init( &options, magma, sequence, request );
    
#if !defined(MORSE_NON_EXPLICIT_COPY)
    magma_zdesc_alloc2( descD, A.mb, A.nb, minMN*A.mb, A.nb, 0, 0, minMN*A.mb, A.nb );
#endif

    ib = MAGMA_IB;
    h_work_size  = A.mb;       /* size of tau                           */
    h_work_size += T.nb * ib;  /* workspace required by unmqr and tsmqr */
    h_work_size *= sizeof(PLASMA_Complex64_t);
    d_work_size  = 0;

    /* magma_scratchpad_handle h_tau; */
    /* magma_scratchpad_handle h_work; */
    /* plagma_scratchpad_handle scratch_h_a, scratch_h_T; */
    /* plagma_scratchpad_handle scratch_h_D, scratch_d_D; */

    /* size_t scratch_h_a_size  = sizeof(PLASMA_Complex64_t)*A->nb*A->mb; */
    /* size_t scratch_h_T_size  = sizeof(PLASMA_Complex64_t)*T->mb*A->mb; */
    /* size_t scratch_D_size    = sizeof(PLASMA_Complex64_t)*T->mb*A->mb; */

    /* magma_alloc_scratchpad(&h_tau,  size_h_tau,  MAGMA_CPU, MAGMA_WORKER_MEM); */
    /* magma_alloc_scratchpad(&h_work, size_h_work, MAGMA_CPU, MAGMA_WORKER_MEM); */

    /* plagma_alloc_scratchpad(&scratch_h_work, scratch_work_size, PLAGMA_CUDA,            PLAGMA_HOST_MEM  ); */
    /* plagma_alloc_scratchpad(&scratch_d_D,    scratch_D_size,    PLAGMA_CUDA,            PLAGMA_WORKER_MEM); */
    /* plagma_alloc_scratchpad(&scratch_h_D,    scratch_D_size,    PLAGMA_CUDA,            PLAGMA_HOST_MEM  ); */
    /* plagma_alloc_scratchpad(&scratch_tau,    scratch_tau_size,  PLAGMA_CPU|PLAGMA_CUDA, PLAGMA_HOST_MEM  ); */
    /* plagma_alloc_scratchpad(&scratch_h_a,    scratch_h_a_size,  PLAGMA_CPU|PLAGMA_CUDA, PLAGMA_HOST_MEM  ); */
    /* plagma_alloc_scratchpad(&scratch_h_T,    scratch_h_T_size,  PLAGMA_CPU|PLAGMA_CUDA, PLAGMA_HOST_MEM  ); */

    morse_options_ws_alloc( &options, h_work_size, d_work_size );
    
    for (k = 0; k < minMN; k++) {
        tempkm = k == A.mt-1 ? A.m-k*A.mb : A.mb;
        tempkn = k == A.nt-1 ? A.n-k*A.nb : A.nb;
        MORSE_zgeqrt(
            &options,
            tempkm, tempkn, ib,
            A(k, k),
            T(k, k));

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
            MORSE_zunmqr(
                &options,
                PlasmaLeft, PlasmaConjTrans,
                tempkm, tempnn, tempkm, ib,
                DIAG(k),
                T(k, k),
                A(k, n));
        }
        for (m = k+1; m < A.mt; m++) {
            tempmm = m == A.mt-1 ? A.m-m*A.mb : A.mb;
            MORSE_ztsqrt(
                &options,
                tempmm, tempkn, ib,
                A(k, k),
                A(m, k),
                T(m, k));

            for (n = k+1; n < A.nt; n++) {
                tempnn = n == A.nt-1 ? A.n-n*A.nb : A.nb;
                MORSE_ztsmqr(
                    &options,
                    PlasmaLeft, PlasmaConjTrans,
                    A.mb, tempnn, tempmm, tempnn, A.nb, ib,
                    A(k, n),
                    A(m, n),
                    A(m, k),
                    T(m, k));
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
