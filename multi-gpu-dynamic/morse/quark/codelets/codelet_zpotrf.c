/**
 *
 *  @file codelet_zpotrf.c
 *
 *  MAGMA codelets kernel
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
#include <lapacke.h>
#include "morse_quark.h"

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_zpotrf_quark = PCORE_zpotrf_quark
#define CORE_zpotrf_quark PCORE_zpotrf_quark
#endif
void CORE_zpotrf_quark(Quark *quark)
{
    int uplo;
    int n;
    PLASMA_Complex64_t *A;
    int lda;
    magma_sequence_t *sequence;
    magma_request_t *request;
    int iinfo;

    int info;

    quark_unpack_args_7(quark, uplo, n, A, lda, sequence, request, iinfo);
    info = LAPACKE_zpotrf_work(
        LAPACK_COL_MAJOR,
        lapack_const(uplo),
        n, A, lda);
    if (sequence->status == MAGMA_SUCCESS && info != 0)
      morse_sequence_flush(quark, sequence, request, iinfo+info);
}

/*
 * Wrapper
 */
void MORSE_zpotrf( MorseOption_t *options, 
                   PLASMA_enum uplo, int n, 
                   magma_desc_t *A, int Am, int An, int iinfo )
{
    int nb = A->desc.nb;
    int lda = BLKLDD( A, Am );

    DAG_CORE_POTRF;
    QUARK_Insert_Task(
        options->quark, CORE_zpotrf_quark, options->task_flags,
        sizeof(PLASMA_enum),                &uplo,      VALUE,
        sizeof(int),                        &n,         VALUE,
        sizeof(PLASMA_Complex64_t)*nb*nb,    BLKADDR(A, PLASMA_Complex64_t, Am, An ),  INOUT,
        sizeof(int),                        &lda,       VALUE,
        sizeof(magma_sequence_t*),          &(options->sequence),  VALUE,
        sizeof(magma_request_t*),           &(options->request),   VALUE,
        sizeof(int),                        &iinfo,     VALUE,
        0);
}
