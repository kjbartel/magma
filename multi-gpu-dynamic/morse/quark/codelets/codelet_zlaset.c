/**
 *
 *  @file codelet_zlaset.c
 *
 *  MAGMA codelets kernel
 *  MAGMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver,
 *  and INRIA Bordeaux Sud-Ouest
 *
 * @version 2.4.2
 * @author Hatem Ltaief
 * @date 2011-11-03
 * @precisions normal z -> c d s
 *
 **/
#include "morse_quark.h"
void MORSE_zlaset(MorseOption_t * options,
                  PLASMA_enum uplo, int M, int N,
                  PLASMA_Complex64_t alpha, PLASMA_Complex64_t beta, 
                  magma_desc_t *A, int Am, int An)
{
    int lda = BLKLDD(A, Am);

    DAG_CORE_LASET;
    QUARK_Insert_Task(options->quark, CORE_zlaset_quark, options->task_flags,
        sizeof(PLASMA_enum),                &uplo,  VALUE,
        sizeof(int),                        &M,     VALUE,
        sizeof(int),                        &N,     VALUE,
        sizeof(PLASMA_Complex64_t),         &alpha, VALUE,
        sizeof(PLASMA_Complex64_t),         &beta,  VALUE,
        sizeof(PLASMA_Complex64_t)*M*N,     BLKADDR(A, PLASMA_Complex64_t, Am, An),      OUTPUT,
        sizeof(int),                        &lda,   VALUE,
        0);
}
