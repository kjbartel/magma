/**
 *
 *  @file codelet_zaxpy.c
 *
 *  MAGMA codelets kernel
 *  MAGMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley, Univ. of Colorado Denver
 *  and INRIA Bordeaux Sud-Ouest
 *
 * @version 2.4.2
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 **/
#include "morse_quark.h"
/***************************************************************************//**
 *
 **/
void MORSE_zaxpy(MorseOption_t * options,
                 int m, int n,  PLASMA_Complex64_t alpha,
                 magma_desc_t *A, int Am, int An, 
                 magma_desc_t *B, int Bm, int Bn)
{
    int lda = BLKLDD(A, Am);
    int ldb = BLKLDD(B, Bm);
    int nb = options->nb;

    DAG_CORE_AXPY;
    QUARK_Insert_Task(options->quark, CORE_zaxpy_quark, options->task_flags,
        sizeof(int),                        &m,     VALUE,
        sizeof(int),                        &n,     VALUE,
        sizeof(PLASMA_Complex64_t),         &alpha, VALUE,
        sizeof(PLASMA_Complex64_t)*nb*nb,   BLKADDR(A, PLASMA_Complex64_t, Am, An),             INPUT,
        sizeof(int),                        &lda,   VALUE,
        sizeof(PLASMA_Complex64_t)*nb*nb,   BLKADDR(B, PLASMA_Complex64_t, Bm, Bn),             INOUT,
        sizeof(int),                        &ldb,   VALUE,
        0);
}

