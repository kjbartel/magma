/**
 *
 *  @file codelet_zplgsy.c
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
#include "morse_quark.h"
/*
 * Wrapper
 */
void MORSE_zplgsy(MorseOption_t *options, 
                          PLASMA_Complex64_t bump, int m, int n, magma_desc_t *A, int Am, int An,
                          int bigM, int m0, int n0, unsigned long long int seed)
{
    int lda = BLKLDD( A, Am );
    DAG_CORE_PLGSY;
    QUARK_Insert_Task(options->quark, CORE_zplgsy_quark, options->task_flags,
        sizeof(PLASMA_Complex64_t),       &bump, VALUE,
        sizeof(int),                      &m,    VALUE,
        sizeof(int),                      &n,    VALUE,
        sizeof(PLASMA_Complex64_t)*lda*n, BLKADDR(A, PLASMA_Complex64_t, Am, An ),         OUTPUT,
        sizeof(int),                      &lda,  VALUE,
        sizeof(int),                      &bigM, VALUE,
        sizeof(int),                      &m0,   VALUE,
        sizeof(int),                      &n0,   VALUE,
        sizeof(unsigned long long int),   &seed, VALUE,
        0);
}
