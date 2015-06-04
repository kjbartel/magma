/**
 *
 *  @file codelet_zunmqr.c
 *
 *  MAGMA codelets kernel
 *  MAGMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver,
 *  and INRIA Bordeaux Sud-Ouest
 *
 *  @version 1.0.0
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
void MORSE_zunmqr( MorseOption_t *options, 
                   int side, int trans,
                   int m, int n, int k, int ib, 
                   magma_desc_t *A, int Am, int An,
                   magma_desc_t *T, int Tm, int Tn,
                   magma_desc_t *C, int Cm, int Cn)
{
    int nb = options->nb;
    int lda = BLKLDD( A, Am );
    int ldt = BLKLDD( T, Tm );
    int ldc = BLKLDD( C, Cm );

    DAG_CORE_UNMQR;
    QUARK_Insert_Task(options->quark, CORE_zunmqr_quark,  options->task_flags,
        sizeof(PLASMA_enum),                &side,  VALUE,
        sizeof(PLASMA_enum),                &trans, VALUE,
        sizeof(int),                        &m,     VALUE,
        sizeof(int),                        &n,     VALUE,
        sizeof(int),                        &k,     VALUE,
        sizeof(int),                        &ib,    VALUE,
        sizeof(PLASMA_Complex64_t)*nb*nb,    BLKADDR(A, PLASMA_Complex64_t, Am, An ),             INPUT | QUARK_REGION_L,
        sizeof(int),                        &lda,   VALUE,
        sizeof(PLASMA_Complex64_t)*ib*nb,    BLKADDR(T, PLASMA_Complex64_t, Tm, Tn ),             INPUT,
        sizeof(int),                        &ldt,   VALUE,
        sizeof(PLASMA_Complex64_t)*nb*nb,    BLKADDR(C, PLASMA_Complex64_t, Cm, Cn ),             INOUT,
        sizeof(int),                        &ldc,   VALUE,
        sizeof(PLASMA_Complex64_t)*ib*nb,    NULL,          SCRATCH,
        sizeof(int),                        &nb,    VALUE,
        0);
}

