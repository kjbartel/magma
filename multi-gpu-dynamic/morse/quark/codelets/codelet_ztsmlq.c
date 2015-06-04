/**
 *
 *  @file codelet_ztsmlq.c
 *
 *  MAGMA codelets kernel
 *  MAGMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver,
 *  and INRIA Bordeaux Sud-Ouest
 *
 * @version 2.4.2
 * @author Hatem Ltaief
 * @author Mathieu Faverge
 * @author Jakub Kurzak
 * @author Azzam Haidar
 * @author Dulceneia Becker
 * @date 2011-11-03
 * @precisions normal z -> c d s
 *
 **/
#include "morse_quark.h"
void MORSE_ztsmlq(MorseOption_t * options,
                  int side, int trans,
                  int m1, int n1, int m2, int n2, int k, int ib, 
                  magma_desc_t *A1, int A1m, int A1n, 
                  magma_desc_t *A2, int A2m, int A2n, 
                  magma_desc_t *V, int Vm, int Vn,
                  magma_desc_t *T, int Tm, int Tn)
{
    int lda1 = BLKLDD(A1, A1m);
    int lda2 = BLKLDD(A2, A2m);
    int ldv = BLKLDD(V, Vm);
    int ldt = BLKLDD(T, Tm);
    int  nb = options->nb;
    int ldwork = side == PlasmaLeft ? ib : nb;

    DAG_CORE_TSMLQ;
    QUARK_Insert_Task(options->quark, CORE_ztsmlq_quark, options->task_flags,
        sizeof(PLASMA_enum),                &side,  VALUE,
        sizeof(PLASMA_enum),                &trans, VALUE,
        sizeof(int),                        &m1,    VALUE,
        sizeof(int),                        &n1,    VALUE,
        sizeof(int),                        &m2,    VALUE,
        sizeof(int),                        &n2,    VALUE,
        sizeof(int),                        &k,     VALUE,
        sizeof(int),                        &ib,    VALUE,
        sizeof(PLASMA_Complex64_t)*nb*nb,   BLKADDR(A1, PLASMA_Complex64_t, A1m, A1n),            INOUT,
        sizeof(int),                        &lda1,  VALUE,
        sizeof(PLASMA_Complex64_t)*nb*nb,   BLKADDR(A2, PLASMA_Complex64_t, A2m, A2n),            INOUT | LOCALITY,
        sizeof(int),                        &lda2,  VALUE,
        sizeof(PLASMA_Complex64_t)*nb*nb,   BLKADDR(V, PLASMA_Complex64_t, Vm, Vn),             INPUT,
        sizeof(int),                        &ldv,   VALUE,
        sizeof(PLASMA_Complex64_t)*ib*nb,   BLKADDR(T, PLASMA_Complex64_t, Tm, Tn),             INPUT,
        sizeof(int),                        &ldt,   VALUE,
        sizeof(PLASMA_Complex64_t)*ib*nb,    NULL,          SCRATCH,
        sizeof(int),                        &ldwork, VALUE,
        0);
}
