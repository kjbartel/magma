/**
 *
 *  @file codelet_zgelqt.c
 *
 *  MAGMA codelets kernel
 *  MAGMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.4.2
 * @author Hatem Ltaief
 * @author Mathieu Faverge
 * @author Jakub Kurzak
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 **/
#include "morse_quark.h"
/***************************************************************************//**
 *
 **/
void MORSE_zgelqt(MorseOption_t * options,
                  int m, int n, int ib, 
                  magma_desc_t *A, int Am, int An,
                  magma_desc_t *T, int Tm, int Tn)
{
    int lda = BLKLDD( A, Am );
    int ldt = BLKLDD( T, Tm );
    int nb  = options->nb;

    DAG_CORE_GELQT;
    QUARK_Insert_Task(options->quark, CORE_zgelqt_quark, options->task_flags,
        sizeof(int),                        &m,     VALUE,
        sizeof(int),                        &n,     VALUE,
        sizeof(int),                        &ib,    VALUE,
        sizeof(PLASMA_Complex64_t)*nb*nb,   BLKADDR(A, PLASMA_Complex64_t, Am, An),             INOUT,
        sizeof(int),                        &lda,   VALUE,
        sizeof(PLASMA_Complex64_t)*ib*nb,   BLKADDR(T, PLASMA_Complex64_t, Tm, Tn),             OUTPUT,
        sizeof(int),                        &ldt,   VALUE,
        sizeof(PLASMA_Complex64_t)*nb,       NULL,          SCRATCH,
        sizeof(PLASMA_Complex64_t)*ib*nb,    NULL,          SCRATCH,
        0);
}
