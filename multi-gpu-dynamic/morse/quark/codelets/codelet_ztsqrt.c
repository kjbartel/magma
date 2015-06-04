/**
 *
 *  @file codelet_ztsqrt.c
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
void MORSE_ztsqrt( MorseOption_t *options, 
                   int m, int n, int ib,
                   magma_desc_t *A1, int A1m, int A1n,
                   magma_desc_t *A2, int A2m, int A2n,
                   magma_desc_t *T, int Tm, int Tn) 
                   /*void *h_tau,
                   void *h_work) , */
/* void *h_work, */
/* void *h_a, */
/* void *h_T, */
/* void *h_D, */
/* void *d_D) */
{
    int nb = options->nb;
    int lda1 = BLKLDD( A1, A1m );
    int lda2 = BLKLDD( A2, A2m );
    int ldt  = BLKLDD( T,  Tm  );

    DAG_CORE_TSQRT;
    QUARK_Insert_Task(options->quark, CORE_ztsqrt_quark, options->task_flags,
                      sizeof(int),                        &m,     VALUE,
                      sizeof(int),                        &n,     VALUE,
                      sizeof(int),                        &ib,    VALUE,
                      sizeof(PLASMA_Complex64_t)*nb*nb,    BLKADDR(A1, PLASMA_Complex64_t, A1m, A1n ),            INOUT | QUARK_REGION_D | QUARK_REGION_U,
                      sizeof(int),                        &lda1,  VALUE,
                      sizeof(PLASMA_Complex64_t)*nb*nb,    BLKADDR(A2, PLASMA_Complex64_t, A2m, A2n ),            INOUT | LOCALITY,
                      sizeof(int),                        &lda2,  VALUE,
                      sizeof(PLASMA_Complex64_t)*ib*nb,    BLKADDR(T, PLASMA_Complex64_t, Tm, Tn ),             OUTPUT,
                      sizeof(int),                        &ldt,   VALUE,
                      sizeof(PLASMA_Complex64_t)*nb,       NULL,          SCRATCH,
                      sizeof(PLASMA_Complex64_t)*ib*nb,    NULL,          SCRATCH,
                      0);
}
