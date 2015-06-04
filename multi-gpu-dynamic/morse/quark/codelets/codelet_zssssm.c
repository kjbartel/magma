/**
 *
 *  @file codelet_zssssm.c
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
void MORSE_zssssm( MorseOption_t *options, 
                          int m1, int n1, int m2, int n2, int k, int ib,
                          magma_desc_t *A1, int A1m, int A1n,
                          magma_desc_t *A2, int A2m, int A2n,
                          magma_desc_t *L1, int L1m, int L1n,
                          magma_desc_t *L2, int L2m, int L2n,
                          int  *IPIV)
{
    int lda1 = BLKLDD( A1, A1m );
    int lda2 = BLKLDD( A2, A2m );
    int ldl1 = BLKLDD( L1, L1m );
    int ldl2 = BLKLDD( L2, L2m );

    DAG_CORE_SSSSM;
    QUARK_Insert_Task(options->quark, CORE_zssssm_quark, options->task_flags,
                      sizeof(int),                        &m1,    VALUE,
                      sizeof(int),                        &n1,    VALUE,
                      sizeof(int),                        &m2,    VALUE,
                      sizeof(int),                        &n2,    VALUE,
                      sizeof(int),                        &k,     VALUE,
                      sizeof(int),                        &ib,    VALUE,
                      sizeof(PLASMA_Complex64_t)*A1->desc.mb*A1->desc.nb, BLKADDR(A1, PLASMA_Complex64_t, A1m, A1n ),            INOUT,
                      sizeof(int),                        &lda1,  VALUE,
                      sizeof(PLASMA_Complex64_t)*A2->desc.mb*A2->desc.nb, BLKADDR(A2, PLASMA_Complex64_t, A2m, A2n ),            INOUT | LOCALITY,
                      sizeof(int),                        &lda2,  VALUE,
                      sizeof(PLASMA_Complex64_t)*L1->desc.mb*L1->desc.nb, BLKADDR(L1, PLASMA_Complex64_t, L1m, L1n ),            INPUT,
                      sizeof(int),                        &ldl1,  VALUE,
                      sizeof(PLASMA_Complex64_t)*L2->desc.mb*L2->desc.nb, BLKADDR(L2, PLASMA_Complex64_t, L2m, L2n ),            INPUT,
                      sizeof(int),                        &ldl2,  VALUE,
                      sizeof(int)*options->nb,                      IPIV,          INPUT,
                      0);
}
