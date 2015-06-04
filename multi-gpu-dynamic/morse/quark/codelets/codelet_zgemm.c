/**
 *
 *  @file codelet_zgemm.c
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
void MORSE_zgemm( MorseOption_t *options, 
                  int transA, int transB,
                  int m, int n, int k, 
                  PLASMA_Complex64_t alpha, magma_desc_t *A, int Am, int An,
                                            magma_desc_t *B, int Bm, int Bn,
                  PLASMA_Complex64_t beta,  magma_desc_t *C, int Cm, int Cn)
{
    int lda = BLKLDD( A, Am );
    int ldb = BLKLDD( B, Bm );
    int ldc = BLKLDD( C, Cm );
    DAG_CORE_GEMM;
    QUARK_Insert_Task(options->quark, CORE_zgemm_quark,  options->task_flags,
                      sizeof(PLASMA_enum),                &transA,    VALUE,
                      sizeof(PLASMA_enum),                &transB,    VALUE,
                      sizeof(int),                        &m,         VALUE,
                      sizeof(int),                        &n,         VALUE,
                      sizeof(int),                        &k,         VALUE,
                      sizeof(PLASMA_Complex64_t),         &alpha,     VALUE,
                      sizeof(PLASMA_Complex64_t)*A->desc.mb*A->desc.nb, BLKADDR(A, PLASMA_Complex64_t, Am, An ), INPUT,
                      sizeof(int),                        &lda,       VALUE,
                      sizeof(PLASMA_Complex64_t)*B->desc.mb*B->desc.nb, BLKADDR(B, PLASMA_Complex64_t, Bm, Bn ), INPUT,
                      sizeof(int),                        &ldb,       VALUE,
                      sizeof(PLASMA_Complex64_t),         &beta,      VALUE,
                      sizeof(PLASMA_Complex64_t)*C->desc.mb*C->desc.nb, BLKADDR(C, PLASMA_Complex64_t, Cm, Cn ), INOUT | LOCALITY,
                      sizeof(int),                        &ldc,       VALUE,
                      0);
}
