/**
 *
 *  @file codelet_ztrsm.c
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
void MORSE_ztrsm( MorseOption_t *options, 
                         int side, int uplo, int transA, int diag,
                         int m, int n, 
                         PLASMA_Complex64_t alpha, 
                         magma_desc_t *A, int Am, int An,
                         magma_desc_t *B, int Bm, int Bn)
{
    int lda = BLKLDD( A, Am );
    int ldb = BLKLDD( B, Bm );

    DAG_CORE_TRSM;
    QUARK_Insert_Task(options->quark, CORE_ztrsm_quark, options->task_flags,
        sizeof(PLASMA_enum),                &side,      VALUE,
        sizeof(PLASMA_enum),                &uplo,      VALUE,
        sizeof(PLASMA_enum),                &transA,    VALUE,
        sizeof(PLASMA_enum),                &diag,      VALUE,
        sizeof(int),                        &m,         VALUE,
        sizeof(int),                        &n,         VALUE,
        sizeof(PLASMA_Complex64_t),         &alpha,     VALUE,
        sizeof(PLASMA_Complex64_t)*A->desc.mb*A->desc.nb, BLKADDR(A, PLASMA_Complex64_t, Am, An ),                 INPUT,
        sizeof(int),                        &lda,       VALUE,
        sizeof(PLASMA_Complex64_t)*B->desc.mb*B->desc.nb, BLKADDR(B, PLASMA_Complex64_t, Bm, Bn ),                 INOUT | LOCALITY,
        sizeof(int),                        &ldb,       VALUE,
    0);
}
