/**
 *
 *  @file codelet_zher2k.c
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
 * @date 2011-11-03
 * @precisions normal z -> c
 *
 **/
#include "morse_quark.h"

void MORSE_zher2k(MorseOption_t * options,
                  int uplo, int trans,
                  int n, int k,
                  PLASMA_Complex64_t alpha, magma_desc_t *A, int Am, int An,
                                            magma_desc_t *B, int Bm, int Bn,
                  double beta, magma_desc_t *C, int Cm, int Cn)
{
    int lda = BLKLDD(A, Am);
    int ldb = BLKLDD(B, Bm);
    int ldc = BLKLDD(C, Cm);
    int nb = options->nb;

    DAG_CORE_HER2K;
    QUARK_Insert_Task(options->quark, CORE_zher2k_quark, options->task_flags,
        sizeof(PLASMA_enum),                &uplo,      VALUE,
        sizeof(PLASMA_enum),                &trans,     VALUE,
        sizeof(int),                        &n,         VALUE,
        sizeof(int),                        &k,         VALUE,
        sizeof(PLASMA_Complex64_t),         &alpha,     VALUE,
        sizeof(PLASMA_Complex64_t)*nb*nb,   BLKADDR(A, PLASMA_Complex64_t, Am, An),                 INPUT,
        sizeof(int),                        &lda,       VALUE,
        sizeof(PLASMA_Complex64_t)*nb*nb,   BLKADDR(B, PLASMA_Complex64_t, Bm, Bn),                 INPUT,
        sizeof(int),                        &ldb,       VALUE,
        sizeof(double),                     &beta,      VALUE,
        sizeof(PLASMA_Complex64_t)*nb*nb,   BLKADDR(C, PLASMA_Complex64_t, Cm, Cn),                 INOUT,
        sizeof(int),                        &ldc,       VALUE,
        0);
}
