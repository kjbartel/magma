/**
 *
 *  @file codelet_zherk.c
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
void MORSE_zherk( MorseOption_t *options, 
                  int uplo, int trans,
                  int n, int k, 
                  double alpha, magma_desc_t *A, int Am, int An,
                  double beta,  magma_desc_t *C, int Cm, int Cn)
{
    int lda = BLKLDD( A, Am );
    int ldc = BLKLDD( C, Cm );
    int nb = options->nb;
    DAG_CORE_HERK;
    QUARK_Insert_Task(options->quark, CORE_zherk_quark, options->task_flags,
        sizeof(PLASMA_enum),                &uplo,      VALUE,
        sizeof(PLASMA_enum),                &trans,     VALUE,
        sizeof(int),                        &n,         VALUE,
        sizeof(int),                        &k,         VALUE,
        sizeof(double),                     &alpha,     VALUE,
        sizeof(PLASMA_Complex64_t)*nb*nb,    BLKADDR(A, PLASMA_Complex64_t, Am, An),                 INPUT,
        sizeof(int),                        &lda,       VALUE,
        sizeof(double),                     &beta,      VALUE,
        sizeof(PLASMA_Complex64_t)*nb*nb,    BLKADDR(C, PLASMA_Complex64_t, Cm, Cn),                 INOUT,
        sizeof(int),                        &ldc,       VALUE,
        0);
}

