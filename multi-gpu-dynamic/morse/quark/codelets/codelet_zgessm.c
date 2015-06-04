/**
 *
 *  @file codelet_zgessm.c
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
 * @precisions normal z -> c d s
 *
 **/
#include "morse_quark.h"
void MORSE_zgessm(MorseOption_t * options,
                       int m, int n, int k, int ib, 
                       int *IPIV,
                       magma_desc_t *L, int Lm, int Ln, 
                       magma_desc_t *D, int Dm, int Dn, 
                       magma_desc_t *A, int Am, int An)
{
    int ldd = BLKLDD(D, Dm);
    int lda = BLKLDD(A, Am);
    int nb  = options->nb;

    DAG_CORE_GESSM;
    QUARK_Insert_Task(options->quark, CORE_zgessm_quark, options->task_flags,
        sizeof(int),                        &m,     VALUE,
        sizeof(int),                        &n,     VALUE,
        sizeof(int),                        &k,     VALUE,
        sizeof(int),                        &ib,    VALUE,
        sizeof(int)*nb,                     IPIV,          INPUT,
        sizeof(PLASMA_Complex64_t)*nb*nb,   BLKADDR(L, PLASMA_Complex64_t, Lm, Ln),             INPUT | QUARK_REGION_L,
        sizeof(int),                        &ldd,   VALUE,
        sizeof(PLASMA_Complex64_t)*nb*nb,   BLKADDR(A, PLASMA_Complex64_t, Am, An),             INOUT,
        sizeof(int),                        &lda,   VALUE,
        0);
}
