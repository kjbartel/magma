/**
 *
 *  @file codelet_zbrdalg.c
 *
 *  MAGMA codelets kernel
 *  MAGMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley, Univ. of Colorado Denver
 *  and INRIA Bordeaux Sud-Ouest
 *
 * @version 2.4.2
 * @author Azzam Haidar
 * @date 2011-05-15
 * @precisions normal z -> c d s
 *
 **/
#include "morse_quark.h"

/***************************************************************************//**
 *
 **/
void MORSE_zbrdalg(MorseOption_t * options,
                        int uplo,
                        int N, int NB,
                        magma_desc_t *A, int Am, int An,
                        magma_desc_t *V, int Vm, int Vn,
                        magma_desc_t *TAU, int TAUm, int TAUn,
                        int i, int j, int m, int grsiz, int BAND,
                        int *PCOL, int *ACOL, int *MCOL)
{
    QUARK_Insert_Task(options->quark, CORE_zbrdalg_quark,   options->task_flags,
        sizeof(int),               &uplo,               VALUE,
        sizeof(int),                  &N,               VALUE,
        sizeof(int),                 &NB,               VALUE,
        sizeof(PLASMA_desc),          BLKADDR(A, PLASMA_Complex64_t, Am, An),               NODEP,
        sizeof(PLASMA_Complex64_t),   BLKADDR(V, PLASMA_Complex64_t, Vm, Vn),               NODEP,
        sizeof(PLASMA_Complex64_t),   BLKADDR(TAU, PLASMA_Complex64_t, TAUm, TAUn),               NODEP,
        sizeof(int),                  &i,               VALUE,
        sizeof(int),                  &j,               VALUE,
        sizeof(int),                  &m,               VALUE,
        sizeof(int),              &grsiz,               VALUE,
        sizeof(int),                PCOL,               INPUT,
        sizeof(int),                ACOL,               INPUT,
        sizeof(int),                MCOL,              OUTPUT | LOCALITY,
        0);
}

