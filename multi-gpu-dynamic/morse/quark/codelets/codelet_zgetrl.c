/**
 *
 *  @file codelet_zgetrl.c
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

#if (PLASMA_VERSION_MAJOR >= 2) && (PLASMA_VERSION_MINOR >= 4)
#define CORE_zgetrf CORE_zgetrf_incpiv
#endif
static void zgetrl_cpu_func(Quark * quark)
{
    int m;
    int n;
    int ib;
    PLASMA_Complex64_t *A;
    int lda;
    PLASMA_Complex64_t *L;
    int ldl;
    int *IPIV;
    magma_sequence_t* sequence;
    magma_request_t*  request;
    PLASMA_bool check_info;
    int iinfo;

    quark_unpack_args_10(quark, m, n, ib, A, lda, IPIV, sequence, request, check_info, iinfo);
    CORE_zgetrf(m, n, ib, A, lda, IPIV, &iinfo);
    if (iinfo != PLASMA_SUCCESS && check_info)
        morse_sequence_flush(quark, sequence, request, MAGMA_ERR_SEQUENCE_FLUSHED);
}

/*
 * Wrapper
 */
void MORSE_zgetrl( MorseOption_t *options, 
                          int m, int n, int ib,
                          magma_desc_t *A, int Am, int An,
                          magma_desc_t *L, int Lm, int Ln,
                          int *IPIV,
                          PLASMA_bool check, int iinfo)
{
    int lda = BLKLDD(A, Am);
    int ldl = BLKLDD(L, Lm);
    int nb  = options->nb;

    DAG_CORE_GETRF;
    QUARK_Insert_Task(options->quark, CORE_zgetrf_incpiv_quark, options->task_flags,
        sizeof(int),                        &m,             VALUE,
        sizeof(int),                        &n,             VALUE,
        sizeof(int),                        &ib,            VALUE,
        sizeof(PLASMA_Complex64_t)*nb*nb,   BLKADDR(A, PLASMA_Complex64_t, Am, An),        INOUT | LOCALITY,
        sizeof(int),                        &lda,           VALUE,
        sizeof(int)*nb,                     IPIV,                                          OUTPUT,
        sizeof(magma_sequence_t*),          &options->sequence,      VALUE,
        sizeof(magma_request_t*),           &options->request,       VALUE,
        sizeof(PLASMA_bool),                &check,         VALUE,
        sizeof(int),                        &iinfo,         VALUE,
        0);
}
