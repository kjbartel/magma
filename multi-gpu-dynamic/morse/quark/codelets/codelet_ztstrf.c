/**
 *
 *  @file codelet_ztstrf.c
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
static void ztstrf_cpu_func(Quark * quark)
{
    int m;
    int n;
    int ib;
    int nb;
    PLASMA_Complex64_t *U;
    int ldu;
    PLASMA_Complex64_t *A;
    int lda;
    PLASMA_Complex64_t *L;
    int ldl;
    int *IPIV;
    PLASMA_Complex64_t *WORK;
    int ldwork;
    PLASMA_bool check_info;
    magma_sequence_t *sequence;
    magma_request_t *request;
    int iinfo;
    int info;

    quark_unpack_args_17(quark, m, n, ib, nb, U, ldu, A, lda, L, ldl, IPIV, WORK,
                         ldwork, sequence, request, check_info, iinfo);

    CORE_ztstrf(m, n, ib, nb, U, ldu, A, lda, L, ldl, IPIV, WORK, ldwork, &info);
    if (info != PLASMA_SUCCESS && check_info)
       magma_sequence_flush(quark, sequence, request, iinfo + info);
}

/*
 * Wrapper
 */
void MORSE_ztstrf( MorseOption_t *options, 
                   int m, int n, int ib, int nb,
                   magma_desc_t *U, int Um, int Un,
                   magma_desc_t *A, int Am, int An,
                   magma_desc_t *L, int Lm, int Ln,
                   int  *IPIV,
                   PLASMA_bool check, int iinfo)
{
    int ldu = BLKLDD( U, Um );
    int lda = BLKLDD( A, Am );
    int ldl = BLKLDD( L, Lm );

    DAG_CORE_TSTRF;
    QUARK_Insert_Task(options->quark, CORE_ztstrf_quark, options->task_flags,
        sizeof(int),                        &m,             VALUE,
        sizeof(int),                        &n,             VALUE,
        sizeof(int),                        &ib,            VALUE,
        sizeof(int),                        &nb,            VALUE,
        sizeof(PLASMA_Complex64_t)*nb*nb,    BLKADDR(U, PLASMA_Complex64_t, Um, Un ),                     INOUT | QUARK_REGION_D | QUARK_REGION_U,
        sizeof(int),                        &ldu,           VALUE,
        sizeof(PLASMA_Complex64_t)*nb*nb,    BLKADDR(A, PLASMA_Complex64_t, Am, An ),                     INOUT | LOCALITY,
        sizeof(int),                        &lda,           VALUE,
        sizeof(PLASMA_Complex64_t)*ib*nb,    BLKADDR(L, PLASMA_Complex64_t, Lm, Ln ),                     OUTPUT,
        sizeof(int),                        &ldl,           VALUE,
        sizeof(int)*nb,                      IPIV,                  OUTPUT,
        sizeof(PLASMA_Complex64_t)*ib*nb,    NULL,                  SCRATCH,
        sizeof(int),                        &nb,            VALUE,
        sizeof(magma_sequence_t*),           &options->sequence,      VALUE,
        sizeof(magma_request_t*),            &options->request,       VALUE,
        sizeof(PLASMA_bool),                &check,         VALUE,
        sizeof(int),                        &iinfo,         VALUE,
        0);
}
