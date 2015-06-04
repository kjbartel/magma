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
#include "morse_starpu.h"

/*
 * Codelet CPU
 */
static void cl_zssssm_cpu_func(void *descr[], void *cl_arg)
{
    int m1;
    int n1;
    int m2;
    int n2;
    int k;
    int ib;
    PLASMA_Complex64_t *A1;
    int lda1;
    PLASMA_Complex64_t *A2;
    int lda2;
    PLASMA_Complex64_t *L1;
    int ldl1;
    PLASMA_Complex64_t *L2;
    int ldl2;
    int *IPIV;

    starpu_unpack_cl_args(cl_arg, &m1, &n1, &m2, &n2, &k, &ib, &lda1, &lda2, &ldl1, &ldl2, &IPIV);

    A1   = (PLASMA_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[0]);
    A2   = (PLASMA_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[1]);
    L1   = (PLASMA_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[2]);
    L2   = (PLASMA_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[3]);

    CORE_zssssm(m1, n1, m2, n2, k, ib, A1, lda1, A2, lda2, L1, ldl1, L2, ldl2, IPIV);
}

/*
 * Codelet Multi-cores
 */
#ifdef MORSE_USE_MULTICORE
static void cl_zssssm_mc_func(void *descr[], void *cl_arg)
{
}
#endif

/*
 * Codelet GPU
 */
#ifdef MORSE_USE_CUDA
static void cl_zssssm_cuda_func(void *descr[], void *cl_arg)
{
    int m1;
    int n1;
    int m2;
    int n2;
    int k;
    int ib;
    cuDoubleComplex *dA1;
    int lda1;
    cuDoubleComplex *dA2;
    int lda2;
    cuDoubleComplex *dL1;
    int ldl1;
    cuDoubleComplex *dL2;
    int ldl2;
    int *IPIV;
    int info;
    
    starpu_unpack_cl_args(cl_arg, &m1, &n1, &m2, &n2, &k, &ib, &lda1, &lda2, &ldl1, &ldl2, &IPIV);

    dA1  = (cuDoubleComplex *)STARPU_MATRIX_GET_PTR(descr[0]);
    dA2  = (cuDoubleComplex *)STARPU_MATRIX_GET_PTR(descr[1]);
    dL1  = (cuDoubleComplex *)STARPU_MATRIX_GET_PTR(descr[2]);
    dL2  = (cuDoubleComplex *)STARPU_MATRIX_GET_PTR(descr[3]);

#if !defined(WITHOUTTRTRI)
    dL1 += ib; /* The kernel is just using the inverted part */
#endif

    magma_zssssm_gpu(
        'C', m1, n1, m2, n2, k, ib, 
        dA1, lda1, dA2, lda2, 
        dL1, ldl1, dL2, ldl2,
        IPIV, &info);

    /* Cedric: same question than getrf ? */
    cudaThreadSynchronize();
}
#endif

/*
 * Codelet definition
 */
CODELETS(zssssm, 4, cl_zssssm_cpu_func, cl_zssssm_cuda_func, cl_zssssm_cpu_func)

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
    starpu_codelet *zssssm_codelet;
    void (*callback)(void*) = options->profiling ? cl_zssssm_callback : NULL;
    int lda1 = BLKLDD( A1, A1m );
    int lda2 = BLKLDD( A2, A2m );
    int ldl1 = BLKLDD( L1, L1m );
    int ldl2 = BLKLDD( L2, L2m );

#ifdef MORSE_USE_MULTICORE
    zssssm_codelet = options->parallel ? &cl_zssssm_mc : &cl_zssssm;
#else
    zssssm_codelet = &cl_zssssm;
#endif

    starpu_Insert_Task(zssssm_codelet,
                       VALUE,  &m1,     sizeof(int),
                       VALUE,  &n1,     sizeof(int),
                       VALUE,  &m2,     sizeof(int),
                       VALUE,  &n2,     sizeof(int),
                       VALUE,  &k,      sizeof(int),
                       VALUE,  &ib,     sizeof(int),
                       INOUT,  BLKADDR( A1, PLASMA_Complex64_t, A1m, A1n ),
                       VALUE,  &lda1,   sizeof(int),
                       INOUT,  BLKADDR( A2, PLASMA_Complex64_t, A2m, A2n ),
                       VALUE,  &lda2,   sizeof(int),
                       INPUT,  BLKADDR( L1, PLASMA_Complex64_t, L1m, L1n ),
                       VALUE,  &ldl1,   sizeof(int),
                       INPUT,  BLKADDR( L2, PLASMA_Complex64_t, L2m, L2n ),
                       VALUE,  &ldl2,   sizeof(int),
                       VALUE,  &IPIV,   sizeof(int*),
                       PRIORITY,     options->priority,
                       CALLBACK,     callback, NULL,
                       0);
}
