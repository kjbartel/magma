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
#include "morse_starpu.h"

/*
 * Codelet CPU
 */
static void cl_zherk_cpu_func(void *descr[], void *cl_arg)
{
    PLASMA_enum uplo;
    PLASMA_enum trans;
    int N;
    int K;
    double alpha;
    PLASMA_Complex64_t *A;
    int LDA;
    double beta;
    PLASMA_Complex64_t *C;
    int LDC;

    A = (PLASMA_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[0]);
    C = (PLASMA_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[1]);

    starpu_codelet_unpack_args(cl_arg, &uplo, &trans, &N, &K, &alpha, &LDA, &beta, &LDC);
    CORE_zherk(uplo, trans,
               N, K,
               alpha, A, LDA,
               beta, C, LDC);
}

/*
 * Codelet Multi-cores
 */
#ifdef MORSE_USE_MULTICORE
static void cl_zherk_mc_func(void *descr[], void *cl_arg)
{
    PLASMA_enum uplo;
    PLASMA_enum trans;
    int N;
    int K;
    double alpha;
    PLASMA_Complex64_t *A;
    int LDA;
    double beta;
    PLASMA_Complex64_t *C;
    int LDC;

    A = (PLASMA_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[0]);
    C = (PLASMA_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[1]);

    starpu_codelet_unpack_args(cl_arg, &uplo, &trans, &N, &K, &alpha, &LDA, &beta, &LDC);

    PLASMA_zherk_Lapack(
        uplo, trans, N, K,
        alpha, A, LDA,
        beta, C, LDC);
}
#else
#define cl_zherk_mc_func cl_zherk_cpu_func
#endif

/*
 * Codelet GPU
 */
#ifdef MORSE_USE_CUDA
static void cl_zherk_cuda_func(void *descr[], void *cl_arg)
{
    int uplo;
    int trans;
    int N;
    int K;
    double alpha;
    cuDoubleComplex *A;
    int LDA;
    double beta;
    cuDoubleComplex *C;
    int LDC;

    A = (cuDoubleComplex *)STARPU_MATRIX_GET_PTR(descr[0]);
    C = (cuDoubleComplex *)STARPU_MATRIX_GET_PTR(descr[1]);

    starpu_codelet_unpack_args(cl_arg, &uplo, &trans, &N, &K, &alpha, &LDA, &beta, &LDC);

    cublasZherk (
                 plasma_lapack_constants[uplo][0], 
                 plasma_lapack_constants[trans][0], 
                 N, K, alpha, A, LDA, beta, C, LDC);

    cudaThreadSynchronize();
}
#endif

/*
 * Codelet definition
 */
CODELETS(zherk, 2, cl_zherk_cpu_func, cl_zherk_cuda_func, cl_zherk_mc_func)

/*
 * Wrapper
 */
void MORSE_zherk( MorseOption_t *options, 
                  int uplo, int trans,
                  int n, int k, 
                  double alpha, magma_desc_t *A, int Am, int An,
                  double beta,  magma_desc_t *C, int Cm, int Cn)
{
    struct starpu_codelet *zherk_codelet;
    void (*callback)(void*) = options->profiling ? cl_zherk_callback : NULL;
    int lda = BLKLDD( A, Am );
    int ldc = BLKLDD( C, Cm );

#ifdef MORSE_USE_MULTICORE
    zherk_codelet = options->parallel ? &cl_zherk_mc : &cl_zherk;
#else
    zherk_codelet = &cl_zherk;
#endif
    
    starpu_insert_task(
        zherk_codelet,
        STARPU_VALUE,         &uplo,  sizeof(PLASMA_enum),             
        STARPU_VALUE,         &trans, sizeof(PLASMA_enum),             
        STARPU_VALUE,         &n,     sizeof(int),                     
        STARPU_VALUE,         &k,     sizeof(int),                     
        STARPU_VALUE,         &alpha, sizeof(double),                  
        STARPU_R,  BLKADDR( A, PLASMA_Complex64_t, Am, An ),
        STARPU_VALUE,         &lda,   sizeof(int),                     
        STARPU_VALUE,         &beta,  sizeof(double),                  
        STARPU_RW,  BLKADDR( C, PLASMA_Complex64_t, Cm, Cn ),
        STARPU_VALUE,         &ldc,   sizeof(int),
        STARPU_PRIORITY,       options->priority,
        STARPU_CALLBACK,       callback, NULL,
        0);
}
