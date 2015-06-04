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
#include "morse_starpu.h"

/*
 * Codelet CPU
 */
static void cl_zgemm_cpu_func(void *descr[], void *cl_arg)
{
    int transA;
    int transB;
    int M;
    int N;
    int K;
    PLASMA_Complex64_t alpha;
    PLASMA_Complex64_t *A;
    int LDA;
    PLASMA_Complex64_t *B;
    int LDB;
    PLASMA_Complex64_t beta;
    PLASMA_Complex64_t *C;
    int LDC;

    A = (PLASMA_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[0]);
    B = (PLASMA_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[1]);
    C = (PLASMA_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[2]);

    starpu_codelet_unpack_args(cl_arg, &transA, &transB, &M, &N, &K, &alpha, &LDA, &LDB, &beta, &LDC);
    CORE_zgemm( transA, transB,
                M, N, K,
                alpha, A, LDA,
                B, LDB,
                beta, C, LDC);
}

/*
 * Codelet Multi-cores
 */
#ifdef MORSE_USE_MULTICORE
static void cl_zgemm_mc_func(void *descr[], void *cl_arg)
{
    int transA;
    int transB;
    int M;
    int N;
    int K;
    PLASMA_Complex64_t alpha;
    PLASMA_Complex64_t *A;
    int LDA;
    PLASMA_Complex64_t *B;
    int LDB;
    PLASMA_Complex64_t beta;
    PLASMA_Complex64_t *C;
    int LDC;

    A = (PLASMA_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[0]);
    B = (PLASMA_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[1]);
    C = (PLASMA_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[2]);

    starpu_codelet_unpack_args(cl_arg, &transA, &transB, &M, &N, &K, &alpha, &LDA, &LDB, &beta, &LDC);

    PLASMA_zgemm_Lapack(transA, transB, 
                        M, N, K, 
                        alpha, A, LDA, 
                        B, LDB, beta, 
                        C, LDC);
}
#else
#define cl_zgemm_mc_func cl_zgemm_cpu_func
#endif

/*
 * Codelet GPU
 */
#ifdef MORSE_USE_CUDA
static void cl_zgemm_cuda_func(void *descr[], void *cl_arg)
{
    int transA;
    int transB;
    int M;
    int N;
    int K;
    cuDoubleComplex alpha;
    cuDoubleComplex *A;
    int LDA;
    cuDoubleComplex *B;
    int LDB;
    cuDoubleComplex beta;
    cuDoubleComplex *C;
    int LDC;

    A = (cuDoubleComplex *)STARPU_MATRIX_GET_PTR(descr[0]);
    B = (cuDoubleComplex *)STARPU_MATRIX_GET_PTR(descr[1]);
    C = (cuDoubleComplex *)STARPU_MATRIX_GET_PTR(descr[2]);

    starpu_codelet_unpack_args(cl_arg, &transA, &transB, &M, &N, &K, &alpha, &LDA, &LDB, &beta, &LDC);

    cublasZgemm (
        plasma_lapack_constants[transA][0],
        plasma_lapack_constants[transB][0],
        M, N, K,
        alpha, A, LDA,
        B, LDB, beta, C, LDC);
    cudaThreadSynchronize();
}
#endif

/*
 * Codelet definition
 */
CODELETS(zgemm, 3, cl_zgemm_cpu_func, cl_zgemm_cuda_func, cl_zgemm_mc_func)


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
    struct starpu_codelet *zgemm_codelet;
    void (*callback)(void*) = options->profiling ? cl_zgemm_callback : NULL;
    int lda = BLKLDD( A, Am );
    int ldb = BLKLDD( B, Bm );
    int ldc = BLKLDD( C, Cm );

#ifdef MORSE_USE_MULTICORE
    zgemm_codelet = options->parallel ? &cl_zgemm_mc : &cl_zgemm;
#else
    zgemm_codelet = &cl_zgemm;
#endif
    
    starpu_insert_task(
        zgemm_codelet,
        STARPU_VALUE,         &transA, sizeof(PLASMA_enum),             
        STARPU_VALUE,         &transB, sizeof(PLASMA_enum),             
        STARPU_VALUE,         &m,      sizeof(int),                     
        STARPU_VALUE,         &n,      sizeof(int),                     
        STARPU_VALUE,         &k,      sizeof(int),                     
        STARPU_VALUE,         &alpha,  sizeof(PLASMA_Complex64_t),      
        STARPU_R,  BLKADDR( A, PLASMA_Complex64_t, Am, An ),
        STARPU_VALUE,         &lda,    sizeof(int),                     
        STARPU_R,  BLKADDR( B, PLASMA_Complex64_t, Bm, Bn ),
        STARPU_VALUE,         &ldb,    sizeof(int),                     
        STARPU_VALUE,         &beta,   sizeof(PLASMA_Complex64_t),      
        STARPU_RW,  BLKADDR( C, PLASMA_Complex64_t, Cm, Cn ),
        STARPU_VALUE,         &ldc,    sizeof(int),                     
        STARPU_PRIORITY,       options->priority,
        STARPU_CALLBACK,       callback, NULL,
        0);
}
