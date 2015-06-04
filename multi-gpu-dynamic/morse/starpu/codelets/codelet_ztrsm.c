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
#include "morse_starpu.h"

/*
 * Codelet CPU
 */
static void cl_ztrsm_cpu_func(void *descr[], void *cl_arg)
{
    int side;
    int uplo;
    int transA;
    int diag;
    int M;
    int N;
    PLASMA_Complex64_t alpha;
    PLASMA_Complex64_t *A;
    int LDA;
    PLASMA_Complex64_t *B;
    int LDB;

    A = (PLASMA_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[0]);
    B = (PLASMA_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[1]);

    starpu_codelet_unpack_args(cl_arg, &side, &uplo, &transA, &diag, &M, &N, &alpha, &LDA, &LDB); 
    CORE_ztrsm( side, uplo, transA, diag,
                M, N,
                alpha, A, LDA,
                B, LDB);
}

/*
 * Codelet Multi-cores
 */
#ifdef MORSE_USE_MULTICORE
static void cl_ztrsm_mc_func(void *descr[], void *cl_arg)
{
    int side;
    int uplo;
    int transA;
    int diag;
    int M;
    int N;
    PLASMA_Complex64_t alpha;
    PLASMA_Complex64_t *A;
    int LDA;
    PLASMA_Complex64_t *B;
    int LDB;

    A = (PLASMA_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[0]);
    B = (PLASMA_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[1]);

    starpu_codelet_unpack_args(cl_arg, &side, &uplo, &transA, &diag, &M, &N, &alpha, &LDA, &LDB); 

    PLASMA_ztrsm_Lapack(
        side, uplo, transA, diag,
        M, N, alpha, A, LDA, B, LDB);
}
#else
#define cl_ztrsm_mc_func cl_ztrsm_cpu_func
#endif

/*
 * Codelet GPU
 */
#ifdef MORSE_USE_CUDA
static void cl_ztrsm_cuda_func(void *descr[], void *cl_arg)
{
    int side;
    int uplo;
    int transA;
    int diag;
    int M;
    int N;
    cuDoubleComplex alpha;
    cuDoubleComplex *A;
    int LDA;
    cuDoubleComplex *B;
    int LDB;

    A = (cuDoubleComplex *)STARPU_MATRIX_GET_PTR(descr[0]);
    B = (cuDoubleComplex *)STARPU_MATRIX_GET_PTR(descr[1]);

    starpu_codelet_unpack_args(cl_arg, &side, &uplo, &transA, &diag, &M, &N, &alpha, &LDA, &LDB); 

    cublasZtrsm ( 
        plasma_lapack_constants[side][0],
        plasma_lapack_constants[uplo][0],
        plasma_lapack_constants[transA][0],
        plasma_lapack_constants[diag][0],
        M, N, alpha, A, LDA, B, LDB);
    cudaThreadSynchronize();
}
#endif

/*
 * Codelet definition
 */
CODELETS(ztrsm, 2, cl_ztrsm_cpu_func, cl_ztrsm_cuda_func, cl_ztrsm_mc_func)


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
    struct starpu_codelet *ztrsm_codelet;
    void (*callback)(void*) = options->profiling ? cl_ztrsm_callback : NULL;
    int lda = BLKLDD( A, Am );
    int ldb = BLKLDD( B, Bm );

#ifdef MORSE_USE_MULTICORE
    ztrsm_codelet = options->parallel ? &cl_ztrsm_mc : &cl_ztrsm;
#else
    ztrsm_codelet = &cl_ztrsm;
#endif

    starpu_insert_task(
        ztrsm_codelet,
        STARPU_VALUE, &side,   sizeof(PLASMA_enum),
        STARPU_VALUE, &uplo,   sizeof(PLASMA_enum),
        STARPU_VALUE, &transA, sizeof(PLASMA_enum),
        STARPU_VALUE, &diag,   sizeof(PLASMA_enum),
        STARPU_VALUE, &m,      sizeof(int),
        STARPU_VALUE, &n,      sizeof(int),
        STARPU_VALUE, &alpha,  sizeof(PLASMA_Complex64_t),
        STARPU_R, BLKADDR( A, PLASMA_Complex64_t, Am, An ),
        STARPU_VALUE, &lda,    sizeof(int),
        STARPU_RW, BLKADDR( B, PLASMA_Complex64_t, Bm, Bn ),
        STARPU_VALUE, &ldb,    sizeof(int),
        STARPU_PRIORITY,       options->priority,
        STARPU_CALLBACK,       callback, NULL,
        0);
}
