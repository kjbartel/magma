/**
 *
 *  @file codelet_zlauum.c
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
static void cl_zlauum_cpu_func(void *descr[], void *cl_arg)
{
    PLASMA_enum uplo;
    int N;
    PLASMA_Complex64_t *A;
    int LDA;

    A = (PLASMA_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[0]);

    starpu_codelet_unpack_args(cl_arg, &uplo, &N, &LDA);
    LAPACKE_zlauum_work(LAPACK_COL_MAJOR, lapack_const(uplo), N, A, LDA);
}

/*
 * Codelet Multi-cores
 */
#ifdef MORSE_USE_MULTICORE
static void cl_zlauum_mc_func(void *descr[], void *cl_arg)
{
    int uplo;
    int N;
    PLASMA_Complex64_t *A;
    int LDA;

    A = (PLASMA_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[0]);

    starpu_codelet_unpack_args(cl_arg, &uplo, &N, &LDA);
    PLASMA_zlauum_Lapack(uplo, N, A, LDA);
}
#else
#define cl_zlauum_mc_func cl_zlauum_cpu_func
#endif

/*
 * Codelet GPU
 */
#ifdef MORSE_USE_CUDA
static void cl_zlauum_cuda_func(void *descr[], void *cl_arg)
{
    int uplo;
    int N;
    cuDoubleComplex *A;
    int LDA;
    int INFO = 0;
    int  ret;

    A = (cuDoubleComplex *)STARPU_MATRIX_GET_PTR(descr[0]);

    starpu_codelet_unpack_args(cl_arg, &uplo, &N, &LDA);

    ret = magma_zlauum_gpu(
        plasma_lapack_constants[uplo][0],
        N, A, LDA, &INFO);
    if (ret != MAGMA_SUCCESS) {
        fprintf(stderr, "Error in Magma: %d\n", ret);
        exit(-1);
    }
    cudaThreadSynchronize();
}
#endif

/*
 * Codelet definition
 */
CODELETS(zlauum, 1, cl_zlauum_cpu_func, cl_zlauum_cuda_func, cl_zlauum_mc_func)

/*
 * Wrapper
 */
void MORSE_zlauum( MorseOption_t *option,
                   PLASMA_enum uplo, int n,
                   magma_desc_t *A, int Am, int An)
{
    struct starpu_codelet *zlauum_codelet;
    void (*callback)(void*) = option->profiling ? cl_zlauum_callback : NULL;
    int lda = BLKLDD( A, Am );

#ifdef MORSE_USE_MULTICORE
    zlauum_codelet = option->parallel ? &cl_zlauum_mc : &cl_zlauum;
#else
    zlauum_codelet = &cl_zlauum;
#endif

    starpu_insert_task(zlauum_codelet,
                       STARPU_VALUE,        &uplo,  sizeof(PLASMA_enum),
                       STARPU_VALUE,        &n,     sizeof(int),
                       STARPU_RW,         BLKADDR( A, PLASMA_Complex64_t, Am, An ),
                       STARPU_VALUE,        &lda,   sizeof(int),
                       STARPU_PRIORITY,     option->priority,
                       STARPU_CALLBACK,     callback, NULL,
                       0);
}
