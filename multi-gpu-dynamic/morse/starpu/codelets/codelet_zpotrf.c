/**
 *
 *  @file codelet_zpotrf.c
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
static void cl_zpotrf_cpu_func(void *descr[], void *cl_arg)
{
    PLASMA_enum uplo;
    int N;
    PLASMA_Complex64_t *A;
    int LDA;
    int info = 0;
    int iinfo;

    A = (PLASMA_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[0]);

    starpu_codelet_unpack_args(cl_arg, &uplo, &N, &LDA, &iinfo);
    CORE_zpotrf( uplo, N, A, LDA, &info );

    /* TODO: get the INFO and cancel if != 0 */
}

/*
 * Codelet Multi-cores
 */
#ifdef MORSE_USE_MULTICORE
static void cl_zpotrf_mc_func(void *descr[], void *cl_arg)
{
    int uplo;
    int N;
    PLASMA_Complex64_t *A;
    int LDA;
    int iinfo;

    A = (PLASMA_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[0]);

    starpu_codelet_unpack_args(cl_arg, &uplo, &N, &LDA, &iinfo);
    PLASMA_zpotrf_Lapack(uplo, N, A, LDA);
}
#else
#define cl_zpotrf_mc_func cl_zpotrf_cpu_func
#endif

/*
 * Codelet GPU
 */
#ifdef MORSE_USE_CUDA
static void cl_zpotrf_cuda_func(void *descr[], void *cl_arg)
{
    int uplo;
    int N;
    cuDoubleComplex *A;
    int LDA;
    int INFO = 0;
    int  ret;
    int iinfo;

    A = (cuDoubleComplex *)STARPU_MATRIX_GET_PTR(descr[0]);

    starpu_codelet_unpack_args(cl_arg, &uplo, &N, &LDA, &iinfo);

    ret = magma_zpotrf_gpu(
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
CODELETS(zpotrf, 1, cl_zpotrf_cpu_func, cl_zpotrf_cuda_func, cl_zpotrf_mc_func)

/*
 * Wrapper
 */
void MORSE_zpotrf( MorseOption_t *option, 
                   PLASMA_enum uplo, int n, 
                   magma_desc_t *A, int Am, int An, int iinfo)
{
    struct starpu_codelet *zpotrf_codelet;
    void (*callback)(void*) = option->profiling ? cl_zpotrf_callback : NULL;
    int lda = BLKLDD( A, Am );

#ifdef MORSE_USE_MULTICORE
    zpotrf_codelet = option->parallel ? &cl_zpotrf_mc : &cl_zpotrf;
#else
    zpotrf_codelet = &cl_zpotrf;
#endif

    starpu_insert_task(zpotrf_codelet,
                       STARPU_VALUE,    &uplo,  sizeof(PLASMA_enum),
                       STARPU_VALUE,    &n,     sizeof(int),
                       STARPU_RW,       BLKADDR( A, PLASMA_Complex64_t, Am, An ),
                       STARPU_VALUE,    &lda,   sizeof(int),
                       STARPU_VALUE,    &iinfo, sizeof(int),
                       STARPU_PRIORITY, option->priority,
                       STARPU_CALLBACK, callback, NULL,
                       0);

    /* TODO: take cancellation into account */
    /* if ( *info != 0 ) */
    /*     return (*info + iinfo); */
}
