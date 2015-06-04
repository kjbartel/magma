/**
 *
 *  @file codelet_ztrtri.c
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
static void cl_ztrtri_cpu_func(void *descr[], void *cl_arg)
{
    PLASMA_enum uplo;
    PLASMA_enum diag;
    int N;
    PLASMA_Complex64_t *A;
    int LDA;
    int iinfo;

    A = (PLASMA_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[0]);

    starpu_codelet_unpack_args(cl_arg, &uplo, &diag, &N, &LDA, &iinfo);
    LAPACKE_ztrtri_work(LAPACK_COL_MAJOR, lapack_const(uplo), 
                        lapack_const(diag), N, A, LDA);

    /* TODO: get the INFO and cancel if != 0 */
}

/*
 * Codelet Multi-cores
 */
#ifdef MORSE_USE_MULTICORE

static void cl_ztrtri_mc_func(void *descr[], void *cl_arg)
{
    int uplo;
    int diag;
    int N;
    PLASMA_Complex64_t *A;
    int LDA;
    int INFO = 0;
    int iinfo;

    A = (PLASMA_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[0]);

    starpu_codelet_unpack_args(cl_arg, &uplo, &diag, &N, &LDA, &iinfo);
    INFO = PLASMA_ztrtri_Lapack(uplo, diag, N, A, LDA);
}
#else
#define cl_ztrtri_mc_func cl_ztrtri_cpu_func
#endif

/*
 * Codelet GPU
 */

#ifdef MORSE_USE_CUDA
static void cl_ztrtri_cuda_func(void *descr[], void *cl_arg)
{
    int uplo;
    int diag;
    int N;
    cuDoubleComplex *A;
    int LDA;
    int INFO = 0;
    int  ret;
    int iinfo;

    A = (cuDoubleComplex *)STARPU_MATRIX_GET_PTR(descr[0]);

    starpu_codelet_unpack_args(cl_arg, &uplo, &diag, &N, &LDA, &iinfo);

    ret = magma_ztrtri_gpu(
        plasma_lapack_constants[uplo][0],
        plasma_lapack_constants[diag][0],
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

CODELETS(ztrtri, 1, cl_ztrtri_cpu_func, cl_ztrtri_cuda_func, cl_ztrtri_mc_func)

/*
 * Wrapper
 */
void MORSE_ztrtri( MorseOption_t *option,
                   PLASMA_enum uplo, PLASMA_enum diag, int n,
                   magma_desc_t *A, int Am, int An, int iinfo)
{
    struct starpu_codelet *ztrtri_codelet;
    void (*callback)(void*) = option->profiling ? cl_ztrtri_callback : NULL;
    int lda = BLKLDD( A, Am );

#ifdef MORSE_USE_MULTICORE
    ztrtri_codelet = option->parallel ? &cl_ztrtri_mc : &cl_ztrtri;
#else
    ztrtri_codelet = &cl_ztrtri;
#endif

    starpu_insert_task(ztrtri_codelet,
                       STARPU_VALUE,    &uplo,  sizeof(PLASMA_enum),
                       STARPU_VALUE,    &diag,  sizeof(PLASMA_enum),
                       STARPU_VALUE,    &n,     sizeof(int),
                       STARPU_RW,       BLKADDR( A, PLASMA_Complex64_t, Am, An ),
                       STARPU_VALUE,    &lda,   sizeof(int),
                       STARPU_VALUE,    &iinfo, sizeof(int),
                       STARPU_PRIORITY, option->priority,
                       STARPU_CALLBACK, callback, NULL,
                       0);
}
