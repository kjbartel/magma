/**
 *
 *  @file codelet_zlacpy.c
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
static void cl_zlacpy_cpu_func(void *descr[], void *cl_arg)
{
    PLASMA_enum uplo;
    int M;
    int N;
    PLASMA_Complex64_t *A;
    int LDA;
    PLASMA_Complex64_t *B;
    int LDB;

    starpu_codelet_unpack_args(cl_arg, &uplo, &M, &N, &LDA, &LDB);

    A = (PLASMA_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[0]);
    B = (PLASMA_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[1]);

    LAPACKE_zlacpy_work(
        LAPACK_COL_MAJOR,
        lapack_const(uplo),
        M, N, A, LDA, B, LDB);
}

/*
 * Codelet Multi-cores
 */
#if defined(MORSE_USE_MULTICORE) && 0
static void cl_zlacpy_mc_func(void *descr[], void *cl_arg)
{
}
#endif

/*
 * Codelet GPU
 */
#ifdef MORSE_USE_CUDA
static void cl_zlacpy_cuda_func(void *descr[], void *cl_arg)
{
    PLASMA_enum uplo;
    int M;
    int N;
    cuDoubleComplex *A;
    int LDA;
    cuDoubleComplex *B;
    int LDB;

    starpu_codelet_unpack_args(cl_arg, &uplo, &M, &N, &LDA, &LDB);

    A = (cuDoubleComplex *)STARPU_MATRIX_GET_PTR(descr[0]);
    B = (cuDoubleComplex *)STARPU_MATRIX_GET_PTR(descr[1]);

    if ( M != LDA || M != LDB )
      fprintf( stderr, "WARNING: lacpy on GPU doesn;t support leading dimension different form M\n");
    cudaMemcpy( B, A, M*N*sizeof(cuDoubleComplex), cudaMemcpyDeviceToDevice );

    /* Not implemented for now, and for LU, we can cpy all the matrix */
    /* magma_zlacpy( lapack_const(uplo)[0], m, n, */
    /*               A, lda,  */
    /*               B, ldb ); */

    /* Cedric: same question than getrf ? */
    cudaThreadSynchronize();
}
#endif

/*
 * Codelet definition
 */
CODELETS(zlacpy, 2, cl_zlacpy_cpu_func, cl_zlacpy_cuda_func, cl_zlacpy_cpu_func)

/*
 * Wrapper
 */
void MORSE_zlacpy( MorseOption_t *options, 
                   PLASMA_enum uplo, int m, int n,
                   magma_desc_t *A, int Am, int An,
                   magma_desc_t *B, int Bm, int Bn)
{
    struct starpu_codelet *zlacpy_codelet;
    void (*callback)(void*) = options->profiling ? cl_zlacpy_callback : NULL;
    int lda = BLKLDD( A, Am );
    int ldb = BLKLDD( B, Bm );

#ifdef MORSE_USE_MULTICORE
    zlacpy_codelet = options->parallel ? &cl_zlacpy_mc : &cl_zlacpy;
#else
    zlacpy_codelet = &cl_zlacpy;
#endif
    
    starpu_insert_task(
            zlacpy_codelet,
            STARPU_VALUE,  &uplo,  sizeof(PLASMA_enum),
            STARPU_VALUE,  &m,     sizeof(int),
            STARPU_VALUE,  &n,     sizeof(int), 
            STARPU_R,  BLKADDR( A, PLASMA_Complex64_t, Am, An ),
            STARPU_VALUE,  &lda,   sizeof(int),
            STARPU_W, BLKADDR( B, PLASMA_Complex64_t, Bm, Bn ),
            STARPU_VALUE,  &ldb,   sizeof(int),
            STARPU_PRIORITY, options->priority,
            STARPU_CALLBACK, callback, NULL,
            0);
}
