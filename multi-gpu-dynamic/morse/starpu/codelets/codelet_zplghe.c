/**
 *
 *  @file codelet_zplghe.c
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
 *  @precisions normal z -> c
 *
 **/
#include "morse_starpu.h"

/*
 * Codelet CPU
 */
static void cl_zplghe_cpu_func(void *descr[], void *cl_arg)
{
    double bump;
    int m;
    int n;
    PLASMA_Complex64_t *A;
    int lda;
    int bigM;
    int m0;
    int n0;
    unsigned long long int seed;

    starpu_unpack_cl_args( cl_arg, &bump, &m, &n, &lda, &bigM, &m0, &n0, &seed );

    A = (PLASMA_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[0]);

    CORE_zplghe( bump, m, n, A, lda, bigM, m0, n0, seed );
}

/*
 * Codelet Multi-cores
 */
#if defined(MORSE_USE_MULTICORE) && 0
static void cl_zplghe_mc_func(void *descr[], void *cl_arg)
{
}
#endif

/*
 * Codelet GPU
 */
#if (defined MORSE_USE_CUDA) && 0
static void cl_zplghe_cuda_func(void *descr[], void *cl_arg)
{
    cudaThreadSynchronize();
}
#endif

/*
 * Codelet definition
 */
CODELETS_CPU(zplghe, 1, cl_zplghe_cpu_func)

/*
 * Wrapper
 */
void MORSE_zplghe( MorseOption_t *options, 
                   double bump, int m, int n, magma_desc_t *A, int Am, int An,
                   int bigM, int m0, int n0, unsigned long long int seed )
{
    starpu_codelet *zplghe_codelet;
    void (*callback)(void*) = options->profiling ? cl_zplghe_callback : NULL;
    int lda = BLKLDD( A, Am );

#ifdef MORSE_USE_MULTICORE
    zplghe_codelet = options->parallel ? &cl_zplghe_mc : &cl_zplghe;
#else
    zplghe_codelet = &cl_zplghe;
#endif
    
    starpu_Insert_Task(
            &cl_zplghe,
            VALUE,  &bump,  sizeof(double),
            VALUE,  &m,     sizeof(int),
            VALUE,  &n,     sizeof(int), 
            OUTPUT,  BLKADDR( A, PLASMA_Complex64_t, Am, An ),
            VALUE,  &lda,   sizeof(int),
            VALUE,  &bigM,  sizeof(int),
            VALUE,  &m0,    sizeof(int), 
            VALUE,  &n0,    sizeof(int),
            VALUE,  &seed,  sizeof(unsigned long long int), 
            PRIORITY, options->priority,
            CALLBACK, callback, NULL,
            0);
}
