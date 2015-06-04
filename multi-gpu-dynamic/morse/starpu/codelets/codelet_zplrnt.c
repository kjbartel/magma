/**
 *
 *  @file codelet_zplrnt.c
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
static void cl_zplrnt_cpu_func(void *descr[], void *cl_arg)
{
    int m;
    int n;
    PLASMA_Complex64_t *A;
    int lda;
    int bigM;
    int m0;
    int n0;
    unsigned long long int seed;

    starpu_codelet_unpack_args( cl_arg, &m, &n, &lda, &bigM, &m0, &n0, &seed );

    A = (PLASMA_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[0]);

    CORE_zplrnt( m, n, A, lda, bigM, m0, n0, seed );
}

/*
 * Codelet Multi-cores
 */
#if defined(MORSE_USE_MULTICORE) && 0
static void cl_zplrnt_mc_func(void *descr[], void *cl_arg)
{
}
#endif

/*
 * Codelet GPU
 */
#if (defined MORSE_USE_CUDA) && 0
static void cl_zplrnt_cuda_func(void *descr[], void *cl_arg)
{
    cudaThreadSynchronize();
}
#endif

/*
 * Codelet definition
 */
CODELETS_CPU(zplrnt, 1, cl_zplrnt_cpu_func)

/*
 * Wrapper
 */
void MORSE_zplrnt( MorseOption_t *options, 
                   int m, int n, magma_desc_t *A, int Am, int An,
                   int bigM, int m0, int n0, unsigned long long int seed )
{
    struct starpu_codelet *zplrnt_codelet;
    void (*callback)(void*) = options->profiling ? cl_zplrnt_callback : NULL;
    int lda = BLKLDD( A, Am );

#ifdef MORSE_USE_MULTICORE
    zplrnt_codelet = options->parallel ? &cl_zplrnt_mc : &cl_zplrnt;
#else
    zplrnt_codelet = &cl_zplrnt;
#endif
    
    starpu_insert_task(
            zplrnt_codelet,
            STARPU_VALUE,  &m,     sizeof(int),
            STARPU_VALUE,  &n,     sizeof(int), 
            STARPU_W,  BLKADDR( A, PLASMA_Complex64_t, Am, An ),
            STARPU_VALUE,  &lda,   sizeof(int),
            STARPU_VALUE,  &bigM,  sizeof(int),
            STARPU_VALUE,  &m0,    sizeof(int), 
            STARPU_VALUE,  &n0,    sizeof(int),
            STARPU_VALUE,  &seed,  sizeof(unsigned long long int), 
            STARPU_PRIORITY, options->priority,
            STARPU_CALLBACK, callback, NULL,
            0);
}
