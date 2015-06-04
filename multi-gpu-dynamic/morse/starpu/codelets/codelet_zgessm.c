/**
 *
 *  @file codelet_zgessm.c
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
static void cl_zgessm_cpu_func(void *descr[], void *cl_arg)
{
    int m;
    int n;
    int k;
    int ib;
    int *IPIV;
    PLASMA_Complex64_t *tinyL;
    int ldtl;
    PLASMA_Complex64_t *L;
    int ldl;
    PLASMA_Complex64_t *A;
    int lda;

    starpu_unpack_cl_args(cl_arg, &m, &n, &k, &ib, &IPIV, &ldtl, &ldl, &lda);

    tinyL = (PLASMA_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[0]);
    L  = (PLASMA_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[1]);
    A  = (PLASMA_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[2]);

    CORE_zgessm(m, n, k, ib, IPIV, L, ldl, A, lda);
}

/*
 * Codelet Multi-cores
 */
#ifdef MORSE_USE_MULTICORE
static void cl_zgessm_mc_func(void *descr[], void *cl_arg)
{
}
#else
#define cl_zgessm_mc_func cl_zgessm_cpu_func
#endif

/*
 * Codelet GPU
 */
#ifdef MORSE_USE_CUDA
static void cl_zgessm_cuda_func(void *descr[], void *cl_arg)
{
    int m;
    int n;
    int k;
    int ib;
    int *IPIV;
    cuDoubleComplex *dtinyL;
    int ldtl;
    cuDoubleComplex *dL;
    int ldl;
    cuDoubleComplex *dA;
    int lda;
    int info;

    starpu_unpack_cl_args(cl_arg, &m, &n, &k, &ib, &IPIV, &ldtl, &ldl, &lda);

    dtinyL = (cuDoubleComplex *)STARPU_MATRIX_GET_PTR(descr[0]);
    dL     = (cuDoubleComplex *)STARPU_MATRIX_GET_PTR(descr[1]);
    dA     = (cuDoubleComplex *)STARPU_MATRIX_GET_PTR(descr[2]);

    dtinyL += ib; /* The kernel is just using the inverted part */

    magma_zgessm_gpu( 
        'C', m, n, k, ib, 
        IPIV, 
        dtinyL, ldtl, 
        dL,  ldl, 
        dA,  lda, 
        &info);

    /* Cedric: same question than getrf ? */
    cudaThreadSynchronize();
}
#endif

/*
 * Codelet definition
 */
CODELETS(zgessm, 3, cl_zgessm_cpu_func, cl_zgessm_cuda_func, cl_zgessm_cpu_func)

/*
 * Wrapper
 */
void MORSE_zgessm( MorseOption_t *options, 
                          int m, int n, int k, int ib,
                          int *IPIV,
                          magma_desc_t *L, int Lm, int Ln,
                          magma_desc_t *D, int Dm, int Dn,
                          magma_desc_t *A, int Am, int An)
{
    starpu_codelet *zgessm_codelet;
    void (*callback)(void*) = options->profiling ? cl_zgessm_callback : NULL;
    int ldl = BLKLDD( L, Lm );
    int ldd = BLKLDD( D, Dm );
    int lda = BLKLDD( A, Am );

#ifdef MORSE_USE_MULTICORE
    zgessm_codelet = options->parallel ? &cl_zgessm_mc : &cl_zgessm;
#else
    zgessm_codelet = &cl_zgessm;
#endif

    starpu_Insert_Task(zgessm_codelet,
                       VALUE,  &m,      sizeof(int),
                       VALUE,  &n,      sizeof(int),
                       VALUE,  &k,      sizeof(int),
                       VALUE,  &ib,     sizeof(int),
                       VALUE,  &IPIV,   sizeof(int*),
                       INPUT,  BLKADDR( L, PLASMA_Complex64_t, Lm, Ln ),
                       VALUE,  &ldl,    sizeof(int),
                       INPUT,  BLKADDR( D, PLASMA_Complex64_t, Dm, Dn ),
                       VALUE,  &ldd,    sizeof(int),
                       INOUT,  BLKADDR( A, PLASMA_Complex64_t, Am, An ),
                       VALUE,  &lda,    sizeof(int),
                       PRIORITY,     options->priority,
                       CALLBACK,     callback, NULL,
                       0);
}
