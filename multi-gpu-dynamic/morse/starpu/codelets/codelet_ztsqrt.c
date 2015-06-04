/**
 *
 *  @file codelet_ztsqrt.c
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
static void cl_ztsqrt_cpu_func(void *descr[], void *cl_arg)
{
    int M;
    int N;
    int IB;
    PLASMA_Complex64_t *A1;
    int LDA1;
    PLASMA_Complex64_t *A2;
    int LDA2;
    PLASMA_Complex64_t *T;
    int LDT;
    PLASMA_Complex64_t *TAU;
    PLASMA_Complex64_t *WORK;

    morse_starpu_ws_t *h_work;

    starpu_unpack_cl_args(cl_arg, &M, &N, &IB, &LDA1, &LDA2, &LDT, 
                          &h_work, NULL);

    A1   = (PLASMA_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[0]);
    A2   = (PLASMA_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[1]);
    T    = (PLASMA_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[2]);

    TAU  = (PLASMA_Complex64_t*)morse_starpu_ws_getlocal(h_work);
    WORK = TAU + min( M, N );

    CORE_ztsqrt(M, N, IB, A1, LDA1, A2, LDA2, T, LDT, TAU, WORK);
}

/*
 * Codelet Multi-cores
 */
#if defined(MORSE_USE_MULTICORE) && 0
static void cl_ztsqrt_mc_func(void *descr[], void *cl_arg)
{
}
#endif

/*
 * Codelet GPU
 */
#if defined(MORSE_USE_CUDA) && 0
static void cl_ztsqrt_cuda_func(void *descr[], void *cl_arg)
{
    int M;
    int N;
    int IB;
    cuDoubleComplex *d_A1;
    int LDA1;
    cuDoubleComplex *d_A2;
    cuDoubleComplex *h_A2;
    int LDA2;
    cuDoubleComplex *h_T, *d_T;
    cuDoubleComplex *h_D, *d_D;
    cuDoubleComplex *h_TAU;
    int LDT;
    cuDoubleComplex *h_WORK, *d_WORK;

    morse_starpu_ws_t *scratch_work, scratch_h_work;
    morse_starpu_ws_t *scratch_h_a, scratch_h_T;
    morse_starpu_ws_t *scratch_h_D, scratch_d_D;
    morse_starpu_ws_t *scratch_tau;

    starpu_unpack_cl_args(cl_arg, &M, &N, &IB, &LDA1, &LDA2, &LDT, 
                          &scratch_tau, &scratch_work,
                          &scratch_h_work, &scratch_h_a, &scratch_h_T, &scratch_h_D, &scratch_d_D);

    d_A1   = (cuDoubleComplex *)STARPU_MATRIX_GET_PTR(descr[0]);
    d_A2   = (cuDoubleComplex *)STARPU_MATRIX_GET_PTR(descr[1]);
    d_T    = (cuDoubleComplex *)STARPU_MATRIX_GET_PTR(descr[2]);
    d_WORK = morse_starpu_ws_getlocal(scratch_work);
    h_WORK = morse_starpu_ws_getlocal(scratch_h_work);
    h_A2   = morse_starpu_ws_getlocal(scratch_h_a);
    h_T    = morse_starpu_ws_getlocal(scratch_h_T);
    h_D    = morse_starpu_ws_getlocal(scratch_h_D);
    d_D    = morse_starpu_ws_getlocal(scratch_d_D);
    h_TAU  = morse_starpu_ws_getlocal(scratch_tau);

    /* Copy A panel */
    memset(h_T,    0, IB*M *sizeof(cuDoubleComplex));
    memset(h_D,    0, IB*M *sizeof(cuDoubleComplex));
    memset(h_TAU,  0, M    *sizeof(cuDoubleComplex));
    memset(h_WORK, 0, 2*M*M*sizeof(cuDoubleComplex));

    cudaMemcpy(h_A2, d_A2, M*IB*sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost);

    /* Could not we just save the entire blocks directly into D and just load
     * the lower parts afterward ? */
    splagma_zsave_d_from_tile(M, IB, d_A1, d_D);

    cudaMemcpy(h_D, d_D, IB*M*sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost);

    int INFO;
    int MxMx2 = M*M*2;
    magma_ztsqrt_gpu( &M, &N, IB,
                        d_A1, d_A2, h_A2, &M,
                        h_D, d_T, &IB, 
                        h_T, &IB,
                        h_TAU, h_WORK, &MxMx2, 
                        d_WORK, &INFO );

    cudaMemcpy(d_D, h_D, IB*M*sizeof(cuDoubleComplex), cudaMemcpyHostToDevice);

    splagma_zload_d_into_tile(M, IB, d_A1, d_D);

    cudaThreadSynchronize();
}
#endif

/*
 * Codelet definition
 */
//CODELETS(ztsqrt, 3, cl_ztsqrt_cpu_func, cl_ztsqrt_cuda_func, cl_ztsqrt_cpu_func)
CODELETS_CPU(ztsqrt, 3, cl_ztsqrt_cpu_func)

/*
 * Wrapper
 */
void MORSE_ztsqrt( MorseOption_t *options, 
                   int m, int n, int ib,
                   magma_desc_t *A1, int A1m, int A1n,
                   magma_desc_t *A2, int A2m, int A2n,
                   magma_desc_t *T,  int Tm,  int Tn)
{
    starpu_codelet *ztsqrt_codelet;
    void (*callback)(void*) = options->profiling ? cl_ztsqrt_callback : NULL;
    int lda1 = BLKLDD( A1, A1m );
    int lda2 = BLKLDD( A2, A2m );
    int ldt  = BLKLDD( T,  Tm  );
    morse_starpu_ws_t *h_work = (morse_starpu_ws_t*)(options->ws_host);
    morse_starpu_ws_t *d_work = (morse_starpu_ws_t*)(options->ws_device);

#ifdef MORSE_USE_MULTICORE
    ztsqrt_codelet = options->parallel ? &cl_ztsqrt_mc : &cl_ztsqrt;
#else
    ztsqrt_codelet = &cl_ztsqrt;
#endif
    
    starpu_Insert_Task(
            &cl_ztsqrt,
            VALUE,  &m,      sizeof(int),
            VALUE,  &n,      sizeof(int),
            VALUE,  &ib,     sizeof(int),
            INOUT,  BLKADDR( A1, PLASMA_Complex64_t, A1m, A1n ),
            VALUE,  &lda1,   sizeof(int),
            INOUT,  BLKADDR( A2, PLASMA_Complex64_t, A2m, A2n ),
            VALUE,  &lda2,   sizeof(int),
            OUTPUT, BLKADDR( T, PLASMA_Complex64_t, Tm, Tn ),
            VALUE,  &ldt,    sizeof(int),
            VALUE,  &h_work, sizeof(morse_starpu_ws_t *),
            VALUE,  &d_work, sizeof(morse_starpu_ws_t *),
            PRIORITY, options->priority,
            CALLBACK, callback, NULL,
            0);
}
