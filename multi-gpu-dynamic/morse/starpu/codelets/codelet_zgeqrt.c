/**
 *
 *  @file codelet_zgeqrt.c
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
static void cl_zgeqrt_cpu_func(void *descr[], void *cl_arg)
{
    int M;
    int N;
    int IB;
    PLASMA_Complex64_t *A;
    int LDA;
    PLASMA_Complex64_t *T;
    int LDT;
    PLASMA_Complex64_t *TAU;
    PLASMA_Complex64_t *WORK;

    morse_starpu_ws_t *h_work;

    starpu_unpack_cl_args(cl_arg, &M, &N, &IB, &LDA, &LDT,
                          &h_work, NULL);

    /* descr[0] : tile from A, descr[1] : tile from T */
    A = (PLASMA_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[0]);
    T = (PLASMA_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[1]);

    TAU  = (PLASMA_Complex64_t*)morse_starpu_ws_getlocal(h_work);
    WORK = TAU + min( M, N );

    CORE_zgeqrt(M, N, IB, A, LDA, T, LDT, TAU, WORK);
}

/*
 * Codelet Multi-cores
 */
#if defined(MORSE_USE_MULTICORE) && 0
static void cl_zgeqrt_mc_func(void *descr[], void *cl_arg)
{
}
#endif

/*
 * Codelet GPU
 */
#if defined(MORSE_USE_CUDA) && 0
static void cl_zgeqrt_cuda_func(void *descr[], void *cl_arg)
{
    int M;
    int N;
    int IB;
    cuDoubleComplex *h_A, *d_A;
    int LDA;
    cuDoubleComplex *h_T, *d_T;
    int LDT;
    cuDoubleComplex *h_D;
    cuDoubleComplex *h_TAU;
    cuDoubleComplex *h_WORK, *d_WORK;
    cuDoubleComplex *d_D;
    int MxMx2;
    int INFO;

    morse_starpu_ws_t *scratch_work;
    morse_starpu_ws_t *scratch_h_work;
    morse_starpu_ws_t *scratch_h_a;
    morse_starpu_ws_t *scratch_h_T;
    morse_starpu_ws_t *scratch_h_D;
    morse_starpu_ws_t *scratch_d_D;
    morse_starpu_ws_t *scratch_tau;

    starpu_unpack_cl_args(cl_arg, &M, &N, &IB, &LDA, &LDT, 
                          &scratch_tau, &scratch_work, 
                          &scratch_h_work, &scratch_h_a, 
                          &scratch_h_T, &scratch_h_D, &scratch_d_D);

    /* descr[0] : tile from A, descr[1] : tile from T */
    d_A    = (cuDoubleComplex *)STARPU_MATRIX_GET_PTR(descr[0]);
    d_T    = (cuDoubleComplex *)STARPU_MATRIX_GET_PTR(descr[1]);
    d_WORK = morse_starpu_ws_getlocal(scratch_work);
    h_A    = morse_starpu_ws_getlocal(scratch_h_a);
    h_T    = morse_starpu_ws_getlocal(scratch_h_T);
    h_D    = morse_starpu_ws_getlocal(scratch_h_D);
    h_WORK = morse_starpu_ws_getlocal(scratch_h_work);
    h_TAU  = morse_starpu_ws_getlocal(scratch_tau);
    d_D    = morse_starpu_ws_getlocal(scratch_d_D);

    /* TODO are the memset really needed ? */
    memset(h_A,    0,  M*N *sizeof(cuDoubleComplex));
    memset(h_T,    0, IB*N *sizeof(cuDoubleComplex));
    memset(h_D,    0, IB*M *sizeof(cuDoubleComplex));
    memset(h_TAU,  0, M    *sizeof(cuDoubleComplex));
    memset(h_WORK, 0, 2*M*M*sizeof(cuDoubleComplex));

    /* Copy A panel */
    cudaMemcpy(h_A, d_A, M*IB*sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost);

    MxMx2 = M*M*2;
    magma_zgeqrt_gpu(
        &M, &N, IB,
        d_A, &LDA,
        h_A, &LDA,
        d_T, &LDT,
        h_T, &LDT,
        h_D, &IB,
        h_TAU,
        h_WORK, &MxMx2, 
        d_WORK, &INFO);

    cudaMemcpy(d_D, h_D, IB*M*sizeof(cuDoubleComplex), cudaMemcpyHostToDevice);
    splagma_zload_d_into_tile(M, IB, d_A, d_D);
    cudaThreadSynchronize();
}
#endif

/*
 * Codelet definition
 */
//CODELETS(zgeqrt, 2, cl_zgeqrt_cpu_func, cl_zgeqrt_cuda_func, cl_zgeqrt_cpu_func)
CODELETS_CPU(zgeqrt, 2, cl_zgeqrt_cpu_func)

/*
 * Wrapper
 */
void MORSE_zgeqrt( MorseOption_t *options, 
                   int m, int n, int ib, 
                   magma_desc_t *A, int Am, int An,
                   magma_desc_t *T, int Tm, int Tn)
{
    starpu_codelet *zgeqrt_codelet;
    void (*callback)(void*) = options->profiling ? cl_zgeqrt_callback : NULL;
    int lda = BLKLDD( A, Am );
    int ldt = BLKLDD( T, Tm );
    morse_starpu_ws_t *h_work = (morse_starpu_ws_t*)(options->ws_host);
    morse_starpu_ws_t *d_work = (morse_starpu_ws_t*)(options->ws_device);

#ifdef MORSE_USE_MULTICORE
    zgeqrt_codelet = options->parallel ? &cl_zgeqrt_mc : &cl_zgeqrt;
#else
    zgeqrt_codelet = &cl_zgeqrt;
#endif
    
    starpu_Insert_Task(
            &cl_zgeqrt,
            VALUE,  &m,      sizeof(int),
            VALUE,  &n,      sizeof(int),
            VALUE,  &ib,     sizeof(int),
            INOUT,  BLKADDR( A, PLASMA_Complex64_t, Am, An ),
            VALUE,  &lda,    sizeof(int),
            OUTPUT, BLKADDR( T, PLASMA_Complex64_t, Tm, Tn ),
            VALUE,  &ldt,    sizeof(int),
            VALUE,  &h_work, sizeof(morse_starpu_ws_t *),
            VALUE,  &d_work, sizeof(morse_starpu_ws_t *),
            PRIORITY, options->priority,
            CALLBACK, callback, NULL,
            0);
}
