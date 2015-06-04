/**
 *
 *  @file codelet_zunmqr.c
 *
 *  MAGMA codelets kernel
 *  MAGMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver,
 *  and INRIA Bordeaux Sud-Ouest
 *
 *  @version 1.0.0
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
static void cl_zunmqr_cpu_func(void *descr[], void *cl_arg)
{
    int side;
    int trans;
    int M;
    int N;
    int K;
    int IB;
    PLASMA_Complex64_t *A;
    int LDA;
    PLASMA_Complex64_t *T;
    int LDT;
    PLASMA_Complex64_t *C;
    int LDC;
    PLASMA_Complex64_t *WORK;
    int LDWORK;
    morse_starpu_ws_t *h_work;

    starpu_unpack_cl_args(cl_arg, &side, &trans, &M, &N, &K, &IB,
                          &LDA, &LDT, &LDC, &LDWORK, &h_work, NULL);

    A = (PLASMA_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[0]);
    T = (PLASMA_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[1]);
    C = (PLASMA_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[2]);

    WORK = (PLASMA_Complex64_t*)morse_starpu_ws_getlocal(h_work);

    CORE_zunmqr(side, trans, M, N, K, IB, A, LDA, T, LDT, C, LDC, WORK, LDWORK);
}

/*
 * Codelet Multi-cores
 */
#if defined(MORSE_USE_MULTICORE) && 0
static void cl_zunmqr_mc_func(void *descr[], void *cl_arg)
{
}
#endif

/*
 * Codelet GPU
 */
#if defined(MORSE_USE_CUDA) && 0
static void cl_zunmqr_cuda_func(void *descr[], void *cl_arg)
{
    int side;
    int trans;
    int M;
    int N;
    int K;
    int IB;
    int i, rows;
    cuDoubleComplex *A;
    int LDA;
    cuDoubleComplex *T;
    int LDT;
    cuDoubleComplex *C;
    int LDC;
    cuDoubleComplex *WORK;
    int LDWORK;
    morse_starpu_ws_t *h_work;

    starpu_unpack_cl_args(cl_arg, &side, &trans, &M, &N, &K, &IB,
                          &LDA, &LDT, &LDC, &LDWORK, &h_work);

    A = (cuDoubleComplex *)STARPU_MATRIX_GET_PTR(descr[0]);
    T = (cuDoubleComplex *)STARPU_MATRIX_GET_PTR(descr[1]);
    C = (cuDoubleComplex *)STARPU_MATRIX_GET_PTR(descr[2]);

    WORK = morse_starpu_ws_getlocal(h_work);

/* if CPY is defined, we already created a valid copy of the triangular before */
#if defined(MORSE_USE_CPY)
    cuDoubleComplex *Acpy = A;
#else
    /* TODO optimize that process by using some precomputed version of that
     * lower triangle. We don't need to compute it again and again.  */
    cuDoubleComplex *Acpy;
    cudaMalloc((void **)&Acpy, M*M*sizeof(cuDoubleComplex));

    /* copy A into Acpy, then blank the (strict) upper part, and put 1.0 in the
     * diagonal. XXX This is of course a dummy implementation ! */
    cudaMemcpy(Acpy, A, M*M*sizeof(cuDoubleComplex), cudaMemcpyDeviceToDevice);
    splagma_zsplit_diag_block_gpu(M, IB, Acpy);
#endif

    /* NB side, trans and K are hardcoded */
    for (i=0;i<M/IB;i++) {
        rows=M-i*IB;
        magma_zlarfb(MagmaLeft, MagmaConjTrans, MagmaForward, MagmaColumnwise,
                     rows, M, IB,
                     Acpy + i*(IB*M+IB), LDA,
                     T + i*IB*IB, LDT,
                     C + i*IB, LDC,
                     WORK, LDWORK);
    }

#if !defined(MORSE_USE_CPY)
    cudaFree(Acpy);
#endif

    cudaThreadSynchronize();
}
#endif

/*
 * Codelet definition
 */
//CODELETS(zunmqr, 3, cl_zunmqr_cpu_func, cl_zunmqr_cuda_func, cl_zunmqr_cpu_func)
CODELETS_CPU(zunmqr, 3, cl_zunmqr_cpu_func)

/*
 * Wrapper
 */
void MORSE_zunmqr( MorseOption_t *options,
                   int side, int trans,
                   int m, int n, int k, int ib,
                   magma_desc_t *A, int Am, int An,
                   magma_desc_t *T, int Tm, int Tn,
                   magma_desc_t *C, int Cm, int Cn)
{
    starpu_codelet *zunmqr_codelet;
    void (*callback)(void*) = options->profiling ? cl_zunmqr_callback : NULL;
    int lda = BLKLDD( A, Am );
    int ldt = BLKLDD( T, Tm );
    int ldc = BLKLDD( C, Cm );
    int ldw = A->desc.mb;
    morse_starpu_ws_t *h_work = (morse_starpu_ws_t*)(options->ws_host);
    morse_starpu_ws_t *d_work = (morse_starpu_ws_t*)(options->ws_device);

#ifdef MORSE_USE_MULTICORE
    zunmqr_codelet = options->parallel ? &cl_zunmqr_mc : &cl_zunmqr;
#else
    zunmqr_codelet = &cl_zunmqr;
#endif

    starpu_Insert_Task(
        zunmqr_codelet,
        VALUE,    &side,   sizeof(PLASMA_enum),
        VALUE,    &trans,  sizeof(PLASMA_enum),
        VALUE,    &m,      sizeof(int),
        VALUE,    &n,      sizeof(int),
        VALUE,    &k,      sizeof(int),
        VALUE,    &ib,     sizeof(int),
        INPUT,     BLKADDR( A, PLASMA_Complex64_t, Am, An ),
        VALUE,    &lda,    sizeof(int),
        INPUT,     BLKADDR( T, PLASMA_Complex64_t, Tm, Tn ),
        VALUE,    &ldt,    sizeof(int),
        INOUT,     BLKADDR( C, PLASMA_Complex64_t, Cm, Cn ),
        VALUE,    &ldc,    sizeof(int),
        VALUE,    &ldw,    sizeof(int),
        VALUE,    &h_work, sizeof(morse_starpu_ws_t *),
        VALUE,    &d_work, sizeof(morse_starpu_ws_t *),
        PRIORITY, options->priority,
        CALLBACK, callback, NULL,
        0);
}
