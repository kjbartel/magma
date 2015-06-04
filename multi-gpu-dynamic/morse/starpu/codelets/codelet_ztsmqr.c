/**
 *
 *  @file codelet_ztsmqr.c
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
static void cl_ztsmqr_cpu_func(void *descr[], void *cl_arg)
{
    int side;
    int trans;
    int M1;
    int M2;
    int N1;
    int N2;
    int K;
    int IB;
    PLASMA_Complex64_t *A1;
    int LDA1;
    PLASMA_Complex64_t *A2;
    int LDA2;
    PLASMA_Complex64_t *V;
    int LDV;
    PLASMA_Complex64_t *T;
    int LDT;
    PLASMA_Complex64_t *WORK;
    int LDWORK;
    morse_starpu_ws_t *h_work;

    starpu_unpack_cl_args(cl_arg, &side, &trans, &M1, &N1, &M2, &N2, &K, &IB, 
                          &LDA1, &LDA2, &LDV, &LDT, &h_work, &LDWORK);

    A1   = (PLASMA_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[0]);
    A2   = (PLASMA_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[1]);
    V    = (PLASMA_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[2]);
    T    = (PLASMA_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[3]);

    WORK = (PLASMA_Complex64_t*)morse_starpu_ws_getlocal(h_work);

    CORE_ztsmqr(side, trans, M1, N1, M2, N2, K, IB, 
                A1, LDA1, A2, LDA2, V, LDV, T, LDT, WORK, LDWORK);
}

/*
 * Codelet Multi-cores
 */
#if defined(MORSE_USE_MULTICORE) && 0
static void cl_ztsmqr_mc_func(void *descr[], void *cl_arg)
{
}
#endif

/*
 * Codelet GPU
 */
#if defined(MORSE_USE_CUDA) && 0
static void cl_ztsmqr_cuda_func(void *descr[], void *cl_arg)
{
    int side;
    int trans;
    int M1;
    int M2;
    int N1;
    int N2;
    int IB;
    int K;
    int i;
    cuDoubleComplex *A1;
    int LDA1;
    cuDoubleComplex *A2;
    int LDA2;
    cuDoubleComplex *V;
    int LDV;
    cuDoubleComplex *T;
    int LDT;
    cuDoubleComplex *WORK;
    int LDWORK;
    morse_starpu_ws_t *scratch_work;

    starpu_unpack_cl_args(cl_arg, &side, &trans, &M1, &N1, &M2, &N2, &K, &IB, 
                          &LDA1, &LDA2, &LDV, &LDT, &LDWORK, &scratch_work);

    A1   = (cuDoubleComplex *)STARPU_MATRIX_GET_PTR(descr[0]);
    A2   = (cuDoubleComplex *)STARPU_MATRIX_GET_PTR(descr[1]);
    V    = (cuDoubleComplex *)STARPU_MATRIX_GET_PTR(descr[2]);
    T    = (cuDoubleComplex *)STARPU_MATRIX_GET_PTR(descr[3]);
    WORK = morse_starpu_ws_getlocal(scratch_work);

   for (i=0;i<M1/IB;i++)
   {
       magma_newzssrfb( M1, M2, &IB, 
                        V + i*M1*IB, &LDV, 
                        T + i*IB*IB, &LDT, 
                        A1 + i*IB, &LDA1, 
                        A2, &LDA2, 
                        WORK, &M2);

   }
   cudaThreadSynchronize();
}
#endif

/*
 * Codelet definition
 */
//CODELETS(ztsmqr, 4, cl_ztsmqr_cpu_func, cl_ztsmqr_cuda_func, cl_ztsmqr_cpu_func)
CODELETS_CPU(ztsmqr, 4, cl_ztsmqr_cpu_func)

/*
 * Wrapper
 */
void MORSE_ztsmqr( MorseOption_t *options, 
                   int side, int trans,
                   int m1, int n1, int m2, int n2, int k, int ib,
                   magma_desc_t *A1, int A1m, int A1n,
                   magma_desc_t *A2, int A2m, int A2n,
                   magma_desc_t *V,  int Vm,  int Vn,
                   magma_desc_t *T,  int Tm,  int Tn)
{
    starpu_codelet *ztsmqr_codelet;
    void (*callback)(void*) = options->profiling ? cl_ztsmqr_callback : NULL;
    int lda1 = BLKLDD( A1, A1m );
    int lda2 = BLKLDD( A2, A2m );
    int ldv  = BLKLDD( V, Vm );
    int ldt  = BLKLDD( T, Tm );
    int ldw  = side == PlasmaLeft ? ib : A1->desc.mb;
    morse_starpu_ws_t *h_work = (morse_starpu_ws_t*)(options->ws_host);

#ifdef MORSE_USE_MULTICORE
    ztsmqr_codelet = options->parallel ? &cl_ztsmqr_mc : &cl_ztsmqr;
#else
    ztsmqr_codelet = &cl_ztsmqr;
#endif
    
    starpu_Insert_Task(
        &cl_ztsmqr,
        VALUE,    &side,   sizeof(PLASMA_enum),
        VALUE,    &trans,  sizeof(PLASMA_enum),
        VALUE,    &m1,     sizeof(int),
        VALUE,    &n1,     sizeof(int), 
        VALUE,    &m2,     sizeof(int),
        VALUE,    &n2,     sizeof(int),
        VALUE,    &k,      sizeof(int),
        VALUE,    &ib,     sizeof(int),
        INOUT,     BLKADDR( A1, PLASMA_Complex64_t, A1m, A1n ),
        VALUE,    &lda1,   sizeof(int),
        INOUT,     BLKADDR( A2, PLASMA_Complex64_t, A2m, A2n ),
        VALUE,    &lda2,   sizeof(int),
        INPUT,     BLKADDR( V, PLASMA_Complex64_t, Vm, Vn ),
        VALUE,    &ldv,    sizeof(int),
        INPUT,     BLKADDR( T, PLASMA_Complex64_t, Tm, Tn ),
        VALUE,    &ldt,    sizeof(int),
        VALUE,    &h_work, sizeof(morse_starpu_ws_t *),
        VALUE,    &ldw,    sizeof(int),
        PRIORITY,  options->priority,
        CALLBACK,  callback, NULL,
        0);
}
