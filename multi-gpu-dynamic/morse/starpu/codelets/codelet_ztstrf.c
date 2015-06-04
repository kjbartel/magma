/**
 *
 *  @file codelet_ztstrf.c
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
static void cl_ztstrf_cpu_func(void *descr[], void *cl_arg)
{
    int m;
    int n;
    int ib;
    int nb;
    PLASMA_Complex64_t *U;
    int ldu;
    PLASMA_Complex64_t *A;
    int lda;
    PLASMA_Complex64_t *L;
    int ldl;
    int *IPIV;
    PLASMA_Complex64_t *WORK;
    int ldwork;
    PLASMA_bool check_info;
    morse_starpu_ws_t *h_work;
    morse_starpu_ws_t *d_work; 
    int iinfo;
    int info;

    starpu_unpack_cl_args(cl_arg, &m, &n, &ib, &nb, &ldu, &lda, &ldl, &IPIV, &ldwork,
                          &h_work, &d_work, &check_info, &iinfo);

    /*
     *  hwork => ib*nb
     */
    WORK = morse_starpu_ws_getlocal(h_work);

    U = (PLASMA_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[0]);
    A = (PLASMA_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[1]);
    L = (PLASMA_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[2]);

    CORE_ztstrf(m, n, ib, nb, U, ldu, A, lda, L, ldl, IPIV, WORK, ldwork, &info);

#if defined(MORSE_USE_CUDA) && !defined(WITHOUTTRTRI)
    /*
     * L stores:
     *      L1     L2    L3     ...
     *      L1^-1  L2^-1 L3^-1  ...
     */
    /* Compute L-1 in lower rectangle of L */
    { 
        int i, sb;
        for (i=0; i<n; i+=ib) {
            sb = min( ib, n-i );
            CORE_zlacpy(PlasmaUpperLower, sb, sb, L+(i*ldl), ldl, L+(i*ldl)+ib, ldl );
            
            CORE_ztrtri( PlasmaLower, PlasmaUnit, sb, L+(i*ldl)+ib, ldl, &info );
            if (info != 0 ) {
                fprintf(stderr, "ERROR, trtri returned with info = %d\n", info);
            }          
        }
    }
#endif

    /* if (info != PLASMA_SUCCESS && check_info) */
    /*     plasma_sequence_flush(quark, sequence, request, iinfo + info); */
}

/*
 * Codelet Multi-cores
 */
#ifdef MORSE_USE_MULTICORE
static void cl_ztstrf_mc_func(void *descr[], void *cl_arg)
{
}
#else
#define cl_ztstrf_mc_func cl_ztstrf_cpu_func
#endif

/*
 * Codelet GPU
 */
#ifdef MORSE_USE_CUDA
static void cl_ztstrf_cuda_func(void *descr[], void *cl_arg)
{
    int m;
    int n;
    int ib;
    int nb;
    cuDoubleComplex *hU, *dU;
    int ldu;
    cuDoubleComplex *hA, *dA;
    int lda;
    cuDoubleComplex *hL, *dL;
    int ldl;
    int *ipiv;
    cuDoubleComplex *hw2, *hw, *dw;
    int ldwork;
    PLASMA_bool check_info;
    morse_starpu_ws_t *h_work;
    morse_starpu_ws_t *d_work; 
    int iinfo;
    int info;

    starpu_unpack_cl_args(cl_arg, &m, &n, &ib, &nb, &ldu, &lda, &ldl, &ipiv, &ldwork,
                          &h_work, &d_work, &check_info, &iinfo);

    dU = (cuDoubleComplex *)STARPU_MATRIX_GET_PTR(descr[0]);
    dA = (cuDoubleComplex *)STARPU_MATRIX_GET_PTR(descr[1]);
    dL = (cuDoubleComplex *)STARPU_MATRIX_GET_PTR(descr[2]);

    /*
     *  hwork => 2*nb*(2*ib+2nb)
     *  dwork => 2*ib*nb
     */
    hw2 = (cuDoubleComplex*)morse_starpu_ws_getlocal(h_work);
    dw  = (cuDoubleComplex*)morse_starpu_ws_getlocal(d_work);

    hU = hw2;
    hA = hU + ldu * nb;
    hL = hA + lda * nb; 
    hw = hL + ldl * nb;
    
    /* Download first panel from A and U */
    cublasGetMatrix( m, n,  sizeof(cuDoubleComplex), dU, ldu, hU, ldu );
    cublasGetMatrix( m, ib, sizeof(cuDoubleComplex), dA, lda, hA, lda );

    /* Initialize L to 0 */
    memset(hL, 0, ldl*nb*sizeof(cuDoubleComplex));

    magma_ztstrf_gpu( 'C', m, n, ib, nb,
                      hU, ldu, dU, ldu, 
                      hA, lda, dA, lda, 
                      hL, ldl, dL, ldl,
                      ipiv, 
                      hw, ldwork, dw, lda,
                      &info );

    cudaThreadSynchronize();
}
#endif

/*
 * Codelet definition
 */
CODELETS(ztstrf, 3, cl_ztstrf_cpu_func, cl_ztstrf_cuda_func, cl_ztstrf_cpu_func)

/*
 * Wrapper
 */
void MORSE_ztstrf( MorseOption_t *options, 
                   int m, int n, int ib, int nb,
                   magma_desc_t *U, int Um, int Un,
                   magma_desc_t *A, int Am, int An,
                   magma_desc_t *L, int Lm, int Ln,
                   int  *IPIV,
                   PLASMA_bool check, int iinfo)
{
    starpu_codelet *ztstrf_codelet;
    void (*callback)(void*) = options->profiling ? cl_ztstrf_callback : NULL;
    int ldu = BLKLDD( U, Um );
    int lda = BLKLDD( A, Am );
    int ldl = BLKLDD( L, Lm );
    morse_starpu_ws_t *h_work = (morse_starpu_ws_t*)(options->ws_host);
    morse_starpu_ws_t *d_work = (morse_starpu_ws_t*)(options->ws_device);

#ifdef MORSE_USE_MULTICORE
    ztstrf_codelet = options->parallel ? &cl_ztstrf_mc : &cl_ztstrf;
#else
    ztstrf_codelet = &cl_ztstrf;
#endif

    starpu_Insert_Task(ztstrf_codelet,
                       VALUE,  &m,      sizeof(int),
                       VALUE,  &n,      sizeof(int),
                       VALUE,  &ib,     sizeof(int),
                       VALUE,  &nb,     sizeof(int),
                       INOUT,  BLKADDR( U, PLASMA_Complex64_t, Um, Un ),
                       VALUE,  &ldu,    sizeof(int),
                       INOUT,  BLKADDR( A, PLASMA_Complex64_t, Am, An ),
                       VALUE,  &lda,    sizeof(int),
                       OUTPUT, BLKADDR( L, PLASMA_Complex64_t, Lm, Ln ),
                       VALUE,  &ldl,    sizeof(int),
                       VALUE,  &IPIV,   sizeof(int*),
                       VALUE,  &nb,     sizeof(int),
                       VALUE,  &h_work, sizeof(morse_starpu_ws_t*),
                       VALUE,  &d_work, sizeof(morse_starpu_ws_t*),
                       VALUE,  &check,  sizeof(PLASMA_bool),
                       VALUE,  &iinfo,  sizeof(int),
                       PRIORITY,     options->priority,
                       CALLBACK,     callback, NULL,
                       0);
}
