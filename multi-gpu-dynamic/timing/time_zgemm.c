/**
 *
 * @precisions normal z -> c d s
 *
 **/
#define _TYPE  PLASMA_Complex64_t
#define _PREC  double
#define _LAMCH LAPACKE_dlamch_work

#define _NAME  "MAGMA_zgemm"
/* See Lawn 41 page 120 */
#define _FMULS FMULS_GEMM(M, N, K)
#define _FADDS FADDS_GEMM(M, N, K)

#include "./timing.c"

static int
RunTest(int *iparam, double *dparam, morse_time_t *t_) 
{
    PLASMA_Complex64_t alpha, beta;
    PASTE_CODE_IPARAM_LOCALS( iparam );
    
    LDB = max(K, iparam[IPARAM_LDB]);
    LDC = max(M, iparam[IPARAM_LDC]);

    /* Allocate Data */
    PASTE_CODE_ALLOCATE_MATRIX( A,      1, PLASMA_Complex64_t, LDA, K   );
    PASTE_CODE_ALLOCATE_MATRIX( B,      1, PLASMA_Complex64_t, LDB, N   );
    PASTE_CODE_ALLOCATE_MATRIX( C,      1, PLASMA_Complex64_t, LDC, N   );
    PASTE_CODE_ALLOCATE_MATRIX( C2, check, PLASMA_Complex64_t, LDC, N   );

    MAGMA_zplrnt( M, K, A, LDA,  453 );
    MAGMA_zplrnt( K, N, B, LDB, 5673 );
    MAGMA_zplrnt( M, N, C, LDC,  740 );

    LAPACKE_zlarnv_work(1, ISEED, 1, &alpha);
    LAPACKE_zlarnv_work(1, ISEED, 1, &beta );

    if (check)
    {
        memcpy(C2, C, LDC*N*sizeof(PLASMA_Complex64_t));
    }

    START_TIMING();
    MAGMA_zgemm( PlasmaNoTrans, PlasmaNoTrans, M, N, K, alpha, A, LDA, B, LDB, beta, C, LDC );
    STOP_TIMING();
    
    /* Check the solution */
    if (check)
    {
        dparam[IPARAM_RES] = zcheck_gemm( PlasmaNoTrans, PlasmaNoTrans, M, N, K,
                                           alpha, A, LDA, B, LDB, beta, C, C2, LDC,
                                           &(dparam[IPARAM_ANORM]), 
                                           &(dparam[IPARAM_BNORM]), 
                                           &(dparam[IPARAM_XNORM]));
        free(C2);
    }

    free( A );
    free( B );
    free( C );

    return 0;
}
