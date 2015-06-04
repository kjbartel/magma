/**
 *
 * @precisions normal z -> c d s
 *
 **/
#define _TYPE  PLASMA_Complex64_t
#define _PREC  double
#define _LAMCH LAPACKE_dlamch_work

#define _NAME  "MAGMA_zgemm_Tile"
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
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descA, 1, PLASMA_Complex64_t, PlasmaComplexDouble, LDA, M, K );
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descB, 1, PLASMA_Complex64_t, PlasmaComplexDouble, LDB, K, N );
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descC, 1, PLASMA_Complex64_t, PlasmaComplexDouble, LDC, M, N );

    /* Initialiaze Data */
    MAGMA_zplrnt_Tile( descA, 5373 );
    MAGMA_zplrnt_Tile( descB, 7672 );
    MAGMA_zplrnt_Tile( descC, 6387 );
    
    LAPACKE_zlarnv_work(1, ISEED, 1, &alpha);
    LAPACKE_zlarnv_work(1, ISEED, 1, &beta);
    
    /* Save C for check */
    PASTE_TILE_TO_LAPACK( descC, C2, check, PLASMA_Complex64_t, LDC, N );

    START_TIMING();
    MAGMA_zgemm_Tile( PlasmaNoTrans, PlasmaNoTrans, alpha, descA, descB, beta, descC );
    STOP_TIMING();
    
    /* Check the solution */
    if (check)
    {
        PASTE_TILE_TO_LAPACK( descA, A, check, PLASMA_Complex64_t, LDA, K );
        PASTE_TILE_TO_LAPACK( descB, B, check, PLASMA_Complex64_t, LDB, N );
        PASTE_TILE_TO_LAPACK( descC, C, check, PLASMA_Complex64_t, LDC, N );

        dparam[IPARAM_RES] = zcheck_gemm( PlasmaNoTrans, PlasmaNoTrans, M, N, K,
                                           alpha, A, LDA, B, LDB, beta, C, C2, LDC,
                                           &(dparam[IPARAM_ANORM]), 
                                           &(dparam[IPARAM_BNORM]), 
                                           &(dparam[IPARAM_XNORM]));
        free(C2);
    }

    PASTE_CODE_FREE_MATRIX( descA );
    PASTE_CODE_FREE_MATRIX( descB );
    PASTE_CODE_FREE_MATRIX( descC );
    return 0;
}
