/**
 *
 * @precisions normal z -> c d s
 *
 **/
#define _TYPE  PLASMA_Complex64_t
#define _PREC  double
#define _LAMCH LAPACKE_dlamch_work

#define _NAME  "MAGMA_zgels_Tile"
/* See Lawn 41 page 120 */
#define _FMULS (FMULS_GEQRF( M, N ) + FMULS_GEQRS( M, N, NRHS ))
#define _FADDS (FADDS_GEQRF( M, N ) + FADDS_GEQRS( M, N, NRHS ))

#include "./timing.c"

static int
RunTest(int *iparam, double *dparam, morse_time_t *t_) 
{
    magma_desc_t *descT;
    PASTE_CODE_IPARAM_LOCALS( iparam );

    if ( M != N ) {
        fprintf(stderr, "This timing works only with M == N\n");
        return -1;
    }

    /* Allocate Data */
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descA, 1, PLASMA_Complex64_t, PlasmaComplexDouble, LDA, M, N    );
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descB, 1, PLASMA_Complex64_t, PlasmaComplexDouble, LDB, M, NRHS );

    MAGMA_zplrnt_Tile( descA, 5373 );
    MAGMA_zplrnt_Tile( descB,  673 );

    /* Allocate Workspace */
    MAGMA_Alloc_Workspace_zgels_Tile(M, N, &descT);

    /* Save A and B for check */
    PASTE_TILE_TO_LAPACK( descA, A, check, PLASMA_Complex64_t, LDA, N    );
    PASTE_TILE_TO_LAPACK( descB, B, check, PLASMA_Complex64_t, LDB, NRHS );

    /* Do the computations */
    START_TIMING();
    MAGMA_zgels_Tile( PlasmaNoTrans, descA, descT, descB );
    STOP_TIMING();
    
    /* Allocate Workspace */
    MAGMA_Dealloc_Handle_Tile(&descT);

    /* Check the solution */
    if ( check )
      {
        /* Copy solution to X */
        PASTE_TILE_TO_LAPACK( descB, X, 1, PLASMA_Complex64_t, LDB, NRHS );

        dparam[IPARAM_RES] = zcheck_solution(M, N, NRHS, A, LDA, B, X, LDB,
                                              &(dparam[IPARAM_ANORM]), 
                                              &(dparam[IPARAM_BNORM]), 
                                              &(dparam[IPARAM_XNORM]));
        free(A); free(B); free(X);
      }

    PASTE_CODE_FREE_MATRIX( descA );
    PASTE_CODE_FREE_MATRIX( descB );

    return 0;
}
