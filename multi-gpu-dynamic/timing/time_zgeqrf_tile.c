/**
 *
 * @precisions normal z -> c d s
 *
 **/
#define _TYPE  PLASMA_Complex64_t
#define _PREC  double
#define _LAMCH LAPACKE_dlamch_work

#define _NAME  "MAGMA_zgeqrf_Tile"
/* See Lawn 41 page 120 */
#define _FMULS FMULS_GEQRF( M, N )
#define _FADDS FADDS_GEQRF( M, N )

#include "./timing.c"

static int
RunTest(int *iparam, double *dparam, morse_time_t *t_) 
{
    magma_desc_t *descT;
    PASTE_CODE_IPARAM_LOCALS( iparam );

    if ( M != N && check ) {
        fprintf(stderr, "Check cannot be perfomed with M != N\n");
        check = 0;
    }

    /* Allocate Data */
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descA, 1, PLASMA_Complex64_t, PlasmaComplexDouble, LDA, M, N );
    MAGMA_zplrnt_Tile( descA, 5373 );

    /* Allocate Workspace */
    MAGMA_Alloc_Workspace_zgels(M, N, &descT);

    /* Save A for check */
    PASTE_TILE_TO_LAPACK( descA, A, check, PLASMA_Complex64_t, LDA, N );
    
    /* Do the computations */
    START_TIMING();
    MAGMA_zgeqrf_Tile( descA, descT );
    STOP_TIMING();
    
    /* Check the solution */
    if ( check )
    {
        /* Allocate B for check */
        PASTE_CODE_ALLOCATE_MATRIX_TILE( descB, 1, PLASMA_Complex64_t, PlasmaComplexDouble, LDB, N, NRHS );
     
        /* Initialize and save B */
        MAGMA_zplrnt_Tile( descB, 2264 );
        PASTE_TILE_TO_LAPACK( descB, B, 1, PLASMA_Complex64_t, LDB, NRHS );

        /* Compute the solution */
#if 0
        MAGMA_zgeqrs_Tile( descA, descT, descB );
#elif 0
        PLASMA_Init(-1);
        PLASMA_Disable(PLASMA_AUTOTUNING);
        PLASMA_Set(PLASMA_TILE_SIZE,        NB );
        PLASMA_Set(PLASMA_INNER_BLOCK_SIZE, IB );
        PLASMA_zgeqrs_Tile( &(descA->desc), &(descT->desc), &(descB->desc) );
        PLASMA_Finalize();
#endif

        /* Copy solution to X */
        PASTE_TILE_TO_LAPACK( descB, X, 1, PLASMA_Complex64_t, LDB, NRHS );

        dparam[IPARAM_RES] = zcheck_solution(M, N, NRHS, A, LDA, B, X, LDB,
                                              &(dparam[IPARAM_ANORM]), 
                                              &(dparam[IPARAM_BNORM]), 
                                              &(dparam[IPARAM_XNORM]));

        /* Free checking structures */
        PASTE_CODE_FREE_MATRIX( descB );
        free( A ); 
        free( B ); 
        free( X );
    }

    /* Free data */
    MAGMA_Dealloc_Workspace( &descT );
    PASTE_CODE_FREE_MATRIX( descA );

    return 0;
}
