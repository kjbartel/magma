/**
 *
 * @precisions normal z -> c d s
 *
 **/
#define _TYPE  PLASMA_Complex64_t
#define _PREC  double
#define _LAMCH LAPACKE_dlamch_work

#define _NAME  "MAGMA_zgels"
/* See Lawn 41 page 120 */
#define _FMULS (FMULS_GEQRF( M, N ) + FMULS_GEQRS( M, N, NRHS ))
#define _FADDS (FADDS_GEQRF( M, N ) + FADDS_GEQRS( M, N, NRHS ))

#include "./timing.c"

static int
RunTest(int *iparam, double *dparam, morse_time_t *t_) 
{
    PLASMA_Complex64_t *T;
    PASTE_CODE_IPARAM_LOCALS( iparam );
    
    if ( M != N ) {
        fprintf(stderr, "This timing works only with M == N\n");
        return -1;
    }

    /* Allocate Data */
    PASTE_CODE_ALLOCATE_MATRIX( A,    1,     PLASMA_Complex64_t, LDA, N   );
    PASTE_CODE_ALLOCATE_MATRIX( x,    1,     PLASMA_Complex64_t, LDB, NRHS);
    PASTE_CODE_ALLOCATE_MATRIX( Acpy, check, PLASMA_Complex64_t, LDA, N   );
    PASTE_CODE_ALLOCATE_MATRIX( b,    check, PLASMA_Complex64_t, LDB, NRHS);

     /* Initialiaze Data */
    MAGMA_zplrnt( M, N,    A, LDA,  453 );
    MAGMA_zplrnt( M, NRHS, x, LDB, 5673 );

    MAGMA_Alloc_Workspace_zgels(M, N, &T);

    /* Save A and b  */
    if (check) {
        LAPACKE_zlacpy_work(LAPACK_COL_MAJOR, 'A', M, N,    A, LDA, Acpy, LDA);
        LAPACKE_zlacpy_work(LAPACK_COL_MAJOR, 'A', M, NRHS, x, LDB, b,    LDB);
    }

    START_TIMING();
    MAGMA_zgels( PlasmaNoTrans, M, N, NRHS, A, LDA, T, x, LDB );
    STOP_TIMING();
    
    /* Check the solution */
    if (check)
    {
        dparam[IPARAM_RES] = zcheck_solution(M, N, NRHS, Acpy, LDA, b, x, LDB,
                                              &(dparam[IPARAM_ANORM]), 
                                              &(dparam[IPARAM_BNORM]), 
                                              &(dparam[IPARAM_XNORM]));
        free(Acpy); free(b);
    }

    free( T );
    free( A );
    free( x );

    return 0;
}
