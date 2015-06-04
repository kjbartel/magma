/**
 *
 * @precisions normal z -> c d s
 *
 **/
#define _TYPE  PLASMA_Complex64_t
#define _PREC  double
#define _LAMCH LAPACKE_dlamch_work

#define _NAME  "MAGMA_ztrsm"
/* See Lawn 41 page 120 */
#define _FMULS FMULS_TRSM( PlasmaLeft, N, NRHS )
#define _FADDS FADDS_TRSM( PlasmaLeft, N, NRHS )

#include "./timing.c"

static int
RunTest(int *iparam, double *dparam, morse_time_t *t_)
{
    PLASMA_Complex64_t alpha;
    PASTE_CODE_IPARAM_LOCALS( iparam );
    
    LDA = max( LDA, N );

    /* Allocate Data */
    PASTE_CODE_ALLOCATE_MATRIX( A,      1, PLASMA_Complex64_t, LDA, N   );
    PASTE_CODE_ALLOCATE_MATRIX( B,      1, PLASMA_Complex64_t, LDB, NRHS);
    PASTE_CODE_ALLOCATE_MATRIX( B2, check, PLASMA_Complex64_t, LDB, NRHS);

     /* Initialiaze Data */
    MAGMA_zplgsy( (PLASMA_Complex64_t)N, N, A, LDA, 453 );
    MAGMA_zplrnt( N, NRHS, B, LDB, 5673 );
    LAPACKE_zlarnv_work(1, ISEED, 1, &alpha);
    alpha = 10.; /*alpha * N  /  2.;*/

    if (check)
    {
        memcpy(B2, B, LDB*NRHS*sizeof(PLASMA_Complex64_t));
    }

    START_TIMING();
    MAGMA_ztrsm( PlasmaLeft, PlasmaUpper, PlasmaNoTrans, PlasmaUnit,
                  N, NRHS, alpha, A, LDA, B, LDB );
    STOP_TIMING();

    /* Check the solution */
    if (check)
    {
        dparam[IPARAM_RES] = zcheck_trsm( PlasmaLeft, PlasmaUpper, PlasmaNoTrans, PlasmaUnit, 
                                           N, NRHS,
                                           alpha, A, LDA, B, B2, LDB,
                                           &(dparam[IPARAM_ANORM]), 
                                           &(dparam[IPARAM_BNORM]),
                                           &(dparam[IPARAM_XNORM]));
        free(B2);
    }

    free( A );
    free( B );

    return 0;
}
