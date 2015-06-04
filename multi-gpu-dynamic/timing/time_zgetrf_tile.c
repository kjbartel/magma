/**
 *
 * @precisions normal z -> c d s
 *
 **/
#define _TYPE  PLASMA_Complex64_t
#define _PREC  double
#define _LAMCH LAPACKE_dlamch_work

#define _NAME  "PLASMA_zgetrf_Tile"
/* See Lawn 41 page 120 */
#define _FMULS (n * (1.0 / 3.0 * n )      * n)
#define _FADDS (n * (1.0 / 3.0 * n - 0.5) * n)

#include "./timing.c"

static int
RunTest(int *iparam, double *dparam, real_Double_t *t_) 
{
    PLASMA_Complex64_t *A, *AT, *b, *bT, *x, *L;
    magma_desc_t        *descA, *descB, *descL;
    real_Double_t       t;
    int                *piv;
    int nb, nb2, nt, ib;
    int n     = iparam[TIMING_N];
    int nrhs  = iparam[TIMING_NRHS];
    int check = iparam[TIMING_CHECK];
    int nocpu = iparam[TIMING_NO_CPU];
    int lda      = n;
    int ldb      = n;

    int peak_profiling = iparam[TIMING_PEAK];
    int profiling      = iparam[TIMING_PROFILE];

    nb  = iparam[TIMING_NB];
    ib  = iparam[TIMING_IB];
    nb2 = nb * nb;
    nt  = n / nb + ((n % nb == 0) ? 0 : 1);
    
    /* Allocate Data */
    AT = (PLASMA_Complex64_t *)malloc(lda*n*sizeof(PLASMA_Complex64_t));

    /* Check if unable to allocate memory */
    if ( !AT ){
        printf("Out of Memory \n ");
        exit(0);
    }

    /* Initialiaze Data */
    MAGMA_Desc_Create(&descA, AT, PlasmaComplexDouble, nb, nb, nb*nb, lda, n, 0, 0, n, n);
    MAGMA_zplrnt_Tile(descA, 51 );

    /* Allocate Workspace */
    //MAGMA_Alloc_Workspace_zgesv_Tile(n, &descL, &piv);
    piv = (int *)malloc ( nt*nt*nb*sizeof(int) );
    L   = (PLASMA_Complex64_t *)malloc( nt*nt*2*ib*nb*sizeof(PLASMA_Complex64_t) );
    MAGMA_Desc_Create(&descL, L, PlasmaComplexDouble, 2*ib, nb, 2*ib*nb, nt*2*ib, nt*nb, 0, 0, nt*2*ib, nt*nb );

    /* Save AT in lapack layout for check */
    if ( check ) {
        A = (PLASMA_Complex64_t *)malloc(lda*n*sizeof(PLASMA_Complex64_t));
        MAGMA_zTile_to_Lapack(descA, (void*)A, lda);
    }

    if ( profiling | peak_profiling )
        MAGMA_Enable( MAGMA_PROFILING_MODE );

    if (nocpu)
        morse_zlocality_allrestrict( MAGMA_CUDA );
    
    /* if (iparam[TIMING_BOUND]) */
    /*     starpu_bound_start(iparam[TIMING_BOUNDDEPS],iparam[TIMING_BOUNDDEPSPRIO]); */
    t = -cWtime();
    MAGMA_zgetrf_Tile( descA, descL, piv );
    t += cWtime();
    /* if (iparam[TIMING_BOUND]) */
    /*     starpu_bound_stop(); */
    *t_ = t;
    
    if (nocpu)
        morse_zlocality_allrestore();

    if ( profiling | peak_profiling )
        MAGMA_Disable( MAGMA_PROFILING_MODE );

    /* Check the solution */
    if ( check )
      {
        b  = (PLASMA_Complex64_t *)malloc(ldb*nrhs*sizeof(PLASMA_Complex64_t));
        bT = (PLASMA_Complex64_t *)malloc(ldb*nrhs*sizeof(PLASMA_Complex64_t));
        x  = (PLASMA_Complex64_t *)malloc(ldb*nrhs*sizeof(PLASMA_Complex64_t));

        MAGMA_Desc_Create(&descB, bT, PlasmaComplexDouble, nb, nb, nb*nb, n, nrhs, 0, 0, n, nrhs);
        MAGMA_zplrnt_Tile(descB, 51 );
        MAGMA_zLapack_to_Tile((void*)b, n, descB);

        MAGMA_zgetrs_Tile( descA, descL, piv, descB );

        MAGMA_zTile_to_Lapack(descB, (void*)x, n);

        dparam[TIMING_RES] = zcheck_solution(n, n, nrhs, A, lda, b, x, ldb,
                                             &(dparam[TIMING_ANORM]), &(dparam[TIMING_BNORM]), 
                                             &(dparam[TIMING_XNORM]));

        MAGMA_Desc_Destroy(&descB);
        free( A ); free( b ); free( bT ); free( x );
      }

    /* Deallocate Workspace */
    //MAGMA_Dealloc_Handle_Tile(&descL);

    MAGMA_Desc_Destroy(&descA);
    MAGMA_Desc_Destroy(&descL);
    free( L );
    free( AT );
    free( piv );

    if (profiling)
    {
        /* Profiling of the scheduler */
        morse_schedprofile_display();
        /* Profile of each kernel */
        morse_zdisplay_allprofile();
    }

    return 0;
}
