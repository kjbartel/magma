/**
 *
 * @precisions normal z -> c d s
 *
 **/
#define _TYPE  PLASMA_Complex64_t
#define _PREC  double
#define _LAMCH LAPACKE_dlamch_work

#define _NAME  "PLASMA_zpotrf_Tile"
/* See Lawn 41 page 120 */
#define _FMULS (n * (1.0 / 6.0 * n + 0.5) * n)
#define _FADDS (n * (1.0 / 6.0 * n )      * n)

#include "./timing.c"

static int
RunTest(int *iparam, double *dparam, real_Double_t *t_) 
{
    PLASMA_Complex64_t *A, *AT, *b, *bT, *x;
    real_Double_t       t;
    magma_desc_t       *descA = NULL;
    magma_desc_t       *descB = NULL;
    int nb, nt;
    int n     = iparam[TIMING_N];
    int nrhs  = iparam[TIMING_NRHS];
    int check = iparam[TIMING_CHECK];
    int nocpu = iparam[TIMING_NO_CPU];
    int lda = n;
    int ldb = n;
    PLASMA_enum uplo = PlasmaLower;

    int peak_profiling = iparam[TIMING_PEAK];
    int profiling      = iparam[TIMING_PROFILE];

    nb  = iparam[TIMING_NB];
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
    MAGMA_zplghe_Tile((double)n, descA, 51 );

    /* Save AT in lapack layout for check */
    if ( check ) {
        A = (PLASMA_Complex64_t *)malloc(lda*n    *sizeof(PLASMA_Complex64_t));
        MAGMA_zTile_to_Lapack( descA, (void*)A, n);
    }

    if ( profiling | peak_profiling )
        MAGMA_Enable( MAGMA_PROFILING_MODE );

    if (nocpu)
        morse_zlocality_allrestrict( MAGMA_CUDA );
    
    /* PLASMA ZPOSV */
    /* if (iparam[TIMING_BOUND]) */
    /*     starpu_bound_start(iparam[TIMING_BOUNDDEPS],iparam[TIMING_BOUNDDEPSPRIO]); */
    t = -cWtime();
    MAGMA_zpotrf_Tile(uplo, descA);
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

        LAPACKE_zlarnv_work(1, ISEED, ldb*nrhs, bT);
        MAGMA_Desc_Create(&descB, bT, PlasmaComplexDouble, nb, nb, nb*nb, ldb, nrhs, 0, 0, n, nrhs);
        MAGMA_zTile_to_Lapack(descB, (void*)b, n);

        MAGMA_zpotrs_Tile( uplo, descA, descB);
        MAGMA_zTile_to_Lapack(descB, (void*)x, n);

        dparam[TIMING_RES] = zcheck_solution(n, n, nrhs, A, lda, b, x, ldb,
                                             &(dparam[TIMING_ANORM]), &(dparam[TIMING_BNORM]), 
                                             &(dparam[TIMING_XNORM]));
        MAGMA_Desc_Destroy(&descB);
        free( A );
        free( b );
        free( bT );
        free( x );
      }

    MAGMA_Desc_Destroy(&descA);
    free(AT);

    if (peak_profiling) {
        real_Double_t peak = 0;
        /*estimate_zgemm_sustained_peak(&peak);*/
        dparam[TIMING_ESTIMATED_PEAK] = (double)peak;
    }
    
    if (profiling)
    {
        /* Profiling of the scheduler */
        morse_schedprofile_display();
        /* Profile of each kernel */
        morse_zdisplay_allprofile();
    }

    return 0;
}
