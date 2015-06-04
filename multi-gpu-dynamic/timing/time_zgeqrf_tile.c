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
#define _FMULS (n * (2.0 / 3.0 * n + 1.5 ) * n)
#define _FADDS (n * (2.0 / 3.0 * n + 0.5 ) * n)

#include "./timing.c"

static int
RunTest(int *iparam, double *dparam, real_Double_t *t_) 
{
    PLASMA_Complex64_t *A, *AT, *b, *bT, *x, *T;
    magma_desc_t        *descA, *descB, *descT;
    real_Double_t       t;
    int nb, ib, nt;
    int n     = iparam[TIMING_N];
    int nrhs  = iparam[TIMING_NRHS];
    int check = iparam[TIMING_CHECK];
    int lda = n;
    int ldb = n;

    nb = iparam[TIMING_NB];
    ib = iparam[TIMING_IB];
    nt  = n / nb + ((n % nb == 0) ? 0 : 1);

    /* Allocate Data */
    AT  = (PLASMA_Complex64_t *)malloc(lda*n*sizeof(PLASMA_Complex64_t));

    /* Check if unable to allocate memory */
    if ( !AT ){
        printf("Out of Memory \n ");
        exit(0);
    }

    /* Initialiaze Data */
    MAGMA_Desc_Create(&descA, AT, PlasmaComplexDouble, nb, nb, nb*nb, lda, n, 0, 0, n, n);
    LAPACKE_zlarnv_work(1, ISEED, lda*n, AT);

    /* Allocate Workspace */
    /*MAGMA_Alloc_Workspace_zgels_Tile(n, n, &descT);*/
    T = (PLASMA_Complex64_t *)malloc( nt*nt*ib*nb*sizeof(PLASMA_Complex64_t) );
    MAGMA_Desc_Create(&descT, T, PlasmaComplexDouble, ib, nb, ib*nb, nt*ib, nt*nb, 0, 0, nt*ib, nt*nb );

    /* Save AT in lapack layout for check */
    if ( check ) {
        A = (PLASMA_Complex64_t *)malloc(lda*n *sizeof(PLASMA_Complex64_t));
        MAGMA_zTile_to_Lapack(descA, (void*)A, lda);
    }

    /* if (iparam[TIMING_BOUND]) */
    /*     starpu_bound_start(iparam[TIMING_BOUNDDEPS],iparam[TIMING_BOUNDDEPSPRIO]); */
    t = -cWtime();
    MAGMA_zgeqrf_Tile( descA, descT );
    t += cWtime();
    /* if (iparam[TIMING_BOUND]) */
    /*     starpu_bound_stop(); */
    *t_ = t;
    
    /* Check the solution */
    if ( check )
      {
        b  = (PLASMA_Complex64_t *)malloc(ldb*nrhs*sizeof(PLASMA_Complex64_t));
        bT = (PLASMA_Complex64_t *)malloc(ldb*nrhs*sizeof(PLASMA_Complex64_t));
        x  = (PLASMA_Complex64_t *)malloc(ldb*nrhs*sizeof(PLASMA_Complex64_t));

        LAPACKE_zlarnv_work(1, ISEED, ldb*nrhs, bT);
        MAGMA_Desc_Create(&descB, bT, PlasmaComplexDouble, nb, nb, nb*nb, ldb, nrhs, 0, 0, n, nrhs);
        MAGMA_zTile_to_Lapack(descB, (void*)b, ldb);

        PLASMA_Init(1);
        PLASMA_Disable(PLASMA_AUTOTUNING);
        PLASMA_Set(PLASMA_TILE_SIZE,        iparam[TIMING_NB] );
        PLASMA_Set(PLASMA_INNER_BLOCK_SIZE, iparam[TIMING_IB] );
        PLASMA_zgeqrs_Tile( &(descA->desc), &(descT->desc), &(descB->desc) );
        PLASMA_Finalize();

        MAGMA_zTile_to_Lapack(descB, (void*)x, ldb);

        dparam[TIMING_RES] = zcheck_solution(n, n, nrhs, A, lda, b, x, ldb,
                                             &(dparam[TIMING_ANORM]), &(dparam[TIMING_BNORM]), 
                                             &(dparam[TIMING_XNORM]));

        MAGMA_Desc_Destroy(&descB);
        free( A ); 
        free( b ); 
        free( bT ); 
        free( x );
      }

    /* Allocate Workspace */
    /*MAGMA_Dealloc_Handle_Tile(&descT);*/
    MAGMA_Desc_Destroy(&descT);
    free( T );

    MAGMA_Desc_Destroy(&descA);

    free( AT );

    return 0;
}
