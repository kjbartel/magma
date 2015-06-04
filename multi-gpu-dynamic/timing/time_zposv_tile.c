/**
 *
 * @precisions normal z -> c d s
 *
 **/
#define _TYPE  PLASMA_Complex64_t
#define _PREC  double
#define _LAMCH LAPACKE_dlamch_work

#define _NAME  "PLASMA_zposv_Tile"
/* See Lawn 41 page 120 */
#define _FMULS (n * (1.0 / 6.0 * n + nrhs + 0.5) * n)
#define _FADDS (n * (1.0 / 6.0 * n + nrhs )      * n)

#include "./timing.c"

static int
RunTest(int *iparam, double *dparam, real_Double_t *t_)
{
    PLASMA_Complex64_t *A, *AT, *b, *bT, *x;
    real_Double_t       t;
    magma_desc_t        *descA, *descB;
    int nb, nb2, nt;
    int n     = iparam[TIMING_N];
    int nrhs  = iparam[TIMING_NRHS];
    int check = iparam[TIMING_CHECK];
    int lda = n;
    int ldb = n;

    nb  = iparam[TIMING_NB];
    nb2 = nb * nb;
    nt  = n / nb + ((n % nb == 0) ? 0 : 1);
    
    /* Allocate Data */
    AT = (PLASMA_Complex64_t *)malloc(nt*nt*nb2*sizeof(PLASMA_Complex64_t));
    bT = (PLASMA_Complex64_t *)malloc(nt*nb2   *sizeof(PLASMA_Complex64_t));

    /* Check if unable to allocate memory */
    if ( (!AT) || (!bT) ) {
        printf("Out of Memory \n ");
        exit(0);
    }

    /* Initialize AT and bT for Symmetric Positif Matrix */
    MAGMA_Desc_Create(&descA, AT, PlasmaComplexDouble, nb, nb, nb*nb, n, n,    0, 0, n, n);
    MAGMA_Desc_Create(&descB, bT, PlasmaComplexDouble, nb, nb, nb*nb, n, nrhs, 0, 0, n, nrhs);
    MAGMA_zplghe_Tile((double)n, descA, 51 );
    LAPACKE_zlarnv_work(1, ISEED, nt*nb2, bT);

    /* Save AT and bT in lapack layout for check */
    if ( check ) {
        A = (PLASMA_Complex64_t *)malloc(lda*n    *sizeof(PLASMA_Complex64_t));
        b = (PLASMA_Complex64_t *)malloc(ldb*nrhs *sizeof(PLASMA_Complex64_t));
        MAGMA_zTile_to_Lapack(descA, (void*)A, n);
        MAGMA_zTile_to_Lapack(descB, (void*)b, n);
    }

    /* PLASMA ZPOSV */
    /* if (iparam[TIMING_BOUND]) */
    /*     starpu_bound_start(iparam[TIMING_BOUNDDEPS],iparam[TIMING_BOUNDDEPSPRIO]); */
    t = -cWtime();
    MAGMA_zposv_Tile(PlasmaUpper, descA, descB);
    t += cWtime();
    /* if (iparam[TIMING_BOUND]) */
    /*     starpu_bound_stop(); */
    *t_ = t;

    /* Check the solution */
    if (check)
      {
        x = (PLASMA_Complex64_t *)malloc(ldb*nrhs *sizeof(PLASMA_Complex64_t));
        MAGMA_zTile_to_Lapack(descB, (void*)x, n);

        dparam[TIMING_RES] = zcheck_solution(n, n, nrhs, A, lda, b, x, ldb,
                                             &(dparam[TIMING_ANORM]), &(dparam[TIMING_BNORM]), 
                                             &(dparam[TIMING_XNORM]));
        free(A); free(b); free(x);
      }

    MAGMA_Desc_Destroy(&descA);
    MAGMA_Desc_Destroy(&descB);

    free(AT); free(bT);

    return 0;
}
