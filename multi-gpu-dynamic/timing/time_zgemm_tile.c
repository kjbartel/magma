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
#define _FMULS (n * n * n )
#define _FADDS (n * n * n )

#include "./timing.c"

static int
RunTest(int *iparam, double *dparam, real_Double_t *t_) 
{
    PLASMA_Complex64_t *AT, *BT, *CT;
    PLASMA_Complex64_t *A, *B, *C1, *C2;
    PLASMA_Complex64_t alpha, beta;
    magma_desc_t        *descA, *descB, *descC;
    real_Double_t       t;
    int nb, nb2, nt;
    int n       = iparam[TIMING_N];
    int check   = iparam[TIMING_CHECK];
    int lda     = n;
    
    /* Allocate Data */
    nb  = iparam[TIMING_NB];
    nb2 = nb * nb;
    nt  = n / nb + ((n % nb == 0) ? 0 : 1);

    AT = (PLASMA_Complex64_t *)malloc(lda*n*sizeof(PLASMA_Complex64_t));
    BT = (PLASMA_Complex64_t *)malloc(lda*n*sizeof(PLASMA_Complex64_t));
    CT = (PLASMA_Complex64_t *)malloc(lda*n*sizeof(PLASMA_Complex64_t));

    /* Check if unable to allocate memory */
    if ( (!AT) || (!BT) || (!CT) ) {
        printf("Out of Memory \n ");
        exit(0);
    }
    
     /* Initialiaze Data */
    LAPACKE_zlarnv_work(1, ISEED, 1, &alpha);
    LAPACKE_zlarnv_work(1, ISEED, 1, &beta);
    LAPACKE_zlarnv_work(1, ISEED, lda*n, AT);
    LAPACKE_zlarnv_work(1, ISEED, lda*n, BT);
    LAPACKE_zlarnv_work(1, ISEED, lda*n, CT);

    /* Initialize AT and bT for Symmetric Positif Matrix */
    MAGMA_Desc_Create(&descA, AT, PlasmaComplexDouble, nb, nb, nb*nb, n, n, 0, 0, n, n);
    MAGMA_Desc_Create(&descB, BT, PlasmaComplexDouble, nb, nb, nb*nb, n, n, 0, 0, n, n);
    MAGMA_Desc_Create(&descC, CT, PlasmaComplexDouble, nb, nb, nb*nb, n, n, 0, 0, n, n);

    if (check)
      {
          C2 = (PLASMA_Complex64_t *)malloc(n*lda*sizeof(PLASMA_Complex64_t));
          MAGMA_zTile_to_Lapack(descC, (void*)C2, n);
      }

    /* if (iparam[TIMING_BOUND]) */
    /*     starpu_bound_start(iparam[TIMING_BOUNDDEPS],iparam[TIMING_BOUNDDEPSPRIO]); */
    t = -cWtime();
    MAGMA_zgemm_Tile( PlasmaNoTrans, PlasmaNoTrans, alpha, descA, descB, beta, descC );
    t += cWtime();
    /* if (iparam[TIMING_BOUND]) */
    /*     starpu_bound_stop(); */
    *t_ = t;
    
    /* Check the solution */
    if (check)
      {
          A = (PLASMA_Complex64_t *)malloc(n*lda*sizeof(PLASMA_Complex64_t));
          MAGMA_zTile_to_Lapack(descA, (void*)A, n);
          free(AT);

          B = (PLASMA_Complex64_t *)malloc(n*lda*sizeof(PLASMA_Complex64_t));
          MAGMA_zTile_to_Lapack(descB, (void*)B, n);
          free(BT);

          C1 = (PLASMA_Complex64_t *)malloc(n*lda*sizeof(PLASMA_Complex64_t));
          MAGMA_zTile_to_Lapack(descC, (void*)C1, n);
          free(CT);

          dparam[TIMING_RES] = zcheck_gemm( PlasmaNoTrans, PlasmaNoTrans, n, n, n, 
                                            alpha, A, lda, B, lda, beta, C1, C2, lda,
                                            &(dparam[TIMING_ANORM]), &(dparam[TIMING_BNORM]), 
                                            &(dparam[TIMING_XNORM]));
          free(C2);
      }
    else {
        free( AT );
        free( BT );
        free( CT );
    }

    MAGMA_Desc_Destroy(&descA);
    MAGMA_Desc_Destroy(&descB);
    MAGMA_Desc_Destroy(&descC);

    return 0;
}
