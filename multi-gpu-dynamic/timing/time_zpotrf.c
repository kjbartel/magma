/**
 *
 * @precisions normal z -> c d s
 *
 **/
#define _TYPE  PLASMA_Complex64_t
#define _PREC  double
#define _LAMCH LAPACKE_dlamch_work

#define _NAME  "PLASMA_zpotrf"
/* See Lawn 41 page 120 */
#define _FMULS (n * (1.0 / 6.0 * n + 0.5) * n)
#define _FADDS (n * (1.0 / 6.0 * n )      * n)

#include "./timing.c"

static int
RunTest(int *iparam, double *dparam, real_Double_t *t_) 
{
    PLASMA_Complex64_t *A, *Acpy, *b, *x;
    real_Double_t       t;
    int n     = iparam[TIMING_N];
    int nrhs  = iparam[TIMING_NRHS];
    int check = iparam[TIMING_CHECK];
    int lda = n;
    int ldb = n;
    int uplo = PlasmaLower;

    /* Allocate Data */
    A = (PLASMA_Complex64_t *)malloc(lda*n*   sizeof(PLASMA_Complex64_t));

    /* Check if unable to allocate memory */
    if ( !A ) {
        printf("Out of Memory \n ");
        exit(0);
    }
    
    /* Initialiaze Data */
    MAGMA_zplghe( (double)n, n, A, lda, 51 );

    /* Save A and b  */
    if (check) {
        Acpy = (PLASMA_Complex64_t *)malloc(lda*n*sizeof(PLASMA_Complex64_t));
        LAPACKE_zlacpy_work(LAPACK_COL_MAJOR,' ', n, n, A, lda, Acpy, lda);
    }

    /* PLASMA ZPOSV */
    /* if (iparam[TIMING_BOUND]) */
    /*     starpu_bound_start(iparam[TIMING_BOUNDDEPS],iparam[TIMING_BOUNDDEPSPRIO]); */
    t = -cWtime();
    MAGMA_zpotrf(uplo, n, A, lda);
    t += cWtime();
    /* if (iparam[TIMING_BOUND]) */
    /*     starpu_bound_stop(); */
    *t_ = t;

    /* Check the solution */
    if (check)
      {
        b = (PLASMA_Complex64_t *)malloc(ldb*nrhs*sizeof(PLASMA_Complex64_t));
        x = (PLASMA_Complex64_t *)malloc(ldb*nrhs*sizeof(PLASMA_Complex64_t));
        LAPACKE_zlarnv_work(1, ISEED, n*nrhs, x);
        LAPACKE_zlacpy_work(LAPACK_COL_MAJOR, 'A', n, nrhs, x, ldb, b, ldb);

        MAGMA_zpotrs(uplo, n, nrhs, A, lda, x, ldb);

        dparam[TIMING_RES] = zcheck_solution(n, n, nrhs, Acpy, lda, b, x, ldb,
                                             &(dparam[TIMING_ANORM]), 
                                             &(dparam[TIMING_BNORM]),
                                             &(dparam[TIMING_XNORM]));

        free(Acpy); free(b); free(x);
      }

    free(A);

    return 0;
}
