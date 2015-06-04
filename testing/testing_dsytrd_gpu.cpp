/*
    -- MAGMA (version 1.2.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       June 2012

       @author Raffaele Solca
       @author Stan Tomov

       @generated d Thu Jun 28 12:31:50 2012

*/

// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <cuda_runtime_api.h>
#include <cublas.h>

// includes, project
#include "flops.h"
#include "magma.h"
#include "magma_lapack.h"
#include "testings.h"

// Flops formula
#define PRECISION_d
#if defined(PRECISION_z) || defined(PRECISION_c)
#define FLOPS(n) ( 6. * FMULS_SYTRD(n) + 2. * FADDS_SYTRD(n))
#else
#define FLOPS(n) (      FMULS_SYTRD(n) +      FADDS_SYTRD(n))
#endif

// This version uses much faster DSYMV (from MAGMA BLAS) but requires extra space 
#define USE_DSYTRD2

/* ////////////////////////////////////////////////////////////////////////////
   -- Testing dsytrd_gpu
*/
int main( int argc, char** argv)
{
    TESTING_CUDA_INIT();

    magma_timestr_t       start, end;
    double           eps, flops, gpu_perf, cpu_perf;
    double *h_A, *h_R, *d_R, *h_Q, *h_work, *work, *dwork;
    double *tau;
    double          *diag, *offdiag, *rwork;
    double           result[2] = {0., 0.};

    /* Matrix size */
    magma_int_t N = 0, n2, lda, lwork;
    magma_int_t size[10] = {1024,2048,3072,4032,5184,6016,7040,8064,9088,10112};

    magma_int_t i, info, nb, checkres, once = 0;
    magma_int_t ione     = 1;
    magma_int_t itwo     = 2;
    magma_int_t ithree   = 3;
    magma_int_t ISEED[4] = {0,0,0,1};
    char *uplo = (char *)MagmaLowerStr;

    if (argc != 1){
        for(i = 1; i<argc; i++){
            if (strcmp("-N", argv[i])==0) {
                N = atoi(argv[++i]);
                once = 1;
            }
            else if (strcmp("-U", argv[i])==0)
                uplo = (char *)MagmaUpperStr;
            else if (strcmp("-L", argv[i])==0)
                uplo = (char *)MagmaLowerStr;
        }
        if ( N > 0 )
            printf("  testing_dsytrd_gpu -L|U -N %d\n\n", (int) N);
        else
        {
            printf("\nUsage: \n");
            printf("  testing_dsytrd_gpu -L|U -N %d\n\n", 1024);
            exit(1);
        }
    }
    else {
        printf("\nUsage: \n");
        printf("  testing_dsytrd_gpu -L|U -N %d\n\n", 1024);
        N = size[9];
    }

    checkres  = getenv("MAGMA_TESTINGS_CHECK") != NULL;

    eps = lapackf77_dlamch( "E" );
    lda = N;
    n2  = lda * N;
    nb  = magma_get_dsytrd_nb(N);
    /* We suppose the magma nb is bigger than lapack nb */
    lwork = N*nb; 

    magma_int_t ldwork = lda*N/16+32;

    /* Allocate host memory for the matrix */
    TESTING_MALLOC(    h_A,    double, lda*N );
    TESTING_DEVALLOC(  d_R,    double, lda*N );
#ifdef USE_DSYTRD2
    TESTING_DEVALLOC(dwork,    double, lda*N );
#endif
    TESTING_HOSTALLOC( h_R,    double, lda*N );
    TESTING_HOSTALLOC( h_work, double, lwork );
    TESTING_MALLOC(    tau,    double, N     );
    TESTING_MALLOC( diag,    double, N   );
    TESTING_MALLOC( offdiag, double, N-1 );

    /* To avoid uninitialized variable warning */
    h_Q   = NULL;
    work  = NULL;
    rwork = NULL; 

    if ( checkres ) {
        TESTING_MALLOC( h_Q,  double, lda*N );
        TESTING_MALLOC( work, double, 2*N*N );
#if defined(PRECISION_z) || defined(PRECISION_c) 
        TESTING_MALLOC( rwork, double, N );
#endif
    }

    printf("  N    CPU GFlop/s    GPU GFlop/s (sec)  |A-QHQ'|/N|A|  |I-QQ'|/N \n");
    printf("===================================================================\n");
    for(i=0; i<10; i++){
        if ( !once ) {
            N = size[i];
        }
        lda  = N;
        n2   = N*lda;
        flops = FLOPS( (double)N ) / 1e6;

        /* ====================================================================
           Initialize the matrix
           =================================================================== */
        lapackf77_dlarnv( &ione, ISEED, &n2, h_A );
        /* Make the matrix hermitian */
        {
            magma_int_t i, j;
            for(i=0; i<N; i++) {
                MAGMA_D_SET2REAL( h_A[i*lda+i], ( MAGMA_D_REAL(h_A[i*lda+i]) ) );
                for(j=0; j<i; j++)
                    h_A[i*lda+j] = (h_A[j*lda+i]);
            }
        }
        magma_dsetmatrix( N, N, h_A, lda, d_R, lda );

        /* ====================================================================
           Performs operation using MAGMA
           =================================================================== */
        start = get_current_time();
#ifdef USE_DSYTRD2
        magma_dsytrd2_gpu(uplo[0], N, d_R, lda, diag, offdiag,
                         tau, h_R, lda, h_work, lwork, dwork , ldwork, &info);
#else
        magma_dsytrd_gpu(uplo[0], N, d_R, lda, diag, offdiag,
                         tau, h_R, lda, h_work, lwork, &info);
#endif
        end = get_current_time();
        if ( info < 0 )
            printf("Argument %d of magma_dsytrd_gpu had an illegal value\n", (int) -info);

        gpu_perf = flops / GetTimerValue(start,end);

        /* =====================================================================
           Check the factorization
           =================================================================== */
        if ( checkres ) {

            magma_dgetmatrix( N, N, d_R, lda, h_R, lda );
            magma_dgetmatrix( N, N, d_R, lda, h_Q, lda );
            lapackf77_dorgtr(uplo, &N, h_Q, &lda, tau, h_work, &lwork, &info);

#if defined(PRECISION_z) || defined(PRECISION_c) 
            lapackf77_dsyt21(&itwo, uplo, &N, &ione, 
                             h_A, &lda, diag, offdiag,
                             h_Q, &lda, h_R, &lda, 
                             tau, work, rwork, &result[0]);

            lapackf77_dsyt21(&ithree, uplo, &N, &ione, 
                             h_A, &lda, diag, offdiag,
                             h_Q, &lda, h_R, &lda, 
                             tau, work, rwork, &result[1]);

#else

            lapackf77_dsyt21(&itwo, uplo, &N, &ione, 
                             h_A, &lda, diag, offdiag,
                             h_Q, &lda, h_R, &lda, 
                             tau, work, &result[0]);

            lapackf77_dsyt21(&ithree, uplo, &N, &ione, 
                             h_A, &lda, diag, offdiag,
                             h_Q, &lda, h_R, &lda, 
                             tau, work, &result[1]);

#endif
        }

        /* =====================================================================
           Performs operation using LAPACK
           =================================================================== */
        start = get_current_time();
        lapackf77_dsytrd(uplo, &N, h_A, &lda, diag, offdiag, tau, 
                         h_work, &lwork, &info);
        end = get_current_time();

        if (info < 0)
            printf("Argument %d of lapackf77_dsytrd had an illegal value.\n", (int) -info);

        cpu_perf = flops / GetTimerValue(start,end);

        /* =====================================================================
           Print performance and error.
           =================================================================== */
        if ( checkres ) {
            printf("%5d   %6.2f        %6.2f (%6.2f)     %e %e\n",
                   (int) N, cpu_perf, gpu_perf, flops/(1e3*gpu_perf),
                   result[0]*eps, result[1]*eps );
        } else {
            printf("%5d   %6.2f        %6.2f (%6.2f)\n",
                   (int) N, cpu_perf, gpu_perf, flops/(1e3*gpu_perf));
        }
 
        if ( once )
            break;
    }

    /* Memory clean up */
    TESTING_FREE( h_A );
    TESTING_FREE( tau );
    TESTING_FREE( diag );
    TESTING_FREE( offdiag );
    TESTING_HOSTFREE( h_R );
    TESTING_HOSTFREE( h_work );
    TESTING_DEVFREE ( d_R ); 
#ifdef USE_DSYTRD2 
    TESTING_DEVFREE ( dwork );
#endif

    if ( checkres ) {
        TESTING_FREE( h_Q );
        TESTING_FREE( work );
#if defined(PRECISION_z) || defined(PRECISION_c) 
        TESTING_FREE( rwork );
#endif
    }

    /* Shutdown */
    TESTING_CUDA_FINALIZE();
    return EXIT_SUCCESS;
}
