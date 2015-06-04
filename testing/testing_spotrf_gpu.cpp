/*
 *  -- MAGMA (version 1.3.0) --
 *     Univ. of Tennessee, Knoxville
 *     Univ. of California, Berkeley
 *     Univ. of Colorado, Denver
 *     November 2012
 *
 * @generated s Wed Nov 14 22:54:13 2012
 *
 **/
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

/* ////////////////////////////////////////////////////////////////////////////
   -- Testing spotrf
*/
int main( int argc, char** argv)
{
    TESTING_CUDA_INIT();

    real_Double_t   gflops, gpu_perf, gpu_time, cpu_perf, cpu_time;
    float *h_A, *h_R;
    float *d_A;
    
    /* Matrix size */
    magma_int_t N = 0, n2, lda, ldda;
    const int MAXTESTS = 10;
    magma_int_t size[MAXTESTS] = { 1024, 2048, 3072, 4032, 5184, 6016, 7040, 8064, 9088, 10112 };
    
    magma_int_t i, info;
    const char *uplo     = MagmaLowerStr;
    float c_neg_one = MAGMA_S_NEG_ONE;
    magma_int_t ione     = 1;
    magma_int_t ISEED[4] = {0,0,0,1};
    float      work[1], error;
    magma_int_t checkres;

    checkres = getenv("MAGMA_TESTINGS_CHECK") != NULL;
    
    // process command line arguments
    printf( "\nUsage: %s -N <n> [-L|-U] -c\n", argv[0] );
    printf( "  -N can be repeated up to %d times.\n", MAXTESTS );
    printf( "  -c or setting $MAGMA_TESTINGS_CHECK runs LAPACK and checks result.\n\n" );
    int ntest = 0;
    for( int i = 1; i < argc; ++i ) {
        if ( strcmp("-N", argv[i]) == 0 && i+1 < argc ) {
            magma_assert( ntest < MAXTESTS, "error: -N repeated more than maximum %d tests\n", MAXTESTS );
            size[ ntest ] = atoi( argv[++i] );
            magma_assert( size[ ntest ] > 0, "error: -N %s is invalid; must be > 0.\n", argv[i] );
            N = max( N, size[ ntest ] );
            ntest++;
        }
        else if ( strcmp("-L", argv[i]) == 0 ) {
            uplo = MagmaLowerStr;
        }
        else if ( strcmp("-U", argv[i]) == 0 ) {
            uplo = MagmaUpperStr;
        }
        else if ( strcmp("-c", argv[i]) == 0 ) {
            checkres = true;
        }
        else {
            printf( "invalid argument: %s\n", argv[i] );
            exit(1);
        }
    }
    if ( ntest == 0 ) {
        ntest = MAXTESTS;
        N = size[ntest-1];
    }

    /* Allocate memory for the matrix */
    n2   = N*N;
    ldda = ((N+31)/32) * 32;
    TESTING_MALLOC(    h_A, float, n2);
    TESTING_HOSTALLOC( h_R, float, n2);
    TESTING_DEVALLOC(  d_A, float, ldda*N );

    printf("  N     CPU GFlop/s (sec)   GPU GFlop/s (sec)   ||R_magma - R_lapack||_F / ||R_lapack||_F\n");
    printf("========================================================\n");
    for( i = 0; i < ntest; ++i ) {
        N   = size[i];
        lda = N;
        n2  = lda*N;
        ldda = ((N+31)/32)*32;
        gflops = FLOPS_SPOTRF( N ) / 1e9;

        /* Initialize the matrix */
        lapackf77_slarnv( &ione, ISEED, &n2, h_A );
        magma_shpd( N, h_A, lda );
        lapackf77_slacpy( MagmaUpperLowerStr, &N, &N, h_A, &lda, h_R, &lda );
        magma_ssetmatrix( N, N, h_A, lda, d_A, ldda );
    
        /* ====================================================================
           Performs operation using MAGMA
           =================================================================== */
        gpu_time = magma_wtime();
        magma_spotrf_gpu(uplo[0], N, d_A, ldda, &info);
        gpu_time = magma_wtime() - gpu_time;
        gpu_perf = gflops / gpu_time;
        if (info != 0)
            printf("magma_spotrf_gpu returned error %d.\n", (int) info);
        
        if ( checkres ) {
            /* =====================================================================
               Performs operation using LAPACK
               =================================================================== */
            cpu_time = magma_wtime();
            lapackf77_spotrf(uplo, &N, h_A, &lda, &info);
            cpu_time = magma_wtime() - cpu_time;
            cpu_perf = gflops / cpu_time;
            if (info != 0)
                printf("lapackf77_spotrf returned error %d.\n", (int) info);
            
            /* =====================================================================
               Check the result compared to LAPACK
               =================================================================== */
            magma_sgetmatrix( N, N, d_A, ldda, h_R, lda );
            error = lapackf77_slange("f", &N, &N, h_A, &lda, work);
            blasf77_saxpy(&n2, &c_neg_one, h_A, &ione, h_R, &ione);
            error = lapackf77_slange("f", &N, &N, h_R, &lda, work) / error;
            
            printf("%5d   %7.2f (%7.2f)   %7.2f (%7.2f)   %8.2e\n",
                   (int) N, cpu_perf, cpu_time, gpu_perf, gpu_time, error );
        }
        else {
            printf("%5d     ---   (  ---  )   %7.2f (%7.2f)     ---  \n",
                   (int) N, gpu_perf, gpu_time );
        }
    }

    /* Memory clean up */
    TESTING_FREE( h_A );
    TESTING_HOSTFREE( h_R );
    TESTING_DEVFREE( d_A );

    TESTING_CUDA_FINALIZE();
    return 0;
}
