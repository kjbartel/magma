/*
 *  -- MAGMA (version 1.3.0) --
 *     Univ. of Tennessee, Knoxville
 *     Univ. of California, Berkeley
 *     Univ. of Colorado, Denver
 *     November 2012
 *
 * @generated c Wed Nov 14 22:54:14 2012
 *
 **/
// includes, system
#include <stdio.h>
#include <stdlib.h>
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
   -- Testing cposv
*/
int main( int argc, char** argv)
{
    TESTING_CUDA_INIT();

    real_Double_t   gflops, gpu_perf, gpu_time;
    float          Rnorm, Anorm, Xnorm, *work;
    cuFloatComplex c_one     = MAGMA_C_ONE;
    cuFloatComplex c_neg_one = MAGMA_C_NEG_ONE;
    cuFloatComplex *h_A, *h_R, *h_B, *h_X;
    const char  *uplo     = MagmaLowerStr;
    magma_int_t lda, ldb;
    magma_int_t i, info, szeA, szeB;
    magma_int_t ione     = 1;
    magma_int_t N        = 0;
    magma_int_t NRHS     = 100;
    magma_int_t ISEED[4] = {0,0,0,1};
    const int MAXTESTS   = 10;
    magma_int_t size[MAXTESTS] = { 1024, 2048, 3072, 4032, 5184, 6016, 7040, 8064, 9088, 10112 };

    // process command line arguments
    printf( "\nUsage: %s -N <matrix size> -R <right hand sides> [-L|-U]\n", argv[0] );
    printf( "  -N can be repeated up to %d times\n\n", MAXTESTS );
    int ntest = 0;
    for( int i = 1; i < argc; ++i ) {
        if ( strcmp("-N", argv[i]) == 0 && i+1 < argc ) {
            magma_assert( ntest < MAXTESTS, "error: -N repeated more than maximum %d tests\n", MAXTESTS );
            size[ntest] = atoi( argv[++i] );
            magma_assert( size[ntest] > 0, "error: -N %s is invalid; must be > 0.\n", argv[i] );
            N = max( N, size[ntest] );
            ntest++;
        }
        else if ( strcmp("-R", argv[i]) == 0 && i+1 < argc ) {
            NRHS = atoi( argv[++i] );
            magma_assert( NRHS > 0, "error: -R %is is invalid; must be > 0.\n", argv[i] );
        }
        else if ( strcmp("-L", argv[i]) == 0 ) {
            uplo = MagmaLowerStr;
        }
        else if ( strcmp("-U", argv[i]) == 0 ) {
            uplo = MagmaUpperStr;
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
    
    // allocate maximum amount of memory required
    lda = ldb = N;
    TESTING_MALLOC( h_A, cuFloatComplex, lda*N    );
    TESTING_MALLOC( h_R, cuFloatComplex, lda*N    );
    TESTING_MALLOC( h_B, cuFloatComplex, ldb*NRHS );
    TESTING_MALLOC( h_X, cuFloatComplex, ldb*NRHS );
    TESTING_MALLOC( work, float,         N        );

    printf("    N   NRHS   GPU GFlop/s (sec)   ||B - AX|| / ||A||*||X||\n");
    printf("===========================================================\n");
    
    for( i = 0; i < ntest; ++i ) {
        N   = size[i];
        lda = ldb = N;
        gflops = ( FLOPS_CPOTRF( (float)N ) +
                   FLOPS_CPOTRS( (float)N, (float)NRHS ) ) / 1e9;

        /* ====================================================================
           Initialize the matrix
           =================================================================== */
        szeA = lda*N;
        szeB = ldb*NRHS;
        lapackf77_clarnv( &ione, ISEED, &szeA, h_A );
        lapackf77_clarnv( &ione, ISEED, &szeB, h_B );
        magma_chpd( N, h_A, lda );
        
        // copy A to R and B to X; save A and B for residual
        lapackf77_clacpy( MagmaUpperLowerStr, &N, &N,    h_A, &lda, h_R, &lda );
        lapackf77_clacpy( MagmaUpperLowerStr, &N, &NRHS, h_B, &ldb, h_X, &ldb );

        /* ====================================================================
           Performs operation using MAGMA
           =================================================================== */
        gpu_time = magma_wtime();
        magma_cposv( uplo[0], N, NRHS, h_R, lda, h_X, ldb, &info );
        gpu_time = magma_wtime() - gpu_time;
        if (info != 0)
            printf("magma_cpotrf returned error %d.\n", (int) info);

        gpu_perf = gflops / gpu_time;

        /* =====================================================================
           Residual
           =================================================================== */
        Anorm = lapackf77_clange("I", &N, &N,    h_A, &lda, work);
        Xnorm = lapackf77_clange("I", &N, &NRHS, h_X, &ldb, work);

        blasf77_cgemm( MagmaNoTransStr, MagmaNoTransStr, &N, &NRHS, &N,
                       &c_one,     h_A, &lda,
                                   h_X, &ldb,
                       &c_neg_one, h_B, &ldb );
        
        Rnorm = lapackf77_clange("I", &N, &NRHS, h_B, &ldb, work);

        printf( "%5d  %5d   %7.2f (%7.2f)   %8.2e\n",
                (int) N, (int) NRHS, gpu_perf, gpu_time, Rnorm/(Anorm*Xnorm) );
    }

    /* Memory clean up */
    TESTING_FREE( h_A );
    TESTING_FREE( h_R );
    TESTING_FREE( h_B );
    TESTING_FREE( h_X );
    TESTING_FREE( work );

    TESTING_CUDA_FINALIZE();
}
