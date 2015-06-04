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
   -- Testing sposv_gpu
*/
int main( int argc, char** argv)
{
    TESTING_CUDA_INIT();

    real_Double_t   gflops, gpu_perf, gpu_time;
    float          Rnorm, Anorm, Xnorm, *work;
    float c_one     = MAGMA_S_ONE;
    float c_neg_one = MAGMA_S_NEG_ONE;
    float *h_A, *h_B, *h_X;
    float *d_A, *d_B;
    const char  *uplo     = MagmaLowerStr;
    magma_int_t lda, ldb, N;
    magma_int_t ldda, lddb;
    magma_int_t i, info, szeA, szeB;
    magma_int_t ione     = 1;
    magma_int_t NRHS     = 100;
    magma_int_t ISEED[4] = {0,0,0,1};
    const int MAXTESTS   = 10;
    magma_int_t size[MAXTESTS] = { 1024, 2048, 3072, 4032, 5184, 6016, 7040, 8064, 9088, 10112 };

    // process command line arguments
    printf( "\nUsage:\n" );
    printf( "  %s -N <matrix size> -R <right hand sides> [-L|-U]\n", argv[0] );
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
    lddb = ldda = ((N+31)/32)*32;
    
    TESTING_MALLOC( h_A, float, lda*N    );
    TESTING_MALLOC( h_B, float, ldb*NRHS );
    TESTING_MALLOC( h_X, float, ldb*NRHS );
    TESTING_MALLOC( work, float,         N        );

    TESTING_DEVALLOC( d_A, float, ldda*N    );
    TESTING_DEVALLOC( d_B, float, lddb*NRHS );

    printf("    N   NRHS   GPU GFlop/s (sec)   ||B - AX|| / ||A||*||X||\n");
    printf("===========================================================\n");
    
    for( i = 0; i < ntest; ++i ) {
        N   = size[i];
        lda = ldb = N;
        ldda = ((N+31)/32)*32;
        lddb = ldda;
        gflops = ( FLOPS_SPOTRF( (float)N ) +
                   FLOPS_SPOTRS( (float)N, (float)NRHS ) ) / 1e9;

        /* ====================================================================
           Initialize the matrix
           =================================================================== */
        szeA = lda*N;
        szeB = ldb*NRHS;
        lapackf77_slarnv( &ione, ISEED, &szeA, h_A );
        lapackf77_slarnv( &ione, ISEED, &szeB, h_B );
        /* Symmetrize and increase the diagonal */
        {
            magma_int_t i, j;
            for(i=0; i<N; i++) {
                MAGMA_S_SET2REAL( h_A[i*lda+i], ( MAGMA_S_REAL(h_A[i*lda+i]) + 1.*N ) );
                for(j=0; j<i; j++)
                    h_A[i*lda+j] = (h_A[j*lda+i]);
            }
        }
        
        magma_ssetmatrix( N, N,    h_A, N, d_A, ldda );
        magma_ssetmatrix( N, NRHS, h_B, N, d_B, lddb );

        /* ====================================================================
           Performs operation using MAGMA
           =================================================================== */
        gpu_time = magma_wtime();
        magma_sposv_gpu( uplo[0], N, NRHS, d_A, ldda, d_B, lddb, &info );
        gpu_time = magma_wtime() - gpu_time;
        if (info != 0)
            printf("magma_spotrf_gpu returned error %d.\n", (int) info);

        gpu_perf = gflops / gpu_time;

        /* =====================================================================
           Residual
           =================================================================== */
        magma_sgetmatrix( N, NRHS, d_B, lddb, h_X, ldb );
        
        Anorm = lapackf77_slange("I", &N, &N,    h_A, &lda, work);
        Xnorm = lapackf77_slange("I", &N, &NRHS, h_X, &ldb, work);

        blasf77_sgemm( MagmaNoTransStr, MagmaNoTransStr, &N, &NRHS, &N,
                       &c_one,     h_A, &lda,
                                   h_X, &ldb,
                       &c_neg_one, h_B, &ldb );
        
        Rnorm = lapackf77_slange("I", &N, &NRHS, h_B, &ldb, work);

        printf( "%5d  %5d   %7.2f (%7.2f)   %8.2e\n",
                (int) N, (int) NRHS, gpu_perf, gpu_time, Rnorm/(Anorm*Xnorm) );
    }

    /* Memory clean up */
    TESTING_FREE( h_A );
    TESTING_FREE( h_B );
    TESTING_FREE( h_X );
    TESTING_FREE( work );

    TESTING_DEVFREE( d_A );
    TESTING_DEVFREE( d_B );

    TESTING_CUDA_FINALIZE();
}
