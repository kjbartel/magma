/*
 *  -- MAGMA (version 1.3.0) --
 *     Univ. of Tennessee, Knoxville
 *     Univ. of California, Berkeley
 *     Univ. of Colorado, Denver
 *     November 2012
 *
 * @generated d Wed Nov 14 22:54:10 2012
 *
 * @author Mark Gates
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
   -- Testing dgeadd
*/
int main( int argc, char** argv)
{
    TESTING_CUDA_INIT();

    real_Double_t   gflops, gpu_perf, gpu_time, cpu_perf, cpu_time;
    double          error;
    double *h_A, *h_B, *d_A, *d_B;
    double alpha = MAGMA_D_MAKE( 3.1415, 2.718 );
    double c_neg_one = MAGMA_D_NEG_ONE;
    double work[1];
    
    /* Matrix size */
    magma_int_t M = 0, N = 0, size, lda, ldda;
    const int MAXTESTS = 10;
    magma_int_t msize[MAXTESTS] = { 1024, 2048, 3072, 4032, 5184, 6016, 7040, 8064, 9088, 10112 };
    magma_int_t nsize[MAXTESTS] = { 1024, 2048, 3072, 4032, 5184, 6016, 7040, 8064, 9088, 10112 };

    magma_int_t i, info;
    magma_int_t ione = 1;
    magma_int_t ISEED[4] = {0,0,0,1};

    // process command line arguments
    printf( "\nUsage: %s -N <m,n> -c -c2 -l\n"
            "  -N  can be repeated up to %d times. If only m is given, then m=n.\n\n",
            argv[0], MAXTESTS );
    
    int ntest = 0;
    for( int i = 1; i < argc; ++i ) {
        if ( strcmp("-N", argv[i]) == 0 && i+1 < argc ) {
            magma_assert( ntest < MAXTESTS, "error: -N repeated more than maximum %d tests\n", MAXTESTS );
            int m, n;
            info = sscanf( argv[++i], "%d,%d", &m, &n );
            if ( info == 2 && m > 0 && n > 0 ) {
                msize[ ntest ] = m;
                nsize[ ntest ] = n;
            }
            else if ( info == 1 && m > 0 ) {
                msize[ ntest ] = m;
                nsize[ ntest ] = m;  // implicitly
            }
            else {
                printf( "error: -N %s is invalid; ensure m > 0, n > 0.\n", argv[i] );
                exit(1);
            }
            M = max( M, msize[ ntest ] );
            N = max( N, nsize[ ntest ] );
            ntest++;
        }
        else if ( strcmp("-M", argv[i]) == 0 ) {
            printf( "-M has been replaced in favor of -N m,n to allow -N to be repeated.\n\n" );
            exit(1);
        }
        else {
            printf( "invalid argument: %s\n", argv[i] );
            exit(1);
        }
    }
    if ( ntest == 0 ) {
        ntest = MAXTESTS;
        M = msize[ntest-1];
        N = nsize[ntest-1];
    }
    
    lda  = M;
    ldda = ((M+31)/32)*32;

    /* Allocate memory for the matrix */
    TESTING_HOSTALLOC( h_A, double, lda *N );
    TESTING_HOSTALLOC( h_B, double, lda *N );
    TESTING_DEVALLOC(  d_A, double, ldda*N );
    TESTING_DEVALLOC(  d_B, double, ldda*N );

    printf("    M     N   CPU GFlop/s (sec)   GPU GFlop/s (sec)   |Bl-Bm|/|Bl|\n");
    printf("=========================================================================\n");
    for( i = 0; i < ntest; ++i ) {
        M = msize[i];
        N = nsize[i];
        lda    = M;
        ldda   = ((M+31)/32)*32;
        size   = lda*N;
        gflops = 2.*M*N / 1e9;

        lapackf77_dlarnv( &ione, ISEED, &size, h_A );
        lapackf77_dlarnv( &ione, ISEED, &size, h_B );
        
        /* ====================================================================
           Performs operation using MAGMA
           =================================================================== */
        magma_dsetmatrix( M, N, h_A, lda, d_A, ldda );
        magma_dsetmatrix( M, N, h_B, lda, d_B, ldda );
        
        gpu_time = magma_wtime();
        magmablas_dgeadd( M, N, alpha, d_A, ldda, d_B, ldda );
        gpu_time = magma_wtime() - gpu_time;
        gpu_perf = gflops / gpu_time;

        /* =====================================================================
           Performs operation using LAPACK
           =================================================================== */
        cpu_time = magma_wtime();
        for( int j = 0; j < N; ++j ) {
            blasf77_daxpy( &M, &alpha, &h_A[j*lda], &ione, &h_B[j*lda], &ione );
        }
        cpu_time = magma_wtime() - cpu_time;
        cpu_perf = gflops / cpu_time;

        /* =====================================================================
           Check result
           =================================================================== */
        magma_dgetmatrix( M, N, d_B, ldda, h_A, lda );
        
        error = lapackf77_dlange( "F", &M, &N, h_B, &lda, work );
        blasf77_daxpy( &size, &c_neg_one, h_A, &ione, h_B, &ione );
        error = lapackf77_dlange( "F", &M, &N, h_B, &lda, work ) / error;
        
        printf("%5d %5d   %7.2f (%7.2f)   %7.2f (%7.2f)   %8.2e\n",
               (int) M, (int) N, cpu_perf, cpu_time, gpu_perf, gpu_time, error );
    }

    /* Memory clean up */
    TESTING_HOSTFREE( h_A );
    TESTING_HOSTFREE( h_B );
    TESTING_DEVFREE(  d_A );
    TESTING_DEVFREE(  d_B );

    TESTING_CUDA_FINALIZE();
    return 0;
}
