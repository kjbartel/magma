/*
 *  -- MAGMA (version 1.3.0) --
 *     Univ. of Tennessee, Knoxville
 *     Univ. of California, Berkeley
 *     Univ. of Colorado, Denver
 *     November 2012
 *
 * @author Mark Gates
 * @generated s Wed Nov 14 22:54:23 2012
 *
 **/

// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <cuda_runtime_api.h>
#include <cublas.h>
#include <assert.h>

// includes, project
#include "flops.h"
#include "magma.h"
#include "magma_lapack.h"
#include "testings.h"

/* ////////////////////////////////////////////////////////////////////////////
   -- Testing sormqr
*/
int main( int argc, char** argv )
{
    TESTING_CUDA_INIT();
    
    real_Double_t   gflops, gpu_perf, gpu_time, cpu_perf, cpu_time;    
    float error, work[1];
    float c_neg_one = MAGMA_S_NEG_ONE;
    magma_int_t ione =  1;
    
    /* Matrix size */
    magma_int_t m, n, k;
    const int MAXTESTS = 10;
    magma_int_t msize[MAXTESTS] = { 1024, 2048, 3072, 4032, 5184, 6016, 7040, 8064, 9088, 10112 };
    magma_int_t nsize[MAXTESTS] = { 1024, 2048, 3072, 4032, 5184, 6016, 7040, 8064, 9088, 10112 };
    magma_int_t ksize[MAXTESTS] = { 1024, 2048, 3072, 4032, 5184, 6016, 7040, 8064, 9088, 10112 };
    magma_int_t size;
    
    magma_int_t info;
    magma_int_t iseed[4] = {0,0,0,1};
    
    printf( "Usage: %s -N m,n,k -c\n"
            "    -N can be repeated %d times. m > 0, n > 0, k > 0 is required.\n"
            "    If only m,n is given, then n=k. If only m is given, then m=n=k.\n"
            "    -c or setting $MAGMA_TESTINGS_CHECK runs LAPACK and checks result.\n\n",
            argv[0], MAXTESTS );

    int checkres = (getenv("MAGMA_TESTINGS_CHECK") != NULL);

    int ntest = 0;
    magma_int_t nmax = 0;
    magma_int_t mmax = 0;
    magma_int_t kmax = 0;
    for( int i = 1; i < argc; i++ ) {
        if ( strcmp("-N", argv[i]) == 0 && i+1 < argc ) {
            magma_assert( ntest < MAXTESTS, "error: -N repeated more than maximum %d tests\n", MAXTESTS );
            info = sscanf( argv[++i], "%d,%d,%d", &m, &n, &k );
            if ( info == 3 && m > 0 && n > 0 && k > 0 ) {
                msize[ ntest ] = m;
                nsize[ ntest ] = n;
                ksize[ ntest ] = k;
            }
            else if ( info == 2 && m > 0 && n > 0 ) {
                msize[ ntest ] = m;
                nsize[ ntest ] = n;
                ksize[ ntest ] = n;  // implicitly
            }
            else if ( info == 1 && m > 0 ) {
                msize[ ntest ] = m;
                nsize[ ntest ] = m;  // implicitly
                ksize[ ntest ] = m;  // implicitly
            }
            else {
                printf( "error: -N %s is invalid; ensure m > 0, n > 0, k > 0.\n", argv[i] );
                exit(1);
            }
            mmax = max( mmax, msize[ntest] );
            nmax = max( nmax, nsize[ntest] );
            kmax = max( kmax, ksize[ntest] );
            ntest++;
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
        nmax = nsize[ntest-1];
        mmax = msize[ntest-1];
        kmax = ksize[ntest-1];
    }
    m = mmax;
    n = nmax;
    k = kmax;
    assert( n > 0 && m > 0 && k > 0 );
    
    magma_int_t nb = magma_get_sgeqrf_nb( m );
    magma_int_t ldc = m;
    magma_int_t lda = max(m,n);
    ldc = ((ldc+31)/32)*32;
    lda = ((lda+31)/32)*32;
    
    // Allocate memory for matrices
    float *C, *R, *A, *W, *tau;
    magma_int_t lwork = max( m*nb, n*nb );
    magma_int_t lwork_max = lwork;
    TESTING_MALLOC( C, float, ldc*n );
    TESTING_MALLOC( R, float, ldc*n );
    TESTING_MALLOC( A, float, lda*k );
    TESTING_MALLOC( W, float, lwork_max );
    TESTING_MALLOC( tau, float, k   );
    
    // test all combinations of input parameters
    const char* side[]   = { MagmaLeftStr,      MagmaRightStr   };
    const char* trans[]  = { MagmaTransStr, MagmaNoTransStr };

    printf("    M     N     K  side   trans      CPU GFlop/s (sec)   GPU GFlop/s (sec)   ||R||_F / ||QC||_F\n");
    printf("===============================================================================================\n");
    for( int i = 0; i < ntest; ++i ) {
        for( int iside = 0; iside < 2; ++iside ) {
        for( int itran = 0; itran < 2; ++itran ) {
            m = msize[i];
            n = nsize[i];
            k = ksize[i];
            
            if ( *side[iside] == 'L' && m < k ) {
                printf( "%5d %5d %5d  %-5s  %-9s   skipping because side=left and m < k\n",
                        (int) m, (int) n, (int) k, side[iside], trans[itran] );
                continue;
            }
            if ( *side[iside] == 'R' && n < k ) {
                printf( "%5d %5d %5d  %-5s  %-9s   skipping because side=right and n < k\n",
                        (int) m, (int) n, (int) k, side[iside], trans[itran] );
                continue;
            }
            
            gflops = FLOPS_SORMQR( m, n, k, *side[iside] ) / 1e9;
            
            // C is full, m x n
            size = ldc*n;
            lapackf77_slarnv( &ione, iseed, &size, C );
            lapackf77_slacpy( "Full", &m, &n, C, &ldc, R, &ldc );
            //magma_ssetmatrix( m,   n, C, ldc, dC, ldc );
            
            // A is m x k (left) or n x k (right)
            lda = (*side[iside] == 'L' ? m : n);
            size = lda*k;
            lapackf77_slarnv( &ione, iseed, &size, A );
            
            // compute QR factorization to get Householder vectors in A, tau
            magma_sgeqrf( lda, k, A, lda, tau, W, lwork_max, &info );
            if ( info != 0 )
                printf("magma_sgeqrf returned error %d\n", info);
            
            /* =====================================================================
               Performs operation using LAPACK
               =================================================================== */
            cpu_time = magma_wtime();
            lapackf77_sormqr( side[iside], trans[itran],
                              &m, &n, &k,
                              A, &lda, tau, C, &ldc, W, &lwork_max, &info );
            cpu_time = magma_wtime() - cpu_time;
            cpu_perf = gflops / cpu_time;
            if (info != 0)
                printf("lapackf77_sormqr returned error %d.\n", (int) info);
            
            /* ====================================================================
               Performs operation using MAGMA
               =================================================================== */
            // query for work size
            lwork = -1;
            magma_sormqr( *side[iside], *trans[itran],
                          m, n, k,
                          A, lda, tau, R, ldc, W, lwork, &info );
            if (info != 0)
                printf("magma_sormqr returned error %d (lwork query).\n", (int) info);
            lwork = (magma_int_t) MAGMA_S_REAL( W[0] );
            if ( lwork < 0 || lwork > lwork_max )
                printf("invalid lwork %d, lwork_max %d\n", lwork, lwork_max );
            
            gpu_time = magma_wtime();
            magma_sormqr( *side[iside], *trans[itran],
                          m, n, k,
                          A, lda, tau, R, ldc, W, lwork, &info );
            gpu_time = magma_wtime() - gpu_time;
            gpu_perf = gflops / gpu_time;
            if (info != 0)
                printf("magma_sormqr returned error %d.\n", (int) info);
            
            //magma_sgetmatrix( m, n, dC, ldc, R, ldc );
            
            /* =====================================================================
               compute relative error |QC_magma - QC_lapack| / |QC_lapack|
               =================================================================== */
            error = lapackf77_slange( "Fro", &m, &n, C, &ldc, work );
            size = ldc*n;
            blasf77_saxpy( &size, &c_neg_one, C, &ione, R, &ione );
            error = lapackf77_slange( "Fro", &m, &n, R, &ldc, work ) / error;
            
            printf( "%5d %5d %5d  %-5s  %-9s  %7.2f (%7.2f)   %7.2f (%7.2f)   %8.2e\n",
                    (int) m, (int) n, (int) k, side[iside], trans[itran],
                    cpu_perf, cpu_time, gpu_perf, gpu_time, error );
        }}  // end iside, itran
        printf( "\n" );
    }  // end i
    
    // Memory clean up
    TESTING_FREE( C );
    TESTING_FREE( R );
    TESTING_FREE( A );
    TESTING_FREE( W );
    TESTING_FREE( tau );
    
    // Shutdown
    TESTING_CUDA_FINALIZE();
    return 0;
}
