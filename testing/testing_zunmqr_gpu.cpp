/*
 *  -- MAGMA (version 1.3.0) --
 *     Univ. of Tennessee, Knoxville
 *     Univ. of California, Berkeley
 *     Univ. of Colorado, Denver
 *     November 2012
 *
 * @author Mark Gates
 * @precisions normal z -> c d s
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
   -- Testing zunmqr_gpu
*/
int main( int argc, char** argv )
{
    TESTING_CUDA_INIT();
    
    real_Double_t   gflops, gpu_perf, gpu_time, cpu_perf, cpu_time;    
    double error, work[1];
    cuDoubleComplex c_neg_one = MAGMA_Z_NEG_ONE;
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
    
    magma_int_t nb  = magma_get_zgeqrf_nb( m );
    magma_int_t ldc = m;
    magma_int_t lda = max(m,n);
    ldc = ((ldc+31)/32)*32;
    lda = ((lda+31)/32)*32;
    
    // if k < m, the max workspace may not occur for max m, so check all sizes
    magma_int_t lwork, lwork_max = 0;
    for( int i = 0; i < ntest; ++i ) {
        m = msize[i];
        n = nsize[i];
        k = ksize[i];
        nb = magma_get_zgeqrf_nb( m );
        lwork = max( (m - k + nb)*(n + nb) + n*nb,    // side = left
                     (n - k + nb)*(m + nb) + m*nb );  // side = right
        lwork_max = max( lwork_max, lwork );
    }
    // reset to max size
    m = mmax;
    n = nmax;
    k = kmax;
    nb = magma_get_zgeqrf_nb( m );
    
    // Allocate memory for matrices
    cuDoubleComplex *C, *R, *A, *W, *tau;
    TESTING_MALLOC( C, cuDoubleComplex, ldc*n );
    TESTING_MALLOC( R, cuDoubleComplex, ldc*n );
    TESTING_MALLOC( A, cuDoubleComplex, lda*k );
    TESTING_MALLOC( W, cuDoubleComplex, lwork_max );
    TESTING_MALLOC( tau, cuDoubleComplex, k   );
    
    // min_pk for geqrf( p, k, ... ) where for left p=m, for right p=n
    magma_int_t min_pk = min( max(m,n), k );
    magma_int_t dt_size = ( 2*min_pk + ((k + 31)/32)*32 )*nb;
    cuDoubleComplex *dC, *dA, *dT;
    TESTING_DEVALLOC( dC, cuDoubleComplex, ldc*n );
    TESTING_DEVALLOC( dA, cuDoubleComplex, lda*k );
    TESTING_DEVALLOC( dT, cuDoubleComplex, dt_size );
    
    // test all combinations of input parameters
    const char* side[]   = { MagmaLeftStr,      MagmaRightStr   };
    const char* trans[]  = { MagmaConjTransStr, MagmaNoTransStr };

    printf("    M     N     K  side   trans      CPU GFlop/s (sec)   GPU GFlop/s (sec)   ||R||_F / ||QC||_F\n");
    printf("===============================================================================================\n");
    for( int i = 0; i < ntest; ++i ) {
        for( int iside = 0; iside < 2; ++iside ) {
        for( int itran = 0; itran < 2; ++itran ) {
            m = msize[i];
            n = nsize[i];
            k = ksize[i];
            nb  = magma_get_zgeqrf_nb( m );

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
            
            gflops = FLOPS_ZUNMQR( m, n, k, *side[iside] ) / 1e9;
            
            // C is full, m x n
            size = ldc*n;
            lapackf77_zlarnv( &ione, iseed, &size, C );
            magma_zsetmatrix( m,   n, C, ldc, dC, ldc );
            
            // A is m x k (left) or n x k (right)
            lda = (*side[iside] == 'L' ? m : n);
            size = lda*k;
            lapackf77_zlarnv( &ione, iseed, &size, A );
            
            // compute QR factorization to get Householder vectors in dA, tau, dT
            magma_zsetmatrix( lda, k, A,  lda, dA, lda );
            magma_zgeqrf_gpu( lda, k, dA, lda, tau, dT, &info );
            magma_zgetmatrix( lda, k, dA, lda, A,  lda );
            if ( info != 0 )
                printf("magma_zgeqrf_gpu returned error %d\n", info);
            
            /* =====================================================================
               Performs operation using LAPACK
               =================================================================== */
            cpu_time = magma_wtime();
            lapackf77_zunmqr( side[iside], trans[itran],
                              &m, &n, &k,
                              A, &lda, tau, C, &ldc, W, &lwork_max, &info );
            cpu_time = magma_wtime() - cpu_time;
            cpu_perf = gflops / cpu_time;
            if (info != 0)
                printf("lapackf77_zunmqr returned error %d.\n", (int) info);
            
            /* ====================================================================
               Performs operation using MAGMA
               =================================================================== */
            // query for work size
            lwork = -1;
            magma_zunmqr_gpu( *side[iside], *trans[itran],
                              m, n, k,
                              dA, lda, tau, dC, ldc, W, lwork, dT, nb, &info );
            if (info != 0)
                printf("magma_zunmqr_gpu returned error %d (lwork query).\n", (int) info);
            lwork = (magma_int_t) MAGMA_Z_REAL( W[0] );
            if ( lwork < 0 || lwork > lwork_max )
                printf("invalid lwork %d, lwork_max %d\n", lwork, lwork_max );
            
            magma_device_sync();  // sync needed for L,N and R,T cases
            gpu_time = magma_wtime();
            magma_zunmqr_gpu( *side[iside], *trans[itran],
                              m, n, k,
                              dA, lda, tau, dC, ldc, W, lwork, dT, nb, &info );
            magma_device_sync();
            gpu_time = magma_wtime() - gpu_time;
            gpu_perf = gflops / gpu_time;
            if (info != 0)
                printf("magma_zunmqr_gpu returned error %d.\n", (int) info);
            
            magma_zgetmatrix( m, n, dC, ldc, R, ldc );
            
            /* =====================================================================
               compute relative error |QC_magma - QC_lapack| / |QC_lapack|
               =================================================================== */
            error = lapackf77_zlange( "Fro", &m, &n, C, &ldc, work );
            size = ldc*n;
            blasf77_zaxpy( &size, &c_neg_one, C, &ione, R, &ione );
            error = lapackf77_zlange( "Fro", &m, &n, R, &ldc, work ) / error;
            
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
    
    TESTING_DEVFREE( dC );
    TESTING_DEVFREE( dA );
    TESTING_DEVFREE( dT );
    
    // Shutdown
    TESTING_CUDA_FINALIZE();
    return 0;
}
