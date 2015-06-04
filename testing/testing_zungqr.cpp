/*
    -- MAGMA (version 1.2.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       June 2012

       @precisions normal z -> s d c

       @author Stan Tomov
       @author Mathieu Faverge
       @author Mark Gates
*/

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
   -- Testing zungqr
*/
int main( int argc, char** argv )
{
    TESTING_CUDA_INIT();

    real_Double_t    gflops, gpu_perf=0., cpu_perf=0., gpu_time=0., cpu_time=0.;
    double           error=0., work[1];
    cuDoubleComplex  c_neg_one = MAGMA_Z_NEG_ONE;
    cuDoubleComplex *hA, *hR, *tau, *hwork;
    cuDoubleComplex *dA, *dT;

    /* Matrix size */
    magma_int_t m=0, n=0, k=0;
    magma_int_t n2, lda, ldda, lwork, min_mn, nb;
    const int MAXTESTS = 10;
    magma_int_t msize[MAXTESTS] = { 1024, 2048, 3072, 4032, 5184, 6016, 7040, 8064, 9088, 9984 };
    magma_int_t nsize[MAXTESTS] = { 1024, 2048, 3072, 4032, 5184, 6016, 7040, 8064, 9088, 9984 };
    magma_int_t ksize[MAXTESTS] = { 1024, 2048, 3072, 4032, 5184, 6016, 7040, 8064, 9088, 9984 };
    
    magma_int_t info;
    magma_int_t ione     = 1;
    magma_int_t iseed[4] = {0,0,0,1};
    
    printf( "Usage: %s -N m,n,k -c\n"
            "    -N can be repeated %d times. m >= n >= k is required.\n"
            "    If only m,n is given, then n=k. If only m is given, then m=n=k.\n"
            "    -c or setting $MAGMA_TESTINGS_CHECK runs LAPACK and checks result.\n",
            argv[0], MAXTESTS );

    int checkres = (getenv("MAGMA_TESTINGS_CHECK") != NULL);

    int ntest = 0;
    magma_int_t nmax = 0;
    magma_int_t mmax = 0;
    magma_int_t kmax = 0;
    for( int i = 1; i < argc; i++ ) {
        if ( strcmp("-N", argv[i]) == 0 ) {
            if ( ++i >= argc ) {
                printf( "error: -N requires an argument\n" );
                exit(1);
            }
            else if ( ntest == MAXTESTS ) {
                printf( "error: -N exceeded maximum %d tests\n", MAXTESTS );
                exit(1);
            }
            else {
                info = sscanf( argv[i], "%d,%d,%d", &m, &n, &k );
                if ( info == 3 and m >= n and n >= k and k > 0 ) {
                    msize[ ntest ] = m;
                    nsize[ ntest ] = n;
                    ksize[ ntest ] = k;
                }
                else if ( info == 2 and m >= n and n > 0 ) {
                    msize[ ntest ] = m;
                    nsize[ ntest ] = n;
                    ksize[ ntest ] = n;  // implicitly
                }
                else if ( info == 1 and m > 0 ) {
                    msize[ ntest ] = m;
                    nsize[ ntest ] = m;  // implicitly
                    ksize[ ntest ] = m;  // implicitly
                }
                else {
                    printf( "error: -N %s is invalid; ensure m >= n >= k.\n", argv[i] );
                    exit(1);
                }
                mmax = max( mmax, msize[ntest] );
                nmax = max( nmax, nsize[ntest] );
                kmax = max( kmax, ksize[ntest] );
                ntest++;
            }
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
    assert( nmax > 0 and mmax > 0 and kmax > 0 );
    
    // allocate memory for largest problem
    lda    = mmax;
    ldda   = ((mmax + 31)/32)*32;
    n2     = lda * nmax;
    min_mn = min(mmax, nmax);
    nb     = magma_get_zgeqrf_nb( mmax );
    lwork  = (mmax + 2*nmax+nb)*nb;

    TESTING_HOSTALLOC( hA,    cuDoubleComplex, lda*nmax  );
    TESTING_HOSTALLOC( hwork, cuDoubleComplex, lwork     );
    TESTING_MALLOC(    hR,    cuDoubleComplex, lda*nmax  );
    TESTING_MALLOC(    tau,   cuDoubleComplex, min_mn    );
    TESTING_DEVALLOC(  dA,    cuDoubleComplex, ldda*nmax );
    TESTING_DEVALLOC(  dT,    cuDoubleComplex, ( 2*min_mn + ((nmax + 31)/32)*32 )*nb );
    
    printf("\n");
    printf("    m     n     k   CPU GFlop/s (sec)   GPU GFlop/s (sec)   ||R|| / ||A||\n");
    printf("=========================================================================\n");
    for( int i = 0; i < ntest; ++i ){
        m = msize[i];
        n = nsize[i];
        k = ksize[i];
        assert( m >= n and n >= k );
        
        lda  = m;
        ldda = ((m + 31)/32)*32;
        n2 = lda*n;
        nb = magma_get_zgeqrf_nb( m );
        gflops = FLOPS_ZUNGQR( (double)m, (double)n, (double)k ) / 1e9;

        lapackf77_zlarnv( &ione, iseed, &n2, hA );
        lapackf77_zlacpy( MagmaUpperLowerStr, &m, &n, hA, &lda, hR, &lda );
        
        /* ====================================================================
           Performs operation using MAGMA
           =================================================================== */
        magma_zsetmatrix( m, n, hA, lda, dA, ldda );
        magma_zgeqrf_gpu( m, n, dA, ldda, tau, dT, &info );
        if ( info != 0 )
            printf("magma_zgeqrf_gpu return error %d\n", info );
        magma_zgetmatrix( m, n, dA, ldda, hR, lda );
        
        gpu_time = magma_wtime();
        magma_zungqr( m, n, k, hR, lda, tau, dT, nb, &info );
        gpu_time = magma_wtime() - gpu_time;
        if ( info != 0 )
            printf("magma_zungqr_gpu return error %d\n", info );
        
        gpu_perf = gflops / gpu_time;
        
        /* =====================================================================
           Performs operation using LAPACK
           =================================================================== */
        if ( checkres ) {
            error = lapackf77_zlange("f", &m, &n, hA, &lda, work );
            
            lapackf77_zgeqrf( &m, &n, hA, &lda, tau, hwork, &lwork, &info );
            if ( info != 0 )
                printf("lapackf77_zgeqrf return error %d\n", info );
            
            cpu_time = magma_wtime();
            lapackf77_zungqr( &m, &n, &k, hA, &lda, tau, hwork, &lwork, &info );
            cpu_time = magma_wtime() - cpu_time;
            if ( info != 0 )
                printf("lapackf77_zungqr return error %d\n", info );
            
            cpu_perf = gflops / cpu_time;

            // compute relative error |R|/|A| := |Q_magma - Q_lapack|/|A|
            blasf77_zaxpy( &n2, &c_neg_one, hA, &ione, hR, &ione );
            error = lapackf77_zlange("f", &m, &n, hR, &lda, work) / error;
        }
        
        printf("%5d %5d %5d   %7.1f (%7.2f)   %7.1f (%7.2f)   %8.2e\n",
               m, n, k, cpu_perf, cpu_time, gpu_perf, gpu_time, error );
    }
    
    /* Memory clean up */
    TESTING_HOSTFREE( hA );
    TESTING_HOSTFREE( hwork );
    TESTING_FREE( hR );
    TESTING_FREE( tau );

    TESTING_DEVFREE( dA );
    TESTING_DEVFREE( dT );

    /* Shutdown */
    TESTING_CUDA_FINALIZE();
    return 0;
}
