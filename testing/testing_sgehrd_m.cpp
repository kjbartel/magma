/*
    -- MAGMA (version 1.4.0-beta2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       June 2013

       @generated s Fri Jun 28 19:34:06 2013

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

#define PRECISION_s

/* ////////////////////////////////////////////////////////////////////////////
   -- Testing sgehrd
*/
int main( int argc, char** argv)
{
    TESTING_INIT();
    magma_setdevice( 0 );  // without this, T -> dT copy fails
    
    real_Double_t    gflops, gpu_perf, gpu_time, cpu_perf, cpu_time;
    float *h_A, *h_R, *h_Q, *h_work, *tau, *twork, *T, *dT;
    #if defined(PRECISION_z) || defined(PRECISION_c)
    float      *rwork;
    #endif
    float      eps, result[2];
    magma_int_t N, n2, lda, nb, lwork, ltwork, info;
    magma_int_t ione     = 1;
    magma_int_t ISEED[4] = {0,0,0,1};
    
    eps   = lapackf77_slamch( "E" );
    
    magma_opts opts;
    parse_opts( argc, argv, &opts );
    
    printf("    N   CPU GFlop/s (sec)   GPU GFlop/s (sec)   |A-QHQ'|/N|A|   |I-QQ'|/N\n");
    printf("=========================================================================\n");
    for( int i = 0; i < opts.ntest; ++i ) {
        for( int iter = 0; iter < opts.niter; ++iter ) {
            N = opts.nsize[i];
            lda    = N;
            n2     = lda*N;
            nb     = magma_get_sgehrd_nb(N);
            // magma needs larger workspace than lapack, esp. multi-gpu verison
            lwork  = N*(nb + nb*MagmaMaxGPUs);
            gflops = FLOPS_SGEHRD( N ) / 1e9;
            
            TESTING_MALLOC   ( h_A,    float, n2    );
            TESTING_MALLOC   ( tau,    float, N     );
            TESTING_MALLOC   ( T,      float, nb*N  );
            TESTING_HOSTALLOC( h_R,    float, n2    );
            TESTING_HOSTALLOC( h_work, float, lwork );
            TESTING_DEVALLOC ( dT,     float, nb*N  );
            
            /* Initialize the matrices */
            lapackf77_slarnv( &ione, ISEED, &n2, h_A );
            lapackf77_slacpy( MagmaUpperLowerStr, &N, &N, h_A, &lda, h_R, &lda );
            
            /* ====================================================================
               Performs operation using MAGMA
               =================================================================== */
            gpu_time = magma_wtime();
            magma_sgehrd_m( N, ione, N, h_R, lda, tau, h_work, lwork, T, &info);
            gpu_time = magma_wtime() - gpu_time;
            gpu_perf = gflops / gpu_time;
            if (info != 0)
                printf("magma_sgehrd_m returned error %d: %s.\n",
                       (int) info, magma_strerror( info ));
            
            /* =====================================================================
               Check the factorization
               =================================================================== */
            if ( opts.check ) {
                ltwork = 2*(N*N);
                TESTING_HOSTALLOC( h_Q,   float, lda*N  );
                TESTING_MALLOC(    twork, float, ltwork );
                #if defined(PRECISION_z) || defined(PRECISION_c)
                TESTING_MALLOC(    rwork, float,          N      );
                #endif
                
                lapackf77_slacpy(MagmaUpperLowerStr, &N, &N, h_R, &lda, h_Q, &lda);
                for( int j = 0; j < N-1; ++j )
                    for( int i = j+2; i < N; ++i )
                        h_R[i+j*lda] = MAGMA_S_ZERO;
                
                magma_ssetmatrix( nb, N, T, nb, dT, nb );
                
                magma_sorghr(N, ione, N, h_Q, lda, tau, dT, nb, &info);
                if ( info != 0 ) {
                    printf("magma_sorghr returned error %d: %s.\n",
                           (int) info, magma_strerror( info ));
                    exit(1);
                }
                #if defined(PRECISION_z) || defined(PRECISION_c)
                lapackf77_shst01(&N, &ione, &N,
                                 h_A, &lda, h_R, &lda,
                                 h_Q, &lda, twork, &ltwork, rwork, result);
                #else
                lapackf77_shst01(&N, &ione, &N,
                                 h_A, &lda, h_R, &lda,
                                 h_Q, &lda, twork, &ltwork, result);
                #endif
                
                TESTING_HOSTFREE( h_Q );
                TESTING_FREE( twork );
                #if defined(PRECISION_z) || defined(PRECISION_c)
                TESTING_FREE( rwork );
                #endif
            }
            
            /* =====================================================================
               Performs operation using LAPACK
               =================================================================== */
            if ( opts.lapack ) {
                cpu_time = magma_wtime();
                lapackf77_sgehrd(&N, &ione, &N, h_R, &lda, tau, h_work, &lwork, &info);
                cpu_time = magma_wtime() - cpu_time;
                cpu_perf = gflops / cpu_time;
                if (info != 0)
                    printf("lapackf77_sgehrd returned error %d: %s.\n",
                           (int) info, magma_strerror( info ));
            }
            
            /* =====================================================================
               Print performance and error.
               =================================================================== */
            if ( opts.lapack ) {
                printf("%5d   %7.2f (%7.2f)   %7.2f (%7.2f)",
                       (int) N, cpu_perf, cpu_time, gpu_perf, gpu_time );
            }
            else {
                printf("%5d     ---   (  ---  )   %7.2f (%7.2f)",
                       (int) N, gpu_perf, gpu_time );
            }
            if ( opts.check ) {
                printf("   %8.2e        %8.2e\n",
                       result[0]*eps, result[1]*eps );
            }
            else {
                printf("     ---             ---\n");
            }
            
            TESTING_FREE    ( h_A  );
            TESTING_FREE    ( tau  );
            TESTING_FREE    ( T    );
            TESTING_HOSTFREE( h_work);
            TESTING_HOSTFREE( h_R  );
            TESTING_DEVFREE ( dT   );
        }
        if ( opts.niter > 1 ) {
            printf( "\n" );
        }
    }
    
    TESTING_FINALIZE();
    return 0;
}
