/*
    -- MAGMA (version 1.4.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       December 2013

       @generated s Tue Dec 17 13:18:56 2013
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

/* ////////////////////////////////////////////////////////////////////////////
   -- Testing spotf2_gpu
*/
int main( int argc, char** argv)
{
    TESTING_INIT();

    real_Double_t   gflops, gpu_perf, gpu_time, cpu_perf, cpu_time;
    float *h_A, *h_R;
    float *d_A;
    magma_int_t N, n2, lda, ldda, info;
    float c_neg_one = MAGMA_S_NEG_ONE;
    magma_int_t ione     = 1;
    magma_int_t ISEED[4] = {0,0,0,1};
    float      work[1], error;

    magma_opts opts;
    parse_opts( argc, argv, &opts );
    opts.lapack |= opts.check;  // check (-c) implies lapack (-l)
    
    printf("    N   CPU GFlop/s (ms)    GPU GFlop/s (ms)    ||R_magma - R_lapack||_F / ||R_lapack||_F\n");
    printf("========================================================\n");
    for( int i = 0; i < opts.ntest; ++i ) {
        for( int iter = 0; iter < opts.niter; ++iter ) {
            N   = opts.nsize[i];
            lda = N;
            n2  = lda*N;
            ldda = ((N+31)/32)*32;
            gflops = FLOPS_SPOTRF( N ) / 1e9;
            
            if ( N > 512 ) {
                fprintf( stderr, "spotf2 does not support N > 512; skipping N=%d.\n", (int) N );
                continue;
            }
            
            TESTING_MALLOC_CPU( h_A, float, n2     );
            TESTING_MALLOC_PIN( h_R, float, n2     );
            TESTING_MALLOC_DEV( d_A, float, ldda*N );
            
            /* Initialize the matrix */
            lapackf77_slarnv( &ione, ISEED, &n2, h_A );
            magma_smake_hpd( N, h_A, lda );
            lapackf77_slacpy( MagmaUpperLowerStr, &N, &N, h_A, &lda, h_R, &lda );
            magma_ssetmatrix( N, N, h_A, lda, d_A, ldda );
            
            /* ====================================================================
               Performs operation using MAGMA
               =================================================================== */
            gpu_time = magma_wtime();
            magma_spotf2_gpu( opts.uplo, N, d_A, ldda, &info );
            gpu_time = magma_wtime() - gpu_time;
            gpu_perf = gflops / gpu_time;
            if (info != 0)
                printf("magma_spotf2_gpu returned error %d: %s.\n",
                       (int) info, magma_strerror( info ));
            
            if ( opts.lapack ) {
                /* =====================================================================
                   Performs operation using LAPACK
                   =================================================================== */
                cpu_time = magma_wtime();
                lapackf77_spotrf( &opts.uplo, &N, h_A, &lda, &info );
                cpu_time = magma_wtime() - cpu_time;
                cpu_perf = gflops / cpu_time;
                if (info != 0)
                    printf("lapackf77_spotrf returned error %d: %s.\n",
                           (int) info, magma_strerror( info ));
                
                /* =====================================================================
                   Check the result compared to LAPACK
                   =================================================================== */
                magma_sgetmatrix( N, N, d_A, ldda, h_R, lda );
                error = lapackf77_slange("f", &N, &N, h_A, &lda, work);
                blasf77_saxpy(&n2, &c_neg_one, h_A, &ione, h_R, &ione);
                error = lapackf77_slange("f", &N, &N, h_R, &lda, work) / error;
                
                printf("%5d   %7.2f (%7.2f)   %7.2f (%7.2f)   %8.2e\n",
                       (int) N, cpu_perf, cpu_time*1000., gpu_perf, gpu_time*1000., error );
            }
            else {
                printf("%5d     ---   (  ---  )   %7.2f (%7.2f)     ---  \n",
                       (int) N, gpu_perf, gpu_time*1000. );
            }
            TESTING_FREE_CPU( h_A );
            TESTING_FREE_PIN( h_R );
            TESTING_FREE_DEV( d_A );
        }
        if ( opts.niter > 1 ) {
            printf( "\n" );
        }
    }

    TESTING_FINALIZE();
    return 0;
}
