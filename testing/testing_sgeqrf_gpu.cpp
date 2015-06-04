/*
    -- MAGMA (version 1.4.0-beta2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       June 2013

       @generated s Fri Jun 28 19:33:55 2013
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
   -- Testing sgeqrf
*/
int main( int argc, char** argv)
{
    TESTING_INIT();

    real_Double_t    gflops, gpu_perf, gpu_time, cpu_perf, cpu_time;
    float           error, work[1];
    float  c_neg_one = MAGMA_S_NEG_ONE;
    float *h_A, *h_R, *tau, *h_work, tmp[1];
    float *d_A;
    magma_int_t M, N, n2, lda, ldda, lwork, info, min_mn, nb, size;
    magma_int_t ione     = 1;
    magma_int_t ISEED[4] = {0,0,0,1};
    
    magma_opts opts;
    parse_opts( argc, argv, &opts );
    opts.lapack |= opts.check;  // check (-c) implies lapack (-l)
    
    printf("  M     N     CPU GFlop/s (sec)   GPU GFlop/s (sec)   ||R||_F / ||A||_F\n");
    printf("=======================================================================\n");
    for( int i = 0; i < opts.ntest; ++i ) {
        for( int iter = 0; iter < opts.niter; ++iter ) {
            M = opts.msize[i];
            N = opts.nsize[i];
            min_mn = min(M, N);
            lda    = M;
            n2     = lda*N;
            ldda   = ((M+31)/32)*32;
            gflops = FLOPS_SGEQRF( M, N ) / 1e9;
            
            lwork = -1;
            lapackf77_sgeqrf(&M, &N, h_A, &M, tau, tmp, &lwork, &info);
            lwork = (magma_int_t)MAGMA_S_REAL( tmp[0] );
            
            TESTING_MALLOC(    tau, float, min_mn );
            TESTING_MALLOC(    h_A, float, n2     );
            TESTING_HOSTALLOC( h_R, float, n2     );
            TESTING_DEVALLOC(  d_A, float, ldda*N );
            TESTING_MALLOC( h_work, float, lwork  );
            
            /* Initialize the matrix */
            lapackf77_slarnv( &ione, ISEED, &n2, h_A );
            lapackf77_slacpy( MagmaUpperLowerStr, &M, &N, h_A, &lda, h_R, &lda );
            magma_ssetmatrix( M, N, h_R, lda, d_A, ldda );
            
            /* ====================================================================
               Performs operation using MAGMA
               =================================================================== */
            gpu_time = magma_wtime();
            if ( opts.version == 2 ) {
                magma_sgeqrf2_gpu( M, N, d_A, ldda, tau, &info);
            }
            else {
                printf( "NOTE: residual checks for this version will be wrong!\n"
                        "Because tester ignores the special structure of MAGMA sgeqrf resuls.\n"
                        "Use --version 2 to see correct residuals.\n" );
                float *dT;
                nb = magma_get_sgeqrf_nb( M );
                size = (2*min(M, N) + (N+31)/32*32 )*nb;
                magma_smalloc( &dT, size );
                if ( opts.version == 3 ) {
                    magma_sgeqrf3_gpu( M, N, d_A, ldda, tau, dT, &info);
                }
                else {
                    magma_sgeqrf_gpu( M, N, d_A, ldda, tau, dT, &info);
                }
                magma_free( dT );
            }
            gpu_time = magma_wtime() - gpu_time;
            gpu_perf = gflops / gpu_time;
            if (info != 0)
                printf("magma_sgeqrf returned error %d: %s.\n",
                       (int) info, magma_strerror( info ));
            
            if ( opts.lapack ) {
                /* =====================================================================
                   Performs operation using LAPACK
                   =================================================================== */
                cpu_time = magma_wtime();
                lapackf77_sgeqrf(&M, &N, h_A, &lda, tau, h_work, &lwork, &info);
                cpu_time = magma_wtime() - cpu_time;
                cpu_perf = gflops / cpu_time;
                if (info != 0)
                    printf("lapackf77_sgeqrf returned error %d: %s.\n",
                           (int) info, magma_strerror( info ));
                
                /* =====================================================================
                   Check the result compared to LAPACK
                   =================================================================== */
                magma_sgetmatrix( M, N, d_A, ldda, h_R, M );
                error = lapackf77_slange("f", &M, &N, h_A, &lda, work);
                blasf77_saxpy(&n2, &c_neg_one, h_A, &ione, h_R, &ione);
                error = lapackf77_slange("f", &M, &N, h_R, &lda, work) / error;
                
                printf("%5d %5d   %7.2f (%7.2f)   %7.2f (%7.2f)   %8.2e\n",
                       (int) M, (int) N, cpu_perf, cpu_time, gpu_perf, gpu_time, error );
            }
            else {
                printf("%5d %5d     ---   (  ---  )   %7.2f (%7.2f)     ---  \n",
                       (int) M, (int) N, gpu_perf, gpu_time);
            }
            
            TESTING_FREE( tau );
            TESTING_FREE( h_A );
            TESTING_FREE( h_work );
            TESTING_HOSTFREE( h_R );
            TESTING_DEVFREE( d_A );
        }
        if ( opts.niter > 1 ) {
            printf( "\n" );
        }
    }
    
    TESTING_FINALIZE();
    return 0;
}
