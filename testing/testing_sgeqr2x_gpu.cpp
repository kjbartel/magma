/*
    -- MAGMA (version 1.6.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date November 2014

       @generated from testing_zgeqr2x_gpu.cpp normal z -> s, Sat Nov 15 19:54:18 2014

*/

// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

// includes, project
#include "flops.h"
#include "magma.h"
#include "magma_lapack.h"
#include "testings.h"
#include "common_magma.h"

/* ////////////////////////////////////////////////////////////////////////////
   -- Testing sgeqrf
*/
int main( int argc, char** argv)
{
    TESTING_INIT();

    real_Double_t    gflops, gpu_perf, gpu_time, cpu_perf, cpu_time;
    float           error, work[1];

    float  c_neg_one = MAGMA_S_NEG_ONE;
    float *h_A, *h_T, *h_R, *tau, *h_work, tmp[1];
    magmaFloat_ptr d_A,  d_T, ddA, dtau;
    magmaFloat_ptr d_A2, d_T2, ddA2, dtau2;
    magmaFloat_ptr dwork, dwork2;

    magma_int_t M, N, lda, ldda, lwork, n2, info, min_mn;
    magma_int_t ione     = 1;
    magma_int_t ISEED[4] = {0,0,0,1};
    magma_int_t status = 0;

    #define BLOCK_SIZE 64

    magma_opts opts;
    parse_opts( argc, argv, &opts );
    
    float tol = 10. * opts.tolerance * lapackf77_slamch("E");
    
    magma_queue_t stream[2];
    magma_queue_create( &stream[0] );
    magma_queue_create( &stream[1] );

    printf("version %d\n", (int) opts.version );
    printf("  M     N     CPU GFlop/s (ms)    GPU GFlop/s (ms)   ||R||_F/||A||_F  ||R_T||\n");
    printf("=============================================================================\n");
    for( int itest = 0; itest < opts.ntest; ++itest ) {
        for( int iter = 0; iter < opts.niter; ++iter ) {
            M     = opts.msize[itest];
            N     = opts.nsize[itest];

            if (N > 128) {
                printf("%5d %5d   skipping because sgeqr2x requires N <= 128\n",
                        (int) M, (int) N);
                continue;
            }
            if (M < N) {
                printf("%5d %5d   skipping because sgeqr2x requires M >= N\n",
                        (int) M, (int) N);
                continue;
            }

            min_mn = min(M, N);
            lda    = M;
            n2     = lda*N;
            ldda   = ((M+31)/32)*32;
            gflops = (FLOPS_SGEQRF( M, N ) + FLOPS_SGEQRT( M, N )) / 1e9;

            /* Allocate memory for the matrix */
            TESTING_MALLOC_CPU( tau,   float, min_mn );
            TESTING_MALLOC_CPU( h_A,   float, n2     );
            TESTING_MALLOC_CPU( h_T,   float, N*N    );
        
            TESTING_MALLOC_PIN( h_R,   float, n2     );
        
            TESTING_MALLOC_DEV( d_A,   float, ldda*N );
            TESTING_MALLOC_DEV( d_T,   float, N*N    );
            TESTING_MALLOC_DEV( ddA,   float, N*N    );
            TESTING_MALLOC_DEV( dtau,  float, min_mn );
        
            TESTING_MALLOC_DEV( d_A2,  float, ldda*N );
            TESTING_MALLOC_DEV( d_T2,  float, N*N    );
            TESTING_MALLOC_DEV( ddA2,  float, N*N    );
            TESTING_MALLOC_DEV( dtau2, float, min_mn );
        
            TESTING_MALLOC_DEV( dwork,  float, max(5*min_mn, (BLOCK_SIZE*2+2)*min_mn) );
            TESTING_MALLOC_DEV( dwork2, float, max(5*min_mn, (BLOCK_SIZE*2+2)*min_mn) );
            
            // todo replace with magma_slaset
            cudaMemset(ddA, 0, N*N*sizeof(float));
            cudaMemset(d_T, 0, N*N*sizeof(float));
        
            cudaMemset(ddA2, 0, N*N*sizeof(float));
            cudaMemset(d_T2, 0, N*N*sizeof(float));
        
            lwork = -1;
            lapackf77_sgeqrf(&M, &N, NULL, &M, NULL, tmp, &lwork, &info);
            lwork = (magma_int_t)MAGMA_S_REAL( tmp[0] );
            lwork = max(lwork, N*N);
        
            TESTING_MALLOC_CPU( h_work, float, lwork );

            /* Initialize the matrix */
            lapackf77_slarnv( &ione, ISEED, &n2, h_A );
            lapackf77_slacpy( MagmaUpperLowerStr, &M, &N, h_A, &lda, h_R, &lda );
            magma_ssetmatrix( M, N, h_R, lda,  d_A, ldda );
            magma_ssetmatrix( M, N, h_R, lda, d_A2, ldda );
    
            /* ====================================================================
               Performs operation using MAGMA
               =================================================================== */
            gpu_time = magma_sync_wtime(0);
    
            if (opts.version == 1)
                magma_sgeqr2x_gpu(M, N, d_A, ldda, dtau, d_T, ddA, dwork, &info);
            else if (opts.version == 2)
                magma_sgeqr2x2_gpu(M, N, d_A, ldda, dtau, d_T, ddA, dwork, &info);
            else if (opts.version == 3)
                magma_sgeqr2x3_gpu(M, N, d_A, ldda, dtau, d_T, ddA, dwork, &info);
            else {
                printf( "call magma_sgeqr2x4_gpu\n" );
                /*
                  Going through NULL stream is faster
                  Going through any stream is slower
                  Doing two streams in parallel is slower than doing them sequentially
                  Queuing happens on the NULL stream - user defined buffers are smaller?
                */
                magma_sgeqr2x4_gpu(M, N, d_A, ldda, dtau, d_T, ddA, dwork, NULL, &info);
                //magma_sgeqr2x4_gpu(M, N, d_A, ldda, dtau, d_T, ddA, dwork, &info, stream[1]);
                //magma_sgeqr2x4_gpu(M, N, d_A2, ldda, dtau2, d_T2, ddA2, dwork2, &info, stream[0]);
                //magma_sgeqr2x4_gpu(M, N, d_A2, ldda, dtau2, d_T2, ddA2, dwork2, &info, NULL);
                //gflops *= 2;
            }
            gpu_time = magma_sync_wtime(0) - gpu_time;
            gpu_perf = gflops / gpu_time;

            if (info != 0) {
                printf("magma_sgeqr2x_gpu version %d returned error %d: %s.\n",
                       (int) opts.version, (int) info, magma_strerror( info ));
            } 
            else {
                if ( opts.check ) {
                    /* =====================================================================
                       Performs operation using LAPACK
                       =================================================================== */
                    cpu_time = magma_wtime();
                    lapackf77_sgeqrf(&M, &N, h_A, &lda, tau, h_work, &lwork, &info);
                    lapackf77_slarft( MagmaForwardStr, MagmaColumnwiseStr,
                                     &M, &N, h_A, &lda, tau, h_work, &N);
                    //magma_sgeqr2(&M, &N, h_A, &lda, tau, h_work, &info);
                    cpu_time = magma_wtime() - cpu_time;
                    cpu_perf = gflops / cpu_time;
                    if (info != 0)
                        printf("lapackf77_sgeqrf returned error %d: %s.\n",
                               (int) info, magma_strerror( info ));
                
                    /* =====================================================================
                       Check the result compared to LAPACK
                       =================================================================== */
                    magma_sgetmatrix( M, N, d_A, ldda, h_R, M );
                    magma_sgetmatrix( N, N, ddA, N,    h_T, N );
    
                    // Restore the upper triangular part of A before the check
                    for(int col=0; col < N; col++){
                        for(int row=0; row <= col; row++)
                            h_R[row + col*M] = h_T[row + col*N];
                    }
                
                    error = lapackf77_slange("M", &M, &N, h_A, &lda, work);
                    blasf77_saxpy(&n2, &c_neg_one, h_A, &ione, h_R, &ione);
                    error = lapackf77_slange("M", &M, &N, h_R, &lda, work) / (N * error);
     
                    // Check if T is the same
                    magma_sgetmatrix( N, N, d_T, N, h_T, N );
    
                    float terr = 0.;
                    for(int col=0; col < N; col++)
                        for(int row=0; row <= col; row++)
                            terr += (  MAGMA_S_ABS(h_work[row + col*N] - h_T[row + col*N])*
                                       MAGMA_S_ABS(h_work[row + col*N] - h_T[row + col*N])  );
                    terr = magma_ssqrt(terr);
    
                    printf("%5d %5d   %7.2f (%7.2f)   %7.2f (%7.2f)     %8.2e     %8.2e   %s\n",
                           (int) M, (int) N, cpu_perf, 1000.*cpu_time, gpu_perf, 1000.*gpu_time,
                           error, terr, (error < tol ? "ok" : "failed") );
                    status += ! (error < tol);
                }
                else {
                    printf("%5d %5d     ---   (  ---  )   %7.2f (%7.2f)     ---  \n",
                           (int) M, (int) N, gpu_perf, 1000.*gpu_time);
                }
            }
            
            TESTING_FREE_CPU( tau    );
            TESTING_FREE_CPU( h_A    );
            TESTING_FREE_CPU( h_T    );
            TESTING_FREE_CPU( h_work );
            
            TESTING_FREE_PIN( h_R    );
        
            TESTING_FREE_DEV( d_A   );
            TESTING_FREE_DEV( d_T   );
            TESTING_FREE_DEV( ddA   );
            TESTING_FREE_DEV( dtau  );
            TESTING_FREE_DEV( dwork );
        
            TESTING_FREE_DEV( d_A2   );
            TESTING_FREE_DEV( d_T2   );
            TESTING_FREE_DEV( ddA2   );
            TESTING_FREE_DEV( dtau2  );
            TESTING_FREE_DEV( dwork2 );
            fflush( stdout );
        }
        if ( opts.niter > 1 ) {
            printf( "\n" );
        }
    }
    
    magma_queue_destroy( stream[0] );
    magma_queue_destroy( stream[1] );

    TESTING_FINALIZE();
    return status;
}
