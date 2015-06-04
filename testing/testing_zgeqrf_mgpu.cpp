/*
    -- MAGMA (version 1.4.0-beta2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       June 2013

       @precisions normal z -> c d s
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

#define PRECISION_z

/* ////////////////////////////////////////////////////////////////////////////
   -- Testing zgeqrf_mgpu
*/
int main( int argc, char** argv )
{
    TESTING_INIT();

    real_Double_t    gflops, gpu_perf, gpu_time, cpu_perf, cpu_time;
    double           error, work[1];
    magmaDoubleComplex c_neg_one = MAGMA_Z_NEG_ONE;
    magmaDoubleComplex *h_A, *h_R, *tau, *h_work, tmp[1];
    magmaDoubleComplex *d_lA[ MagmaMaxGPUs ];
    magma_int_t M, N, n2, lda, ldda, n_local, ngpu;
    magma_int_t info, min_mn, nb, lhwork;
    magma_int_t ione     = 1;
    magma_int_t ISEED[4] = {0,0,0,1};
    
    magma_opts opts;
    parse_opts( argc, argv, &opts );
    opts.lapack |= opts.check;  // check (-c) implies lapack (-l)
    
    printf("ngpu %d\n", opts.ngpu );
    printf("    M     N   CPU GFlop/s (sec)   GPU GFlop/s (sec)   ||R||_F / ||A||_F\n");
    printf("=======================================================================\n");
    for( int i = 0; i < opts.ntest; ++i ) {
        for( int iter = 0; iter < opts.niter; ++iter ) {
            M = opts.msize[i];
            N = opts.nsize[i];
            min_mn = min(M, N);
            lda    = M;
            n2     = lda*N;
            ldda   = ((M+31)/32)*32;
            nb     = magma_get_zgeqrf_nb( M );
            gflops = FLOPS_ZGEQRF( M, N ) / 1e9;
            
            // ngpu must be at least the number of blocks
            ngpu = min( opts.ngpu, int((N+nb-1)/nb) );
            if ( ngpu < opts.ngpu ) {
                printf( " * too many GPUs for the matrix size, using %d GPUs\n", (int) ngpu );
            }
            
            // query for workspace size
            lhwork = -1;
            lapackf77_zgeqrf( &M, &N, h_A, &M, tau, tmp, &lhwork, &info );
            lhwork = (magma_int_t) MAGMA_Z_REAL( tmp[0] );
            
            // Allocate host memory for the matrix
            TESTING_MALLOC(    tau,    magmaDoubleComplex, min_mn );
            TESTING_MALLOC(    h_A,    magmaDoubleComplex, n2     );
            TESTING_HOSTALLOC( h_R,    magmaDoubleComplex, n2     );
            TESTING_MALLOC(    h_work, magmaDoubleComplex, lhwork );
            
            // Allocate device memory
            for( int dev = 0; dev < ngpu; dev++ ) {
                n_local = ((N/nb)/ngpu)*nb;
                if (dev < (N/nb) % ngpu)
                    n_local += nb;
                else if (dev == (N/nb) % ngpu)
                    n_local += N % nb;
                magma_setdevice( dev );
                TESTING_DEVALLOC(  d_lA[dev], magmaDoubleComplex, ldda*n_local );
            }
            
            /* Initialize the matrix */
            lapackf77_zlarnv( &ione, ISEED, &n2, h_A );
            lapackf77_zlacpy( MagmaUpperLowerStr, &M, &N, h_A, &lda, h_R, &lda );
            
            /* =====================================================================
               Performs operation using LAPACK
               =================================================================== */
            if ( opts.lapack ) {
                cpu_time = magma_wtime();
                lapackf77_zgeqrf( &M, &N, h_A, &M, tau, h_work, &lhwork, &info );
                cpu_time = magma_wtime() - cpu_time;
                cpu_perf = gflops / cpu_time;
                if (info != 0)
                    printf("lapack_zgeqrf returned error %d: %s.\n",
                           (int) info, magma_strerror( info ));
            }
            
            /* ====================================================================
               Performs operation using MAGMA
               =================================================================== */
            magma_zsetmatrix_1D_col_bcyclic( M, N, h_R, lda, d_lA, ldda, ngpu, nb );
            
            gpu_time = magma_wtime();
            magma_zgeqrf2_mgpu( ngpu, M, N, d_lA, ldda, tau, &info );
            gpu_time = magma_wtime() - gpu_time;
            gpu_perf = gflops / gpu_time;
            if (info != 0)
                printf("magma_zgeqrf2 returned error %d: %s.\n",
                       (int) info, magma_strerror( info ));
            
            magma_zgetmatrix_1D_col_bcyclic( M, N, d_lA, ldda, h_R, lda, ngpu, nb );
            magma_queue_sync( NULL );
            
            /* =====================================================================
               Check the result compared to LAPACK
               =================================================================== */
            if ( opts.lapack ) {
                error = lapackf77_zlange("f", &M, &N, h_A, &lda, work );
                blasf77_zaxpy( &n2, &c_neg_one, h_A, &ione, h_R, &ione );
                error = lapackf77_zlange("f", &M, &N, h_R, &lda, work ) / error;
                
                printf("%5d %5d   %7.2f (%7.2f)   %7.2f (%7.2f)   %8.2e\n",
                       (int) M, (int) N, cpu_perf, cpu_time, gpu_perf, gpu_time, error );
            }
            else {
                printf("%5d %5d     ---   (  ---  )   %7.2f (%7.2f)     ---\n",
                       (int) M, (int) N, gpu_perf, gpu_time );
            }
            
            TESTING_FREE( tau );
            TESTING_FREE( h_A );
            TESTING_FREE( h_work );
            TESTING_HOSTFREE( h_R );
            for( int dev=0; dev < ngpu; dev++ ){
                magma_setdevice( dev );
                TESTING_DEVFREE( d_lA[dev] );
            }
        }
        if ( opts.niter > 1 ) {
            printf( "\n" );
        }
    }
    
    TESTING_FINALIZE();
    return 0;
}
