/*
    -- MAGMA (version 1.4.0-beta2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       June 2013

       @generated s Fri Jun 28 19:33:43 2013
       @author Chongxiao Cao
*/

// make sure that asserts are enabled
#ifdef NDEBUG
#undef NDEBUG
#endif
#include <assert.h>

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

#define h_A(i,j) (h_A + (i) + (j)*lda)

/* ////////////////////////////////////////////////////////////////////////////
   -- Testing strsm
*/
int main( int argc, char** argv)
{
    TESTING_INIT();

    real_Double_t   gflops, cublas_perf, cublas_time, cpu_perf, cpu_time;
    float          cublas_error, normA, normx, normr, work[1];
    magma_int_t N, info;
    magma_int_t sizeA;
    magma_int_t lda, ldda;
    magma_int_t ione     = 1;
    magma_int_t ISEED[4] = {0,0,0,1};
    magma_int_t *ipiv;

    float *h_A, *h_b, *h_x, *h_xcublas;
    float *d_A, *d_x;
    float c_neg_one = MAGMA_S_NEG_ONE;
    
    magma_opts opts;
    parse_opts( argc, argv, &opts );
    
    printf("uplo = %c, transA = %c, diag = %c\n", opts.uplo, opts.transA, opts.diag );
    printf("    N  CUBLAS Gflop/s (ms)   CPU Gflop/s (ms)   CUBLAS error\n");
    printf("============================================================\n");
    for( int i = 0; i < opts.ntest; ++i ) {
        for( int iter = 0; iter < opts.niter; ++iter ) {
            N = opts.nsize[i];
            gflops = FLOPS_STRSM(opts.side, N, 1) / 1e9;
            lda    = N;
            ldda   = ((lda+31)/32)*32;
            sizeA  = lda*N;
            
            TESTING_MALLOC( ipiv, magma_int_t, N );
            TESTING_MALLOC( h_A,  float, lda*N );
            TESTING_MALLOC( h_b,  float, N );
            TESTING_MALLOC( h_x,  float, N );
            TESTING_MALLOC( h_xcublas, float, N  );
            
            TESTING_DEVALLOC( d_A, float, ldda*N );
            TESTING_DEVALLOC( d_x, float, N  );
            
            /* Initialize the matrices */
            /* Factor A into LU to get well-conditioned triangular matrix.
             * Copy L to U, since L seems okay when used with non-unit diagonal
             * (i.e., from U), while U fails when used with unit diagonal. */
            lapackf77_slarnv( &ione, ISEED, &sizeA, h_A );
            lapackf77_sgetrf( &N, &N, h_A, &lda, ipiv, &info );
            for( int j = 0; j < N; ++j ) {
                for( int i = 0; i < j; ++i ) {
                    *h_A(i,j) = *h_A(j,i);
                }
            }
            
            lapackf77_slarnv( &ione, ISEED, &N, h_b );
            blasf77_scopy( &N, h_b, &ione, h_x, &ione );
            
            /* =====================================================================
               Performs operation using CUDA-BLAS
               =================================================================== */
            magma_ssetmatrix( N, N, h_A, lda, d_A, ldda );
            magma_ssetvector( N, h_x, 1, d_x, 1 );
            
            cublas_time = magma_sync_wtime( NULL );
            cublasStrsv( opts.uplo, opts.transA, opts.diag,
                         N,
                         d_A, ldda,
                         d_x, 1 );
            cublas_time = magma_sync_wtime( NULL ) - cublas_time;
            cublas_perf = gflops / cublas_time;
            
            magma_sgetvector( N, d_x, 1, h_xcublas, 1 );
            
            /* =====================================================================
               Performs operation using CPU BLAS
               =================================================================== */
            if ( opts.lapack ) {
                cpu_time = magma_wtime();
                blasf77_strsv( &opts.uplo, &opts.transA, &opts.diag,
                               &N,
                               h_A, &lda,
                               h_x, &ione );
                cpu_time = magma_wtime() - cpu_time;
                cpu_perf = gflops / cpu_time;
            }
            
            /* =====================================================================
               Check the result
               =================================================================== */
            // ||b - Ax|| / (||A||*||x||)
            // error for CUBLAS
            normA = lapackf77_slange( "F", &N, &N, h_A, &lda, work );
            
            normx = lapackf77_slange( "F", &N, &ione, h_xcublas, &ione, work );
            blasf77_strmv( &opts.uplo, &opts.transA, &opts.diag,
                           &N,
                           h_A, &lda,
                           h_xcublas, &ione );
            blasf77_saxpy( &N, &c_neg_one, h_b, &ione, h_xcublas, &ione );
            normr = lapackf77_slange( "F", &N, &ione, h_xcublas, &N, work );
            cublas_error = normr / (normA*normx);

            if ( opts.lapack ) {
                printf("%5d   %7.2f (%7.2f)   %7.2f (%7.2f)   %8.2e\n",
                        (int) N,
                        cublas_perf, 1000.*cublas_time,
                        cpu_perf,    1000.*cpu_time,
                        cublas_error );
            }
            else {
                printf("%5d   %7.2f (%7.2f)     ---  (  ---  )   %8.2e\n",
                        (int) N,
                        cublas_perf, 1000.*cublas_time,
                        cublas_error );
            }
            
            TESTING_FREE( h_A );
            TESTING_FREE( h_x );
            TESTING_FREE( h_xcublas );
            
            TESTING_DEVFREE( d_A );
            TESTING_DEVFREE( d_x );
        }
        if ( opts.niter > 1 ) {
            printf( "\n" );
        }
    }

    TESTING_FINALIZE();
    return 0;
}
