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

double get_LU_error(magma_int_t M, magma_int_t N, 
                    magmaDoubleComplex *A,  magma_int_t lda, 
                    magmaDoubleComplex *LU, magma_int_t *IPIV)
{
    magma_int_t min_mn = min(M,N);
    magma_int_t ione   = 1;
    magma_int_t i, j;
    magmaDoubleComplex alpha = MAGMA_Z_ONE;
    magmaDoubleComplex beta  = MAGMA_Z_ZERO;
    magmaDoubleComplex *L, *U;
    double work[1], matnorm, residual;
                       
    TESTING_MALLOC( L, magmaDoubleComplex, M*min_mn);
    TESTING_MALLOC( U, magmaDoubleComplex, min_mn*N);
    memset( L, 0, M*min_mn*sizeof(magmaDoubleComplex) );
    memset( U, 0, min_mn*N*sizeof(magmaDoubleComplex) );

    lapackf77_zlaswp( &N, A, &lda, &ione, &min_mn, IPIV, &ione);
    lapackf77_zlacpy( MagmaLowerStr, &M, &min_mn, LU, &lda, L, &M      );
    lapackf77_zlacpy( MagmaUpperStr, &min_mn, &N, LU, &lda, U, &min_mn );

    for(j=0; j<min_mn; j++)
        L[j+j*M] = MAGMA_Z_MAKE( 1., 0. );
    
    matnorm = lapackf77_zlange("f", &M, &N, A, &lda, work);

    blasf77_zgemm("N", "N", &M, &N, &min_mn,
                  &alpha, L, &M, U, &min_mn, &beta, LU, &lda);

    for( j = 0; j < N; j++ ) {
        for( i = 0; i < M; i++ ) {
            LU[i+j*lda] = MAGMA_Z_SUB( LU[i+j*lda], A[i+j*lda] );
        }
    }
    residual = lapackf77_zlange("f", &M, &N, LU, &lda, work);

    TESTING_FREE(L);
    TESTING_FREE(U);

    return residual / (matnorm * N);
}

/* ////////////////////////////////////////////////////////////////////////////
   -- Testing zgetrf_mgpu
*/
int main( int argc, char** argv )
{
    TESTING_INIT();

    real_Double_t    gflops, gpu_perf, gpu_time, cpu_perf, cpu_time;
    double           error;
    magmaDoubleComplex *h_A, *h_R;
    magmaDoubleComplex *d_lA[ MagmaMaxGPUs ];
    magma_int_t *ipiv;
    magma_int_t M, N, n2, lda, ldda, n_local, ngpu;
    magma_int_t info, min_mn, nb, ldn_local;
    magma_int_t ione     = 1;
    magma_int_t ISEED[4] = {0,0,0,1};

    magma_opts opts;
    parse_opts( argc, argv, &opts );

    printf("ngpu %d\n", opts.ngpu );
    printf("    M     N   CPU GFlop/s (sec)   GPU GFlop/s (sec)   ||PA-LU||/(||A||*N)\n");
    printf("=========================================================================\n");
    for( int i = 0; i < opts.ntest; ++i ) {
        for( int iter = 0; iter < opts.niter; ++iter ) {
            M = opts.msize[i];
            N = opts.nsize[i];
            min_mn = min(M, N);
            lda    = M;
            n2     = lda*N;
            ldda   = ((M+31)/32)*32;
            nb     = magma_get_zgetrf_nb( M );
            gflops = FLOPS_ZGETRF( M, N ) / 1e9;
            
            // ngpu must be at least the number of blocks
            ngpu = min( opts.ngpu, int((N+nb-1)/nb) );
            if ( ngpu < opts.ngpu ) {
                printf( " * too many GPUs for the matrix size, using %d GPUs\n", (int) ngpu );
            }
            
            // Allocate host memory for the matrix
            TESTING_MALLOC(    ipiv, magma_int_t,        min_mn );
            TESTING_MALLOC(    h_A,  magmaDoubleComplex, n2     );
            TESTING_HOSTALLOC( h_R,  magmaDoubleComplex, n2     );
            
            // Allocate device memory
            for( int dev=0; dev < ngpu; dev++){
                n_local = ((N/nb)/ngpu)*nb;
                if (dev < (N/nb) % ngpu)
                    n_local += nb;
                else if (dev == (N/nb) % ngpu)
                    n_local += N % nb;
                ldn_local = ((n_local+31)/32)*32;  // TODO why?
                magma_setdevice( dev );
                TESTING_DEVALLOC( d_lA[dev], magmaDoubleComplex, ldda*ldn_local );
            }
        
            /* Initialize the matrix */
            lapackf77_zlarnv( &ione, ISEED, &n2, h_A );
            lapackf77_zlacpy( MagmaUpperLowerStr, &M, &N, h_A, &lda, h_R, &lda );
    
            /* =====================================================================
               Performs operation using LAPACK
               =================================================================== */
            if ( opts.lapack ) {
                cpu_time = magma_wtime();
                lapackf77_zgetrf( &M, &N, h_A, &lda, ipiv, &info );
                cpu_time = magma_wtime() - cpu_time;
                cpu_perf = gflops / cpu_time;
                if (info != 0)
                    printf("lapackf77_zgetrf returned error %d: %s.\n",
                           (int) info, magma_strerror( info ));
                // restore A for error check
                lapackf77_zlacpy( MagmaUpperLowerStr, &M, &N, h_R, &lda, h_A, &lda );
            }
            
            /* ====================================================================
               Performs operation using MAGMA
               =================================================================== */
            magma_zsetmatrix_1D_col_bcyclic( M, N, h_R, lda, d_lA, ldda, ngpu, nb );
    
            gpu_time = magma_wtime();
            magma_zgetrf_mgpu( ngpu, M, N, d_lA, ldda, ipiv, &info );
            gpu_time = magma_wtime() - gpu_time;
            gpu_perf = gflops / gpu_time;
            if (info != 0)
                printf("magma_zgetrf_mgpu returned error %d: %s.\n",
                       (int) info, magma_strerror( info ));
                       
            magma_zgetmatrix_1D_col_bcyclic( M, N, d_lA, ldda, h_R, lda, ngpu, nb );
    
            /* =====================================================================
               Check the factorization
               =================================================================== */
            if ( opts.lapack ) {
                printf("%5d %5d  %7.2f (%7.2f)   %7.2f (%7.2f)",
                       (int) M, (int) N, cpu_perf, cpu_time, gpu_perf, gpu_time );
            }
            else {
                printf("%5d %5d    ---   (  ---  )   %7.2f (%7.2f)",
                       (int) M, (int) N, gpu_perf, gpu_time );
            }
            if ( opts.check ) {
                error = get_LU_error( M, N, h_A, lda, h_R, ipiv );
                printf( "   %8.2e\n", error );
            }
            else {
                printf( "     ---\n" );
            }
            
            TESTING_FREE( ipiv );
            TESTING_FREE( h_A );
            TESTING_HOSTFREE( h_R );
            for( int dev=0; dev < ngpu; dev++ ) {
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
