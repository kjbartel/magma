/*
    -- MAGMA (version 1.6.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date November 2014

       @author Mark Gates
       @generated from testing_zunmqr_gpu.cpp normal z -> s, Sat Nov 15 19:54:18 2014
*/
// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>

// includes, project
#include "flops.h"
#include "magma.h"
#include "magma_lapack.h"
#include "testings.h"

/* ////////////////////////////////////////////////////////////////////////////
   -- Testing sormqr_gpu
*/
int main( int argc, char** argv )
{
    TESTING_INIT();
    
    real_Double_t   gflops, gpu_perf, gpu_time, cpu_perf, cpu_time;
    float error, work[1];
    float c_neg_one = MAGMA_S_NEG_ONE;
    magma_int_t ione = 1;
    magma_int_t m, n, k, size, info;
    magma_int_t ISEED[4] = {0,0,0,1};
    magma_int_t nb, ldc, lda, lwork, lwork_max, dt_size;
    float *C, *R, *A, *W, *tau;
    magmaFloat_ptr dC, dA, dT;
    magma_int_t status = 0;
    
    magma_opts opts;
    parse_opts( argc, argv, &opts );
    
    float tol = 2. * opts.tolerance * lapackf77_slamch("E");
    
    // test all combinations of input parameters
    magma_side_t  side [] = { MagmaLeft,       MagmaRight   };
    magma_trans_t trans[] = { MagmaTrans, MagmaNoTrans };

    printf("    M     N     K   side   trans   CPU GFlop/s (sec)   GPU GFlop/s (sec)   ||R||_F / ||QC||_F\n");
    printf("===============================================================================================\n");
    for( int itest = 0; itest < opts.ntest; ++itest ) {
      for( int iside = 0; iside < 2; ++iside ) {
      for( int itran = 0; itran < 2; ++itran ) {
        for( int iter = 0; iter < opts.niter; ++iter ) {        
            m = opts.msize[itest];
            n = opts.nsize[itest];
            k = opts.ksize[itest];
            nb  = magma_get_sgeqrf_nb( m );
            ldc = ((m + 31)/32)*32;
            lda = ((max(m,n) + 31)/32)*32;
            gflops = FLOPS_SORMQR( m, n, k, side[iside] ) / 1e9;
            
            if ( side[iside] == MagmaLeft && m < k ) {
                printf( "%5d %5d %5d   %4c   %5c   skipping because side=left  and m < k\n",
                        (int) m, (int) n, (int) k,
                        lapacke_side_const( side[iside] ),
                        lapacke_trans_const( trans[itran] ) );
                continue;
            }
            if ( side[iside] == MagmaRight && n < k ) {
                printf( "%5d %5d %5d   %4c   %5c   skipping because side=right and n < k\n",
                        (int) m, (int) n, (int) k,
                        lapacke_side_const( side[iside] ),
                        lapacke_trans_const( trans[itran] ) );
                continue;
            }
            
            if ( side[iside] == MagmaLeft ) {
                // side = left
                lwork_max = (m - k + nb)*(n + nb) + n*nb;
                dt_size = ( 2*min(m,k) + ((max(m,n) + 31)/32)*32 )*nb;
            }
            else {
                // side = right
                lwork_max = (n - k + nb)*(m + nb) + m*nb;
                dt_size = ( 2*min(n,k) + ((max(m,n) + 31)/32)*32 )*nb;
            }
            
            TESTING_MALLOC_CPU( C,   float, ldc*n );
            TESTING_MALLOC_CPU( R,   float, ldc*n );
            TESTING_MALLOC_CPU( A,   float, lda*k );
            TESTING_MALLOC_CPU( W,   float, lwork_max );
            TESTING_MALLOC_CPU( tau, float, k );
            
            TESTING_MALLOC_DEV( dC, float, ldc*n );
            TESTING_MALLOC_DEV( dA, float, lda*k );
            TESTING_MALLOC_DEV( dT, float, dt_size );
            
            // C is full, m x n
            size = ldc*n;
            lapackf77_slarnv( &ione, ISEED, &size, C );
            magma_ssetmatrix( m, n, C, ldc, dC, ldc );
            
            // A is m x k (left) or n x k (right)
            lda = (side[iside] == MagmaLeft ? m : n);
            size = lda*k;
            lapackf77_slarnv( &ione, ISEED, &size, A );
            
            // compute QR factorization to get Householder vectors in dA, tau, dT
            magma_ssetmatrix( lda, k, A,  lda, dA, lda );
            magma_sgeqrf_gpu( lda, k, dA, lda, tau, dT, &info );
            magma_sgetmatrix( lda, k, dA, lda, A,  lda );
            if (info != 0)
                printf("magma_sgeqrf_gpu returned error %d: %s.\n",
                       (int) info, magma_strerror( info ));
            
            /* =====================================================================
               Performs operation using LAPACK
               =================================================================== */
            cpu_time = magma_wtime();
            lapackf77_sormqr( lapack_side_const( side[iside] ), lapack_trans_const( trans[itran] ),
                              &m, &n, &k,
                              A, &lda, tau, C, &ldc, W, &lwork_max, &info );
            cpu_time = magma_wtime() - cpu_time;
            cpu_perf = gflops / cpu_time;
            if (info != 0)
                printf("lapackf77_sormqr returned error %d: %s.\n",
                       (int) info, magma_strerror( info ));
            
            /* ====================================================================
               Performs operation using MAGMA
               =================================================================== */
            // query for workspace size
            lwork = -1;
            magma_sormqr_gpu( side[iside], trans[itran],
                              m, n, k,
                              dA, lda, tau, dC, ldc, W, lwork, dT, nb, &info );
            if (info != 0)
                printf("magma_sormqr_gpu (lwork query) returned error %d: %s.\n",
                       (int) info, magma_strerror( info ));
            lwork = (magma_int_t) MAGMA_S_REAL( W[0] );
            if ( lwork < 0 || lwork > lwork_max )
                printf("invalid lwork %d, lwork_max %d\n", (int) lwork, (int) lwork_max );
            
            gpu_time = magma_sync_wtime( 0 );  // sync needed for L,N and R,T cases
            magma_sormqr_gpu( side[iside], trans[itran],
                              m, n, k,
                              dA, lda, tau, dC, ldc, W, lwork, dT, nb, &info );
            gpu_time = magma_sync_wtime( 0 ) - gpu_time;
            gpu_perf = gflops / gpu_time;
            if (info != 0)
                printf("magma_sormqr_gpu returned error %d: %s.\n",
                       (int) info, magma_strerror( info ));
            
            magma_sgetmatrix( m, n, dC, ldc, R, ldc );
            
            /* =====================================================================
               compute relative error |QC_magma - QC_lapack| / |QC_lapack|
               =================================================================== */
            error = lapackf77_slange( "Fro", &m, &n, C, &ldc, work );
            size = ldc*n;
            blasf77_saxpy( &size, &c_neg_one, C, &ione, R, &ione );
            error = lapackf77_slange( "Fro", &m, &n, R, &ldc, work ) / error;
            
            printf( "%5d %5d %5d   %4c   %5c   %7.2f (%7.2f)   %7.2f (%7.2f)   %8.2e   %s\n",
                    (int) m, (int) n, (int) k,
                    lapacke_side_const( side[iside] ),
                    lapacke_trans_const( trans[itran] ),
                    cpu_perf, cpu_time, gpu_perf, gpu_time,
                    error, (error < tol ? "ok" : "failed") );
            status += ! (error < tol);
            
            TESTING_FREE_CPU( C );
            TESTING_FREE_CPU( R );
            TESTING_FREE_CPU( A );
            TESTING_FREE_CPU( W );
            TESTING_FREE_CPU( tau );
            
            TESTING_FREE_DEV( dC );
            TESTING_FREE_DEV( dA );
            TESTING_FREE_DEV( dT );
            fflush( stdout );
        }
        if ( opts.niter > 1 ) {
            printf( "\n" );
        }
      }}  // end iside, itran
      printf( "\n" );
    }
    
    TESTING_FINALIZE();
    return status;
}
