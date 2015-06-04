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

#define PRECISION_s

/* ////////////////////////////////////////////////////////////////////////////
   -- Testing sgels
*/
int main( int argc, char** argv)
{
    TESTING_INIT();
    
    real_Double_t    gflops, gpu_perf, gpu_time, cpu_perf, cpu_time;
    float           gpu_error, cpu_error, matnorm, work[1];
    float  c_one     = MAGMA_S_ONE;
    float  c_neg_one = MAGMA_S_NEG_ONE;
    float *h_A, *h_A2, *h_B, *h_X, *h_R, *tau, *h_work, tmp[1];
    float *d_A, *d_B;
    magma_int_t M, N, n2, nrhs, lda, ldb, ldda, lddb, min_mn, max_mn, nb, info;
    magma_int_t lworkgpu, lhwork, lhwork2;
    magma_int_t ione     = 1;
    magma_int_t ISEED[4] = {0,0,0,1};

    magma_opts opts;
    parse_opts( argc, argv, &opts );
 
    magma_int_t status = 0;
    float tol = opts.tolerance * lapackf77_slamch("E");

    nrhs = opts.nrhs;
    
    printf("                                                            ||b-Ax|| / (N||A||)\n");
    printf("    M     N  NRHS   CPU GFlop/s (sec)   GPU GFlop/s (sec)   CPU        GPU     \n");
    printf("===============================================================================\n");
    for( int i = 0; i < opts.ntest; ++i ) {
        for( int iter = 0; iter < opts.niter; ++iter ) {
            M = opts.msize[i];
            N = opts.nsize[i];
            if ( M < N ) {
                printf( "skipping M=%d, N=%d because M < N is not yet supported.\n", (int) M, (int) N );
                continue;
            }
            min_mn = min(M, N);
            max_mn = max(M, N);
            lda    = M;
            ldb    = max_mn;
            n2     = lda*N;
            ldda   = ((M+31)/32)*32;
            lddb   = ((max_mn+31)/32)*32;
            nb     = magma_get_sgeqrf_nb(M);
            gflops = (FLOPS_SGEQRF( M, N ) + FLOPS_SGEQRS( M, N, nrhs )) / 1e9;
            
            lworkgpu = (M - N + nb)*(nrhs + nb) + nrhs*nb;
            
            // query for workspace size
            lhwork = -1;
            lapackf77_sgeqrf(&M, &N, NULL, &M, NULL, tmp, &lhwork, &info);
            lhwork2 = (magma_int_t) MAGMA_S_REAL( tmp[0] );
            
            lhwork = -1;
            lapackf77_sormqr( MagmaLeftStr, MagmaTransStr,
                              &M, &nrhs, &min_mn, NULL, &lda, NULL,
                              NULL, &ldb, tmp, &lhwork, &info);
            lhwork = (magma_int_t) MAGMA_S_REAL( tmp[0] );
            lhwork = max( max( lhwork, lhwork2 ), lworkgpu );
            
            TESTING_MALLOC_CPU( tau,    float, min_mn    );
            TESTING_MALLOC_CPU( h_A,    float, lda*N     );
            TESTING_MALLOC_CPU( h_A2,   float, lda*N     );
            TESTING_MALLOC_CPU( h_B,    float, ldb*nrhs  );
            TESTING_MALLOC_CPU( h_X,    float, ldb*nrhs  );
            TESTING_MALLOC_CPU( h_R,    float, ldb*nrhs  );
            TESTING_MALLOC_CPU( h_work, float, lhwork    );
            
            TESTING_MALLOC_DEV( d_A,    float, ldda*N    );
            TESTING_MALLOC_DEV( d_B,    float, lddb*nrhs );
            
            /* Initialize the matrices */
            lapackf77_slarnv( &ione, ISEED, &n2, h_A );
            lapackf77_slacpy( MagmaUpperLowerStr, &M, &N, h_A, &lda, h_A2, &lda );
            
            // make random RHS
            n2 = M*nrhs;
            lapackf77_slarnv( &ione, ISEED, &n2, h_B );
            lapackf77_slacpy( MagmaUpperLowerStr, &M, &nrhs, h_B, &ldb, h_R, &ldb );
            
            // make consistent RHS
            //n2 = N*nrhs;
            //lapackf77_slarnv( &ione, ISEED, &n2, h_X );
            //blasf77_sgemm( MagmaNoTransStr, MagmaNoTransStr, &M, &nrhs, &N,
            //               &c_one,  h_A, &lda,
            //                        h_X, &ldb,
            //               &c_zero, h_B, &ldb );
            //lapackf77_slacpy( MagmaUpperLowerStr, &M, &nrhs, h_B, &ldb, h_R, &ldb );
            
            /* ====================================================================
               Performs operation using MAGMA
               =================================================================== */
            magma_ssetmatrix( M, N,    h_A, lda, d_A, ldda );
            magma_ssetmatrix( M, nrhs, h_B, ldb, d_B, lddb );
            
            gpu_time = magma_wtime();
            magma_sgels3_gpu( MagmaNoTrans, M, N, nrhs, d_A, ldda,
                              d_B, lddb, h_work, lworkgpu, &info);
            gpu_time = magma_wtime() - gpu_time;
            gpu_perf = gflops / gpu_time;
            if (info != 0)
                printf("magma_sgels3_gpu returned error %d: %s.\n",
                       (int) info, magma_strerror( info ));
            
            // Get the solution in h_X
            magma_sgetmatrix( N, nrhs, d_B, lddb, h_X, ldb );
            
            // compute the residual
            blasf77_sgemm( MagmaNoTransStr, MagmaNoTransStr, &M, &nrhs, &N,
                           &c_neg_one, h_A, &lda,
                                       h_X, &ldb,
                           &c_one,     h_R, &ldb);
            matnorm = lapackf77_slange("f", &M, &N, h_A, &lda, work);
            
            /* =====================================================================
               Performs operation using LAPACK
               =================================================================== */
            lapackf77_slacpy( MagmaUpperLowerStr, &M, &nrhs, h_B, &ldb, h_X, &ldb );
            
            cpu_time = magma_wtime();
            lapackf77_sgels( MagmaNoTransStr, &M, &N, &nrhs,
                             h_A, &lda, h_X, &ldb, h_work, &lhwork, &info);
            cpu_time = magma_wtime() - cpu_time;
            cpu_perf = gflops / cpu_time;
            if (info != 0)
                printf("lapackf77_sgels returned error %d: %s.\n",
                       (int) info, magma_strerror( info ));
            
            blasf77_sgemm( MagmaNoTransStr, MagmaNoTransStr, &M, &nrhs, &N,
                           &c_neg_one, h_A2, &lda,
                                       h_X,  &ldb,
                           &c_one,     h_B,  &ldb);
            
            cpu_error = lapackf77_slange("f", &M, &nrhs, h_B, &ldb, work) / (min_mn*matnorm);
            gpu_error = lapackf77_slange("f", &M, &nrhs, h_R, &ldb, work) / (min_mn*matnorm);
            
            printf("%5d %5d %5d   %7.2f (%7.2f)   %7.2f (%7.2f)   %8.2e   %8.2e",
                   (int) M, (int) N, (int) nrhs,
                   cpu_perf, cpu_time, gpu_perf, gpu_time, cpu_error, gpu_error );
            printf("%s\n", (gpu_error < tol ? "" : "  failed"));
            status |= ! (gpu_error < tol);
            
            TESTING_FREE_CPU( tau    );
            TESTING_FREE_CPU( h_A    );
            TESTING_FREE_CPU( h_A2   );
            TESTING_FREE_CPU( h_B    );
            TESTING_FREE_CPU( h_X    );
            TESTING_FREE_CPU( h_R    );
            TESTING_FREE_CPU( h_work );
            
            TESTING_FREE_DEV( d_A    );
            TESTING_FREE_DEV( d_B    );
        }
        if ( opts.niter > 1 ) {
            printf( "\n" );
        }
    }

    TESTING_FINALIZE();
    return status;
}
