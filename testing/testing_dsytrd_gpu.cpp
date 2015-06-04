/*
    -- MAGMA (version 1.4.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       August 2013

       @author Raffaele Solca
       @author Stan Tomov
       @author Azzam Haidar

       @generated d Wed Aug 14 12:18:08 2013

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

#define PRECISION_d

/* ////////////////////////////////////////////////////////////////////////////
   -- Testing dsytrd_gpu
*/
int main( int argc, char** argv)
{
    TESTING_INIT();

    real_Double_t    gflops, gpu_perf, cpu_perf, gpu_time, cpu_time;
    double           eps;
    double *h_A, *h_R, *d_R, *h_Q, *h_work, *work, *dwork;
    double *tau;
    double          *diag, *offdiag;
    double           result[2] = {0., 0.};
    magma_int_t N, n2, lda, lwork, info, nb, ldwork;
    magma_int_t ione     = 1;
    magma_int_t itwo     = 2;
    magma_int_t ithree   = 3;
    magma_int_t ISEED[4] = {0,0,0,1};
    magma_int_t status = 0;
    
    #if defined(PRECISION_z) || defined(PRECISION_c)
    double *rwork;
    #endif

    eps = lapackf77_dlamch( "E" );

    magma_opts opts;
    parse_opts( argc, argv, &opts );

    double tol = opts.tolerance * lapackf77_dlamch("E");
    
    printf("Running version %d; available are (specified through --version num):\n", 
           (int) opts.version);
    printf("1 - uses DSYMV from CUBLAS (default)\n");
    printf("2 - uses DSYMV from MAGMA BLAS that requires extra space\n\n");

    printf("  N     CPU GFlop/s (sec)   GPU GFlop/s (sec)   |A-QHQ'|/N|A|   |I-QQ'|/N\n");
    printf("===========================================================================\n");
    for( int i = 0; i < opts.ntest; ++i ) {
        for( int iter = 0; iter < opts.niter; ++iter ) {
            N = opts.nsize[i];
            lda    = N;
            n2     = N*lda;
            nb     = magma_get_dsytrd_nb(N);
            lwork  = N*nb;  /* We suppose the magma nb is bigger than lapack nb */
            gflops = FLOPS_DSYTRD( N ) / 1e9;
            ldwork = (N*N+64-1)/64 + 2*N*nb;
            
            TESTING_MALLOC(    h_A,     double, lda*N );
            TESTING_DEVALLOC(  d_R,     double, lda*N );
            TESTING_DEVALLOC(dwork,     double, ldwork);
            TESTING_HOSTALLOC( h_R,     double, lda*N );
            TESTING_HOSTALLOC( h_work,  double, lwork );
            TESTING_MALLOC(    tau,     double, N     );
            TESTING_MALLOC(    diag,    double, N   );
            TESTING_MALLOC(    offdiag, double, N-1 );
            
            if ( opts.check ) {
                TESTING_MALLOC( h_Q,  double, lda*N );
                TESTING_MALLOC( work, double, 2*N*N );
                #if defined(PRECISION_z) || defined(PRECISION_c)
                TESTING_MALLOC( rwork, double, N );
                #endif
            }
            
            /* ====================================================================
               Initialize the matrix
               =================================================================== */
            lapackf77_dlarnv( &ione, ISEED, &n2, h_A );
            magma_dmake_symmetric( N, h_A, lda );
            magma_dsetmatrix( N, N, h_A, lda, d_R, lda );
            
            /* ====================================================================
               Performs operation using MAGMA
               =================================================================== */
            gpu_time = magma_wtime();
            if (opts.version == 1)
                magma_dsytrd_gpu( opts.uplo, N, d_R, lda, diag, offdiag,
                                  tau, h_R, lda, h_work, lwork, &info );
            else
                magma_dsytrd2_gpu( opts.uplo, N, d_R, lda, diag, offdiag,
                                   tau, h_R, lda, h_work, lwork, dwork, ldwork, &info );
            gpu_time = magma_wtime() - gpu_time;
            gpu_perf = gflops / gpu_time;
            if (info != 0)
                printf("magma_dsytrd_gpu returned error %d: %s.\n",
                       (int) info, magma_strerror( info ));
            
            /* =====================================================================
               Check the factorization
               =================================================================== */
            if ( opts.check ) {
                magma_dgetmatrix( N, N, d_R, lda, h_R, lda );
                magma_dgetmatrix( N, N, d_R, lda, h_Q, lda );
                lapackf77_dorgtr( &opts.uplo, &N, h_Q, &lda, tau, h_work, &lwork, &info );
                
                #if defined(PRECISION_z) || defined(PRECISION_c)
                lapackf77_dsyt21( &itwo, &opts.uplo, &N, &ione,
                                  h_A, &lda, diag, offdiag,
                                  h_Q, &lda, h_R, &lda,
                                  tau, work, rwork, &result[0] );
                
                lapackf77_dsyt21( &ithree, &opts.uplo, &N, &ione,
                                  h_A, &lda, diag, offdiag,
                                  h_Q, &lda, h_R, &lda,
                                  tau, work, rwork, &result[1] );
                #else
                lapackf77_dsyt21( &itwo, &opts.uplo, &N, &ione,
                                  h_A, &lda, diag, offdiag,
                                  h_Q, &lda, h_R, &lda,
                                  tau, work, &result[0] );
                
                lapackf77_dsyt21( &ithree, &opts.uplo, &N, &ione,
                                  h_A, &lda, diag, offdiag,
                                  h_Q, &lda, h_R, &lda,
                                  tau, work, &result[1] );
                #endif
            }
            
            /* =====================================================================
               Performs operation using LAPACK
               =================================================================== */
            if ( opts.lapack ) {
                cpu_time = magma_wtime();
                lapackf77_dsytrd( &opts.uplo, &N, h_A, &lda, diag, offdiag, tau,
                                  h_work, &lwork, &info );
                cpu_time = magma_wtime() - cpu_time;
                cpu_perf = gflops / cpu_time;
                if (info != 0)
                    printf("lapackf77_dsytrd returned error %d: %s.\n",
                           (int) info, magma_strerror( info ));
            }
            
            /* =====================================================================
               Print performance and error.
               =================================================================== */
            if ( opts.lapack ) {
                printf("%5d   %7.2f (%7.2f)   %7.2f (%7.2f)",
                       (int) N, cpu_perf, cpu_time, gpu_perf, gpu_time );
            } else {
                printf("%5d     ---   (  ---  )   %7.2f (%7.2f)",
                       (int) N, gpu_perf, gpu_time );
            }
            if ( opts.check ) {
                printf("   %8.2e        %8.2e%s", result[0]*eps, result[1]*eps,
                        ( ( (result[0]*eps < tol) && (result[1]*eps < tol) ) ? "" : "  failed")  );
                status |= ! (result[0]*eps < tol);
                status |= ! (result[1]*eps < tol);
            } else {
                printf("     ---             ---\n" );
            }
            
            TESTING_FREE( h_A );
            TESTING_FREE( tau );
            TESTING_FREE( diag );
            TESTING_FREE( offdiag );
            TESTING_HOSTFREE( h_R );
            TESTING_HOSTFREE( h_work );
            TESTING_DEVFREE ( d_R );
            TESTING_DEVFREE ( dwork );
            
            if ( opts.check ) {
                TESTING_FREE( h_Q );
                TESTING_FREE( work );
                #if defined(PRECISION_z) || defined(PRECISION_c)
                TESTING_FREE( rwork );
                #endif
            }
        }
    }

    TESTING_FINALIZE();
    return status;
}
