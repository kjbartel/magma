/*
    -- MAGMA (version 1.4.0-beta2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       June 2013

       @generated d Fri Jun 28 19:34:07 2013

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
   -- Testing dgebrd
*/
int main( int argc, char** argv)
{
    TESTING_INIT();

    real_Double_t    gflops, gpu_perf, gpu_time, cpu_perf, cpu_time;
    double *h_A, *h_Q, *h_PT, *h_work, *chkwork;
    double *taup, *tauq;
    #if defined(PRECISION_z) || defined(PRECISION_c)
    double      *rwork;
    #endif
    double      *diag, *offdiag;
    double      eps, result[3] = {0., 0., 0.};
    magma_int_t M, N, n2, lda, lhwork, lchkwork, info, minmn, nb;
    magma_int_t ione     = 1;
    magma_int_t ISEED[4] = {0,0,0,1};

    eps = lapackf77_dlamch( "E" );

    magma_opts opts;
    parse_opts( argc, argv, &opts );
    
    printf("    M     N   CPU GFlop/s (sec)   GPU GFlop/s (sec)   |A-QBP'|/N|A|  |I-QQ'|/N  |I-PP'|/N\n");
    printf("=========================================================================================\n");
    for( int i = 0; i < opts.ntest; ++i ) {
        for( int iter = 0; iter < opts.niter; ++iter ) {
            M = opts.msize[i];
            N = opts.nsize[i];
            minmn  = min(M, N);
            nb     = magma_get_dgebrd_nb(N);
            lda    = M;
            n2     = lda*N;
            lhwork = (M + N)*nb;
            gflops = FLOPS_DGEBRD( M, N ) / 1e9;

            TESTING_MALLOC( h_A,     double, lda*N );
            TESTING_MALLOC( tauq,    double, minmn  );
            TESTING_MALLOC( taup,    double, minmn  );
            TESTING_MALLOC( diag,    double, minmn   );
            TESTING_MALLOC( offdiag, double, (minmn-1) );
            TESTING_HOSTALLOC( h_Q, double, lda*N );
            TESTING_HOSTALLOC( h_work, double, lhwork );
            
            /* Initialize the matrices */
            lapackf77_dlarnv( &ione, ISEED, &n2, h_A );
            lapackf77_dlacpy( MagmaUpperLowerStr, &M, &N, h_A, &lda, h_Q, &lda );
            
            /* ====================================================================
               Performs operation using MAGMA
               =================================================================== */
            gpu_time = magma_wtime();
            magma_dgebrd( M, N, h_Q, lda,
                          diag, offdiag, tauq, taup,
                          h_work, lhwork, &info);
            gpu_time = magma_wtime() - gpu_time;
            gpu_perf = gflops / gpu_time;
            if (info != 0)
                printf("magma_dgebrd returned error %d: %s.\n",
                       (int) info, magma_strerror( info ));
            
            /* =====================================================================
               Check the factorization
               =================================================================== */
            if ( opts.check ) {
                lchkwork = max( minmn * nb, M+N );
                /* For optimal performance in dort01 */
                lchkwork = max( lchkwork, minmn*minmn );
                TESTING_MALLOC( h_PT,    double, lda*N   );
                TESTING_MALLOC( chkwork, double, lchkwork );
                #if defined(PRECISION_z) || defined(PRECISION_c)
                TESTING_MALLOC( rwork, double, 5*minmn );
                #endif

                lapackf77_dlacpy(MagmaUpperLowerStr, &M, &N, h_Q, &lda, h_PT, &lda);
                
                // generate Q & P'
                lapackf77_dorgbr("Q", &M, &minmn, &N, h_Q,  &lda, tauq, chkwork, &lchkwork, &info);
                if (info != 0)
                    printf("lapackf77_dorgbr returned error %d: %s.\n",
                           (int) info, magma_strerror( info ));
                lapackf77_dorgbr("P", &minmn, &N, &M, h_PT, &lda, taup, chkwork, &lchkwork, &info);
                if (info != 0)
                    printf("lapackf77_dorgbr (2) returned error %d: %s.\n",
                           (int) info, magma_strerror( info ));
                
                // Test 1:  Check the decomposition A := Q * B * PT
                //      2:  Check the orthogonality of Q
                //      3:  Check the orthogonality of PT
                #if defined(PRECISION_z) || defined(PRECISION_c)
                lapackf77_dbdt01(&M, &N, &ione,
                                 h_A, &lda, h_Q, &lda,
                                 diag, offdiag, h_PT, &lda,
                                 chkwork, rwork, &result[0]);
                lapackf77_dort01("Columns", &M, &minmn, h_Q,  &lda, chkwork, &lchkwork, rwork, &result[1]);
                lapackf77_dort01("Rows",    &minmn, &N, h_PT, &lda, chkwork, &lchkwork, rwork, &result[2]);
                #else
                lapackf77_dbdt01(&M, &N, &ione,
                                 h_A, &lda, h_Q, &lda,
                                 diag, offdiag, h_PT, &lda,
                                 chkwork, &result[0]);
                lapackf77_dort01("Columns", &M, &minmn, h_Q,  &lda, chkwork, &lchkwork, &result[1]);
                lapackf77_dort01("Rows",    &minmn, &N, h_PT, &lda, chkwork, &lchkwork, &result[2]);
                #endif
                
                TESTING_FREE( h_PT );
                TESTING_FREE( chkwork );
                #if defined(PRECISION_z) || defined(PRECISION_c)
                TESTING_FREE( rwork );
                #endif
            }
            
            /* =====================================================================
               Performs operation using LAPACK
               =================================================================== */
            if ( opts.lapack ) {
                cpu_time = magma_wtime();
                lapackf77_dgebrd(&M, &N, h_A, &lda,
                                 diag, offdiag, tauq, taup,
                                 h_work, &lhwork, &info);
                cpu_time = magma_wtime() - cpu_time;
                cpu_perf = gflops / cpu_time;
                if (info != 0)
                    printf("lapackf77_dgebrd returned error %d: %s.\n",
                           (int) info, magma_strerror( info ));
            }
            
            /* =====================================================================
               Print performance and error.
               =================================================================== */
            if ( opts.lapack ) {
                printf("%5d %5d   %7.2f (%7.2f)   %7.2f (%7.2f)",
                       (int) M, (int) N, cpu_perf, cpu_time, gpu_perf, gpu_time );
            }
            else {
                printf("%5d %5d     ---   (  ---  )   %7.2f (%7.2f)",
                       (int) M, (int) N, gpu_perf, gpu_time );
            }
            if ( opts.check ) {
                printf("   %8.2e       %8.2e   %8.2e\n",
                       result[0]*eps, result[1]*eps, result[2]*eps );
            } else {
                printf("     ---            ---      ---\n");
            }
            
            TESTING_FREE( h_A );
            TESTING_FREE( tauq );
            TESTING_FREE( taup );
            TESTING_FREE( diag );
            TESTING_FREE( offdiag );
            TESTING_HOSTFREE( h_Q );
            TESTING_HOSTFREE( h_work );
        }
        if ( opts.niter > 1 ) {
            printf( "\n" );
        }
    }

    TESTING_FINALIZE();
    return 0;
}
