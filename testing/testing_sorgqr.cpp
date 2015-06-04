/*
    -- MAGMA (version 1.4.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       August 2013

       @generated s Wed Aug 14 12:18:07 2013

       @author Stan Tomov
       @author Mathieu Faverge
       @author Mark Gates
*/

// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <cuda_runtime_api.h>
#include <cublas.h>
#include <assert.h>

// includes, project
#include "flops.h"
#include "magma.h"
#include "magma_lapack.h"
#include "testings.h"


/* ////////////////////////////////////////////////////////////////////////////
   -- Testing sorgqr
*/
int main( int argc, char** argv )
{
    TESTING_INIT();

    real_Double_t    gflops, gpu_perf, gpu_time, cpu_perf, cpu_time;
    float           error, work[1];
    float  c_neg_one = MAGMA_S_NEG_ONE;
    float *hA, *hR, *tau, *h_work;
    float *dA, *dT;
    magma_int_t m, n, k;
    magma_int_t n2, lda, ldda, lwork, min_mn, nb, info;
    magma_int_t ione     = 1;
    magma_int_t ISEED[4] = {0,0,0,1};
    
    magma_opts opts;
    parse_opts( argc, argv, &opts );
    opts.lapack |= opts.check;  // check (-c) implies lapack (-l)
    
    printf("Running version %d; available are (specified through --version num):\n",
           (int) opts.version);
    printf("1 - uses precomputed slarft matrices (default)\n");
    printf("2 - recomputes the slarft matrices on the fly\n\n");

    printf("    m     n     k   CPU GFlop/s (sec)   GPU GFlop/s (sec)   ||R|| / ||A||\n");
    printf("=========================================================================\n");
    for( int i = 0; i < opts.ntest; ++i ){
        for( int iter = 0; iter < opts.niter; ++iter ) {
            m = opts.msize[i];
            n = opts.nsize[i];
            k = opts.ksize[i];
            if ( m < n || n < k ) {
                printf( "skipping m %d, n %d, k %d because m < n or n < k\n", (int) m, (int) n, (int) k );
                continue;
            }
            
            lda  = m;
            ldda = ((m + 31)/32)*32;
            n2 = lda*n;
            min_mn = min(m, n);
            nb = magma_get_sgeqrf_nb( m );
            lwork  = (m + 2*n+nb)*nb;
            gflops = FLOPS_SORGQR( m, n, k ) / 1e9;
            
            TESTING_MALLOC(    hA,     float, lda*n  );
            TESTING_HOSTALLOC( h_work, float, lwork  );
            TESTING_HOSTALLOC( hR,     float, lda*n  );
            TESTING_MALLOC(    tau,    float, min_mn );
            TESTING_DEVALLOC(  dA,     float, ldda*n );
            TESTING_DEVALLOC(  dT,     float, ( 2*min_mn + ((n + 31)/32)*32 )*nb );
            
            lapackf77_slarnv( &ione, ISEED, &n2, hA );
            lapackf77_slacpy( MagmaUpperLowerStr, &m, &n, hA, &lda, hR, &lda );
            
            /* ====================================================================
               Performs operation using MAGMA
               =================================================================== */
            // first, get QR factors
            magma_ssetmatrix( m, n, hA, lda, dA, ldda );
            magma_sgeqrf_gpu( m, n, dA, ldda, tau, dT, &info );
            if (info != 0)
                printf("magma_sgeqrf_gpu returned error %d: %s.\n",
                       (int) info, magma_strerror( info ));
            magma_sgetmatrix( m, n, dA, ldda, hR, lda );
            
            gpu_time = magma_wtime();
            if (opts.version == 1)
                magma_sorgqr( m, n, k, hR, lda, tau, dT, nb, &info );
            else
                magma_sorgqr2(m, n, k, hR, lda, tau, &info );
            gpu_time = magma_wtime() - gpu_time;
            gpu_perf = gflops / gpu_time;
            if (info != 0)
                printf("magma_sorgqr_gpu returned error %d: %s.\n",
                       (int) info, magma_strerror( info ));
            
            /* =====================================================================
               Performs operation using LAPACK
               =================================================================== */
            if ( opts.lapack ) {
                error = lapackf77_slange("f", &m, &n, hA, &lda, work );
                
                lapackf77_sgeqrf( &m, &n, hA, &lda, tau, h_work, &lwork, &info );
                if (info != 0)
                    printf("lapackf77_sgeqrf returned error %d: %s.\n",
                           (int) info, magma_strerror( info ));
                
                cpu_time = magma_wtime();
                lapackf77_sorgqr( &m, &n, &k, hA, &lda, tau, h_work, &lwork, &info );
                cpu_time = magma_wtime() - cpu_time;
                cpu_perf = gflops / cpu_time;
                if (info != 0)
                    printf("lapackf77_sorgqr returned error %d: %s.\n",
                           (int) info, magma_strerror( info ));
                
                // compute relative error |R|/|A| := |Q_magma - Q_lapack|/|A|
                blasf77_saxpy( &n2, &c_neg_one, hA, &ione, hR, &ione );
                error = lapackf77_slange("f", &m, &n, hR, &lda, work) / error;
                
                printf("%5d %5d %5d   %7.1f (%7.2f)   %7.1f (%7.2f)   %8.2e\n",
                       (int) m, (int) n, (int) k,
                       cpu_perf, cpu_time, gpu_perf, gpu_time, error );
            }
            else {
                printf("%5d %5d %5d     ---   (  ---  )   %7.1f (%7.2f)     ---  \n",
                       (int) m, (int) n, (int) k,
                       gpu_perf, gpu_time );
            }
            
            TESTING_FREE(     hA     );
            TESTING_HOSTFREE( h_work );
            TESTING_HOSTFREE( hR     );
            TESTING_FREE(     tau    );
            TESTING_DEVFREE(  dA     );
            TESTING_DEVFREE(  dT     );
        }
        if ( opts.niter > 1 ) {
            printf( "\n" );
        }
    }
    
    TESTING_FINALIZE();
    return 0;
}
