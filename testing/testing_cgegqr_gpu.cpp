/*
    -- MAGMA (version 1.4.1-beta2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       December 2013

       @generated c Mon Dec  9 16:19:29 2013
       @author Stan Tomov

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

extern "C" magma_int_t
magma_cgegqr_gpu( magma_int_t m, magma_int_t n,
                  magmaFloatComplex *dA,   magma_int_t ldda,
                  magmaFloatComplex *dwork, magmaFloatComplex *work,
                  magma_int_t *info );


/* ////////////////////////////////////////////////////////////////////////////
   -- Testing cgegqr
*/
int main( int argc, char** argv)
{
    TESTING_INIT();

    real_Double_t    gflops, gpu_perf, gpu_time, cpu_perf, cpu_time;
    float           e1, e2, work[1];
    magmaFloatComplex *h_A, *h_R, *tau, *dtau, *h_work, tmp[1];
    magmaFloatComplex *d_A, *dwork, *ddA, *d_T;
    magma_int_t M, N, n2, lda, ldda, lwork, info, min_mn;
    magma_int_t ione     = 1;
    magma_int_t ISEED[4] = {0,0,0,1};

    magma_opts opts;
    parse_opts( argc, argv, &opts );
    opts.lapack |= opts.check;  // check (-c) implies lapack (-l)
    
    printf("  M     N     CPU GFlop/s (ms)    GPU GFlop/s (ms)    ||I - Q'Q||_F    \n");
    printf("=======================================================================\n");
    for( int i = 0; i < opts.ntest; ++i ) {
        for( int iter = 0; iter < opts.niter; ++iter ) {
            M = opts.msize[i];
            N = opts.nsize[i];
            min_mn = min(M, N);
            lda    = M;
            n2     = lda*N;
            ldda   = ((M+31)/32)*32;
            gflops = FLOPS_CGEQRF( M, N ) / 1e9 +  FLOPS_CUNGQR( M, N, N ) / 1e9;
            
            // query for workspace size
            lwork = -1;
            lapackf77_cgeqrf(&M, &N, NULL, &M, NULL, tmp, &lwork, &info);
            lwork = (magma_int_t)MAGMA_C_REAL( tmp[0] );
            lwork = max(lwork, 3*N*N);
            
            TESTING_MALLOC_PIN( tau,    magmaFloatComplex, min_mn );
            TESTING_MALLOC_PIN( h_work, magmaFloatComplex, lwork  );
            
            TESTING_MALLOC_CPU( h_A,   magmaFloatComplex, n2     );
            TESTING_MALLOC_CPU( h_R,   magmaFloatComplex, n2     );
            
            TESTING_MALLOC_DEV( d_A,   magmaFloatComplex, ldda*N );
            TESTING_MALLOC_DEV( dtau,  magmaFloatComplex, min_mn );
            TESTING_MALLOC_DEV( dwork, magmaFloatComplex, N*N    );
            TESTING_MALLOC_DEV( ddA,   magmaFloatComplex, N*N    );
            TESTING_MALLOC_DEV( d_T,   magmaFloatComplex, N*N    );
            
            cudaMemset( ddA, 0, N*N*sizeof(magmaFloatComplex) );
            cudaMemset( d_T, 0, N*N*sizeof(magmaFloatComplex) );

            /* Initialize the matrix */
            lapackf77_clarnv( &ione, ISEED, &n2, h_A );
            lapackf77_clacpy( MagmaUpperLowerStr, &M, &N, h_A, &lda, h_R, &lda );
            magma_csetmatrix( M, N, h_R, lda, d_A, ldda );
            
            // warmup
            magma_cgegqr_gpu( M, N, d_A, ldda, dwork, h_work, &info );
            magma_csetmatrix( M, N, h_R, lda, d_A, ldda );
            
            /* ====================================================================
               Performs operation using MAGMA
               =================================================================== */
            gpu_time = magma_sync_wtime( 0 );
            if (opts.version == 2) {
                int min_mn = min(M, N);
                int     nb = N;

                cuFloatComplex *dtau = dwork;
                
                magma_cgeqr2x3_gpu(&M, &N, d_A, &ldda, dtau, d_T, ddA, 
                                   (float *)(dwork+min_mn), &info);
                magma_cgetmatrix( min_mn, 1, dtau, min_mn, tau, min_mn);  
                magma_cungqr_gpu( M, N, N, d_A, ldda, tau, d_T, nb, &info );
            }
            else
               magma_cgegqr_gpu( M, N, d_A, ldda, dwork, h_work, &info );
            gpu_time = magma_sync_wtime( 0 ) - gpu_time;

            gpu_perf = gflops / gpu_time;
            if (info != 0)
                printf("magma_cgegqr returned error %d: %s.\n",
                       (int) info, magma_strerror( info ));
            
            if ( opts.lapack ) {
                /* =====================================================================
                   Performs operation using LAPACK
                   =================================================================== */
                cpu_time = magma_wtime();

                /* Orthogonalize on the CPU */
                lapackf77_cgeqrf(&M, &N, h_A, &lda, tau, h_work, &lwork, &info);
                lapackf77_cungqr(&M, &N, &N, h_A, &lda, tau, h_work, &lwork, &info );

                cpu_time = magma_wtime() - cpu_time;
                cpu_perf = gflops / cpu_time;
                if (info != 0)
                    printf("lapackf77_cungqr returned error %d: %s.\n",
                           (int) info, magma_strerror( info ));
                
                /* =====================================================================
                   Check the result compared to LAPACK
                   =================================================================== */
                magma_cgetmatrix( M, N, d_A, ldda, h_R, M );

                magmaFloatComplex one = MAGMA_C_ONE, zero = MAGMA_C_ZERO;
                blasf77_cgemm("t", "n", &N, &N, &M, &one, h_R, &M, h_R, &M, &zero, h_work, &N);
                for(int ii=0; ii<N*N; ii+=(N+1)) h_work[ii] = MAGMA_C_SUB(h_work[ii], one);

                e1    = lapackf77_clange("f", &N, &N, h_work, &N, work);

                blasf77_cgemm("t", "n", &N, &N, &M, &one, h_A, &M, h_A, &M, &zero, h_work, &N);
                for(int ii=0; ii<N*N; ii+=(N+1)) h_work[ii] = MAGMA_C_SUB(h_work[ii], one);
                e2    = lapackf77_clange("f", &N, &N, h_work, &N, work);
                
                printf("%5d %5d   %7.2f (%7.2f)   %7.2f (%7.2f)   %8.2e  %8.2e\n",
                       (int) M, (int) N, cpu_perf, 1000.*cpu_time, gpu_perf, 1000.*gpu_time, e1, e2 );
            }
            else {
                printf("%5d %5d     ---   (  ---  )   %7.2f (%7.2f)     ---  \n",
                       (int) M, (int) N, gpu_perf, 1000.*gpu_time );
            }
            
            TESTING_FREE_PIN( tau    );
            TESTING_FREE_PIN( h_work );
            
            TESTING_FREE_CPU( h_A  );
            TESTING_FREE_CPU( h_R  );
            
            TESTING_FREE_DEV( d_A   );
            TESTING_FREE_DEV( dtau  );
            TESTING_FREE_DEV( dwork );
            TESTING_FREE_DEV( ddA   );
            TESTING_FREE_DEV( d_T   );
        }
        if ( opts.niter > 1 ) {
            printf( "\n" );
        }
    }
    
    TESTING_FINALIZE();
    return 0;
}
