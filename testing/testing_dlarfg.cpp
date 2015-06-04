/*
 *  -- MAGMA (version 1.3.0) --
 *     Univ. of Tennessee, Knoxville
 *     Univ. of California, Berkeley
 *     Univ. of Colorado, Denver
 *     November 2012
 *
 * @precisions normal d -> s
 *
 **/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <cuda_runtime_api.h>
#include <cublas.h>

#include "flops.h"
#include "magma.h"
#include "magma_lapack.h"
#include "testings.h"

#define PRECISION_d

#define FMULS_LARFG(n) (2*n)
#define FADDS_LARFG(n) (  n)

#define FLOPS_ZLARFG(n) (6. * FMULS_LARFG((double)n) + 2. * FADDS_LARFG((double)n) )
#define FLOPS_CLARFG(n) (6. * FMULS_LARFG((double)n) + 2. * FADDS_LARFG((double)n) )
#define FLOPS_DLARFG(n) (     FMULS_LARFG((double)n) +      FADDS_LARFG((double)n) )
#define FLOPS_SLARFG(n) (     FMULS_LARFG((double)n) +      FADDS_LARFG((double)n) )

// TODO move prototype to header
extern "C"
void magma_dlarfg( int n, double* dx0, double* dx, int incx, double* dtau );  //, double* beta );


int main( int argc, char** argv)
{
    TESTING_CUDA_INIT();

    real_Double_t   gflops, gpu_perf, gpu_time, cpu_perf, cpu_time;
    double      error, work[1];

    magma_int_t ione     = 1;
    magma_int_t ISEED[4] = {0,0,0,1};
    magma_int_t count    = 1;
    
    /* Matrix size */
    magma_int_t N = 0, size;
    const int MAXTESTS = 10;
    magma_int_t nsize[MAXTESTS] = { 1024, 2048, 3072, 4032, 5184, 6016, 7040, 8064, 9088, 10112 };
    
    double *h_x, *h_x1, *h_x2, *h_tau;
    double *d_x, *d_tau;
    double c_neg_one = MAGMA_D_NEG_ONE;

    // process command line arguments
    printf( "\nUsage: %s -N n\n"
            "  -N can be repeated up to %d times.\n\n",
            argv[0], MAXTESTS );
    
    int ntest = 0;
    for( int i = 1; i < argc; ++i ) {
        if ( strcmp("-N", argv[i]) == 0 && i+1 < argc ) {
            magma_assert( ntest < MAXTESTS, "error: -N repeated more than maximum %d tests\n", MAXTESTS );
            nsize[ ntest ] = atoi( argv[++i] );
            N = max( N, nsize[ ntest ] );
            ntest++;
        }
        else if ( strcmp("-count", argv[i]) == 0 && i+1 < argc ) {
            count = atoi( argv[++i] );
            magma_assert( count > 0, "error: -count %s is invalid; must be > 0.\n", argv[i] );
        }
        else {
            printf( "invalid argument: %s\n", argv[i] );
            exit(1);
        }
    }
    if ( ntest == 0 ) {
        ntest = MAXTESTS;
        N = nsize[ntest-1];
    }

    const int CNT = 64;
    TESTING_MALLOC( h_x,   double, N*CNT );
    TESTING_MALLOC( h_x1,  double, N*CNT );
    TESTING_MALLOC( h_x2,  double, N*CNT );
    TESTING_MALLOC( h_tau, double, CNT   );

    TESTING_DEVALLOC( d_x,   double, N*CNT );
    TESTING_DEVALLOC( d_tau, double,   CNT );
    
    magma_queue_t queue = 0;

    printf("    N    CPU GFLop/s (sec)   GPU GFlop/s (sec)   error\n");
    printf("======================================================\n");
    for( int i = 0; i < ntest; ++i ) {
    for( int cnt = 0; cnt < count; ++cnt ) {
        N = nsize[i];
        gflops = FLOPS_ZLARFG( N ) / 1e9 * CNT;

        /* Initialize the vector */
        size = N*CNT;
        lapackf77_dlarnv( &ione, ISEED, &size, h_x );
        blasf77_dcopy( &size, h_x, &ione, h_x1, &ione );
        
        /* =====================================================================
           Performs operation using MAGMA-BLAS
           =================================================================== */
        magma_dsetvector( size, h_x, ione, d_x, ione );

        gpu_time = magma_sync_wtime( queue );
        for( int j = 0; j < CNT; ++j ) {
            magma_dlarfg( N, &d_x[0+j*N], &d_x[1+j*N], ione, &d_tau[j] );
        }
        gpu_time = magma_sync_wtime( queue ) - gpu_time;
        gpu_perf = gflops / gpu_time;
        
        magma_dgetvector( size, d_x, ione, h_x2, ione );
        
        /* =====================================================================
           Performs operation using LAPACK
           =================================================================== */
        cpu_time = magma_wtime();
        for( int j = 0; j < CNT; ++j ) {
            lapackf77_dlarfg( &N, &h_x1[0+j*N], &h_x1[1+j*N], &ione, &h_tau[j] );
        }
        cpu_time = magma_wtime() - cpu_time;
        cpu_perf = gflops / cpu_time;
        
        /* =====================================================================
           Error Computation and Performance Compariosn
           =================================================================== */
        blasf77_daxpy( &size, &c_neg_one, h_x1, &ione, h_x2, &ione);
        error = lapackf77_dlange( "F", &N, &CNT, h_x2, &N, work );
        printf("%5d   %7.2f (%7.2f)   %7.2f (%7.2f)   %8.2g\n",
               (int) N, cpu_perf, cpu_time, gpu_perf, gpu_time, error );
    }
    if ( count > 1 ) {
        printf( "\n" );
    }
    }

    /* Memory clean up */
    TESTING_FREE( h_x   );
    TESTING_FREE( h_x1  );
    TESTING_FREE( h_x2  );
    TESTING_FREE( h_tau );

    TESTING_DEVFREE( d_x   );
    TESTING_DEVFREE( d_tau );

    TESTING_CUDA_FINALIZE();
    return 0;
}
