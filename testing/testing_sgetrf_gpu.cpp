/*
 *  -- MAGMA (version 1.3.0) --
 *     Univ. of Tennessee, Knoxville
 *     Univ. of California, Berkeley
 *     Univ. of Colorado, Denver
 *     November 2012
 *
 * @generated s Wed Nov 14 22:54:16 2012
 *
 **/
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


float get_LU_error(magma_int_t M, magma_int_t N, 
                    float *A,  magma_int_t lda, 
                    float *LU, magma_int_t *IPIV)
{
    magma_int_t min_mn = min(M,N);
    magma_int_t ione   = 1;
    magma_int_t i, j;
    float alpha = MAGMA_S_ONE;
    float beta  = MAGMA_S_ZERO;
    float *L, *U;
    float work[1], matnorm, residual;
                       
    TESTING_MALLOC( L, float, M*min_mn);
    TESTING_MALLOC( U, float, min_mn*N);
    memset( L, 0, M*min_mn*sizeof(float) );
    memset( U, 0, min_mn*N*sizeof(float) );

    lapackf77_slaswp( &N, A, &lda, &ione, &min_mn, IPIV, &ione);
    lapackf77_slacpy( MagmaLowerStr, &M, &min_mn, LU, &lda, L, &M      );
    lapackf77_slacpy( MagmaUpperStr, &min_mn, &N, LU, &lda, U, &min_mn );

    for(j=0; j<min_mn; j++)
        L[j+j*M] = MAGMA_S_MAKE( 1., 0. );
    
    matnorm = lapackf77_slange("f", &M, &N, A, &lda, work);

    blasf77_sgemm("N", "N", &M, &N, &min_mn,
                  &alpha, L, &M, U, &min_mn, &beta, LU, &lda);

    for( j = 0; j < N; j++ ) {
        for( i = 0; i < M; i++ ) {
            LU[i+j*lda] = MAGMA_S_SUB( LU[i+j*lda], A[i+j*lda] );
        }
    }
    residual = lapackf77_slange("f", &M, &N, LU, &lda, work);

    TESTING_FREE(L);
    TESTING_FREE(U);

    return residual / (matnorm * N);
}

/* ////////////////////////////////////////////////////////////////////////////
   -- Testing sgetrf
*/
int main( int argc, char** argv)
{
    TESTING_CUDA_INIT();

    real_Double_t   gflops, gpu_perf, gpu_time, cpu_perf, cpu_time;
    float          error;
    float *h_A, *h_R;
    float *d_A;
    magma_int_t     *ipiv;

    /* Matrix size */
    magma_int_t M = 0, N = 0, n2, lda, ldda;
    const int MAXTESTS = 10;
    magma_int_t msize[MAXTESTS] = { 1024, 2048, 3072, 4032, 5184, 6016, 7040, 8064, 9088, 10112 };
    magma_int_t nsize[MAXTESTS] = { 1024, 2048, 3072, 4032, 5184, 6016, 7040, 8064, 9088, 10112 };

    magma_int_t i, info, min_mn, nb;
    magma_int_t ione     = 1;
    magma_int_t ISEED[4] = {0,0,0,1};
    magma_int_t checkres;

    checkres = getenv("MAGMA_TESTINGS_CHECK") != NULL;

    // process command line arguments
    printf( "\nUsage: %s -N <m,n> -c\n", argv[0] );
    printf( "  -N can be repeated up to %d times. If only m is given, then m=n.\n", MAXTESTS );
    printf( "  -c or setting $MAGMA_TESTINGS_CHECK runs LAPACK and checks result.\n\n" );
    int ntest = 0;
    for( int i = 1; i < argc; ++i ) {
        if ( strcmp("-N", argv[i]) == 0 && i+1 < argc ) {
            magma_assert( ntest < MAXTESTS, "error: -N repeated more than maximum %d tests\n", MAXTESTS );
            int m, n;
            info = sscanf( argv[++i], "%d,%d", &m, &n );
            if ( info == 2 && m > 0 && n > 0 ) {
                msize[ ntest ] = m;
                nsize[ ntest ] = n;
            }
            else if ( info == 1 && m > 0 ) {
                msize[ ntest ] = m;
                nsize[ ntest ] = m;  // implicitly
            }
            else {
                printf( "error: -N %s is invalid; ensure m > 0, n > 0.\n", argv[i] );
                exit(1);
            }
            M = max( M, msize[ ntest ] );
            N = max( N, nsize[ ntest ] );
            ntest++;
        }
        else if ( strcmp("-M", argv[i]) == 0 ) {
            printf( "-M has been replaced in favor of -N m,n to allow -N to be repeated.\n\n" );
            exit(1);
        }
        else if ( strcmp("-c", argv[i]) == 0 ) {
            checkres = true;
        }
        else {
            printf( "invalid argument: %s\n", argv[i] );
            exit(1);
        }
    }
    if ( ntest == 0 ) {
        ntest = MAXTESTS;
        M = msize[ntest-1];
        N = nsize[ntest-1];
    }

    ldda   = ((M+31)/32)*32;
    n2     = M * N;
    min_mn = min(M, N);
    nb     = magma_get_sgetrf_nb(min_mn);

    /* Allocate memory for the matrix */
    TESTING_MALLOC(ipiv, magma_int_t, min_mn);
    TESTING_MALLOC(    h_A, float, n2     );
    TESTING_HOSTALLOC( h_R, float, n2     );
    TESTING_DEVALLOC(  d_A, float, ldda*N );

    printf("  M     N     CPU GFlop/s (sec)   GPU GFlop/s (sec)   ||PA-LU||/(||A||*N)\n");
    printf("=========================================================================\n");
    for( i = 0; i < ntest; ++i ) {
        M = msize[i];
        N = nsize[i];
        min_mn= min(M, N);
        lda   = M;
        n2    = lda*N;
        ldda  = ((M+31)/32)*32;
        gflops = FLOPS_SGETRF( M, N ) / 1e9;

        /* Initialize the matrix */
        if ( checkres ) {
            lapackf77_slarnv( &ione, ISEED, &n2, h_A );
            lapackf77_slacpy( MagmaUpperLowerStr, &M, &N, h_A, &lda, h_R, &lda );
        }
        else
            lapackf77_slarnv( &ione, ISEED, &n2, h_R );
        magma_ssetmatrix( M, N, h_R, lda, d_A, ldda );

        /* =====================================================================
           Performs operation using LAPACK
           =================================================================== */
        if ( checkres ) {
            cpu_time = magma_wtime();
            lapackf77_sgetrf(&M, &N, h_A, &lda, ipiv, &info);
            cpu_time = magma_wtime() - cpu_time;
            cpu_perf = gflops / cpu_time;
            if (info != 0)
                printf("lapackf77_sgetrf returned error %d.\n", (int) info);    
        }
        
        /* ====================================================================
           Performs operation using MAGMA
           =================================================================== */
        gpu_time = magma_wtime();
        magma_sgetrf_gpu( M, N, d_A, ldda, ipiv, &info);
        gpu_time = magma_wtime() - gpu_time;
        gpu_perf = gflops / gpu_time;
        if (info != 0)
            printf("magma_sgetrf_gpu returned error %d.\n", (int) info);

        /* =====================================================================
           Check the factorization
           =================================================================== */
        if ( checkres ) {
            magma_sgetmatrix( M, N, d_A, ldda, h_A, lda );
            error = get_LU_error(M, N, h_R, lda, h_A, ipiv);
            
            printf("%5d %5d   %7.2f (%7.2f)   %7.2f (%7.2f)   %8.2e\n",
                   (int) M, (int) N, cpu_perf, cpu_time, gpu_perf, gpu_time, error);
        }
        else {
            printf("%5d %5d     ---   (  ---  )   %7.2f (%7.2f)     ---  \n",
                   (int) M, (int) N, gpu_perf, gpu_time);
        }
    }

    /* Memory clean up */
    TESTING_FREE( ipiv );
    TESTING_FREE( h_A );
    TESTING_HOSTFREE( h_R );
    TESTING_DEVFREE( d_A );

    TESTING_CUDA_FINALIZE();
    return 0;
}
