/*
 *  -- MAGMA (version 1.2.1) --
 *     Univ. of Tennessee, Knoxville
 *     Univ. of California, Berkeley
 *     Univ. of Colorado, Denver
 *     June 2012
 *
 * @generated d Thu Jun 28 12:31:43 2012
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

// Flops formula
#define PRECISION_d
#if defined(PRECISION_z) || defined(PRECISION_c)
#define FLOPS(m, n) ( 6. * FMULS_GETRF(m, n) + 2. * FADDS_GETRF(m, n) )
#else
#define FLOPS(m, n) (      FMULS_GETRF(m, n) +      FADDS_GETRF(m, n) )
#endif

double get_LU_error(magma_int_t M, magma_int_t N, 
                    double *A,  magma_int_t lda, 
                    double *LU, magma_int_t *IPIV)
{
    magma_int_t min_mn = min(M,N);
    magma_int_t ione   = 1;
    magma_int_t i, j;
    double alpha = MAGMA_D_ONE;
    double beta  = MAGMA_D_ZERO;
    double *L, *U;
    double work[1], matnorm, residual;
                       
    TESTING_MALLOC( L, double, M*min_mn);
    TESTING_MALLOC( U, double, min_mn*N);
    memset( L, 0, M*min_mn*sizeof(double) );
    memset( U, 0, min_mn*N*sizeof(double) );

    lapackf77_dlaswp( &N, A, &lda, &ione, &min_mn, IPIV, &ione);
    lapackf77_dlacpy( MagmaLowerStr, &M, &min_mn, LU, &lda, L, &M      );
    lapackf77_dlacpy( MagmaUpperStr, &min_mn, &N, LU, &lda, U, &min_mn );

    for(j=0; j<min_mn; j++)
        L[j+j*M] = MAGMA_D_MAKE( 1., 0. );
    
    matnorm = lapackf77_dlange("f", &M, &N, A, &lda, work);

    blasf77_dgemm("N", "N", &M, &N, &min_mn,
                  &alpha, L, &M, U, &min_mn, &beta, LU, &lda);

    for( j = 0; j < N; j++ ) {
        for( i = 0; i < M; i++ ) {
            LU[i+j*lda] = MAGMA_D_SUB( LU[i+j*lda], A[i+j*lda] );
        }
    }
    residual = lapackf77_dlange("f", &M, &N, LU, &lda, work);

    TESTING_FREE(L);
    TESTING_FREE(U);

    return residual / (matnorm * N);
}

/* ////////////////////////////////////////////////////////////////////////////
   -- Testing dgetrf
*/
int main( int argc, char** argv)
{
    TESTING_CUDA_INIT();

    magma_timestr_t       start, end;
    double           flops, gpu_perf, cpu_perf, error;
    double *h_A, *h_R;
    double *d_A;
    magma_int_t     *ipiv;

    /* Matrix size */
    magma_int_t M = 0, N = 0, n2, lda, ldda;
    magma_int_t size[10] = {960,1920,3072,4032,4992,5952,7104,8064,9024,9984};

    magma_int_t i, info, min_mn, nb, maxn;
    magma_int_t ione     = 1;
    magma_int_t ISEED[4] = {0,0,0,1};

    if (argc != 1){
        for(i = 1; i<argc; i++){
            if (strcmp("-N", argv[i])==0)
                N = atoi(argv[++i]);
            else if (strcmp("-M", argv[i])==0)
                M = atoi(argv[++i]);
        }
        if (M>0 && N>0)
            printf("  testing_dgetrf -M %d -N %d\n\n", (int) M, (int) N);
        else
            {
                printf("\nUsage: \n");
                printf("  testing_dgetrf -M %d -N %d\n\n", 1024, 1024);
                exit(1);
            }
    }
    else {
        printf("\nUsage: \n");
        printf("  testing_dgetrf_gpu -M %d -N %d\n\n", 1024, 1024);
        M = N = size[9];
    }

    ldda   = ((M+31)/32)*32;
    maxn   = ((N+31)/32)*32;
    n2     = M * N;
    min_mn = min(M, N);
    nb     = magma_get_dgetrf_nb(min_mn);

    /* Allocate host memory for the matrix */
    TESTING_MALLOC(ipiv, magma_int_t, min_mn);
    TESTING_MALLOC(    h_A, double, n2     );
    TESTING_HOSTALLOC( h_R, double, n2     );
    TESTING_DEVALLOC(  d_A, double, ldda*N );

    printf("  M     N   CPU GFlop/s    GPU GFlop/s   ||PA-LU||/(||A||*N)\n");
    printf("============================================================\n");
    for(i=0; i<10; i++){
        if (argc == 1){
            M = N = size[i];
        }
        min_mn= min(M, N);
        lda   = M;
        n2    = lda*N;
        ldda  = ((M+31)/32)*32;
        flops = FLOPS( (double)M, (double)N ) / 1000000;

        /* Initialize the matrix */
        lapackf77_dlarnv( &ione, ISEED, &n2, h_A );
        lapackf77_dlacpy( MagmaUpperLowerStr, &M, &N, h_A, &lda, h_R, &lda );

        /* =====================================================================
           Performs operation using LAPACK
           =================================================================== */
        start = get_current_time();
        lapackf77_dgetrf(&M, &N, h_A, &lda, ipiv, &info);
        end = get_current_time();
        if (info < 0)
            printf("Argument %d of dgetrf had an illegal value.\n", (int) -info);

        cpu_perf = flops / GetTimerValue(start, end);

        /* ====================================================================
           Performs operation using MAGMA
           =================================================================== */
        magma_dsetmatrix( M, N, h_R, lda, d_A, ldda );
        start = get_current_time();
        magma_dgetrf_gpu( M, N, d_A, ldda, ipiv, &info);
        end = get_current_time();
        if (info < 0)
            printf("Argument %d of dgetrf had an illegal value.\n", (int) -info);

        gpu_perf = flops / GetTimerValue(start, end);

        /* =====================================================================
           Check the factorization
           =================================================================== */
        magma_dgetmatrix( M, N, d_A, ldda, h_A, lda );
        error = get_LU_error(M, N, h_R, lda, h_A, ipiv);
        
        printf("%5d %5d  %6.2f         %6.2f         %e\n",
               (int) M, (int) N, cpu_perf, gpu_perf, error);

        if (argc != 1)
            break;
    }

    /* Memory clean up */
    TESTING_FREE( ipiv );
    TESTING_FREE( h_A );
    TESTING_HOSTFREE( h_R );
    TESTING_DEVFREE( d_A );

    /* Shutdown */
    TESTING_CUDA_FINALIZE();
}
