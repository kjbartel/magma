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

/* ////////////////////////////////////////////////////////////////////////////
   -- Testing sgetrf
*/
int main( int argc, char** argv)
{
    TESTING_CUDA_INIT();

    magma_timestr_t  start, end;
    float      flops, gpu_perf, cpu_perf;
    float *h_A, *h_R;
    float *d_A, *dwork;
    magma_int_t N = 0, n2, lda, ldda;
    magma_int_t size[10] = { 1024, 2048, 3072, 4032, 5184, 6016, 7040, 8064, 9088, 10112 };
    magma_int_t ntest = 10;
    
    magma_int_t i, info;
    float c_neg_one = MAGMA_S_NEG_ONE;
    magma_int_t ione     = 1;
    magma_int_t ISEED[4] = {0, 0, 0, 1};
    float *work;
    float tmp;
    float rwork[1];
    magma_int_t *ipiv;
    magma_int_t lwork, ldwork;
    float A_norm, R_norm;
    
    if (argc != 1){
        for(i = 1; i<argc; i++){
            if (strcmp("-N", argv[i])==0)
                N = atoi(argv[++i]);
        }
        if (N>0) size[0] = size[ntest-1] = N;
        else exit(1);
    }
    else {
        printf("\nUsage: \n");
        printf("  testing_sgetri_gpu -N %d\n\n", 1024);
    }
    
    /* query for Lapack workspace size */
    N     = size[ntest-1];
    lda   = N;
    work  = &tmp;
    lwork = -1;
    lapackf77_sgetri( &N, h_A, &lda, ipiv, work, &lwork, &info );
    if (info != 0)
        printf("lapackf77_sgetri returned error %d\n", (int) info);
    lwork = int( MAGMA_S_REAL( *work ));

    /* query for Magma workspace size */
    ldwork = N * magma_get_sgetri_nb( N );

    /* Allocate memory */
    n2   = N * N;
    ldda = ((N+31)/32) * 32;
    TESTING_MALLOC(    ipiv,  magma_int_t,     N      );
    TESTING_MALLOC(    work,  float, lwork  );
    TESTING_MALLOC(    h_A,   float, n2     );
    TESTING_HOSTALLOC( h_R,   float, n2     );
    TESTING_DEVALLOC(  d_A,   float, ldda*N );
    TESTING_DEVALLOC(  dwork, float, ldwork );

    printf("  N    CPU GFlop/s    GPU GFlop/s    ||R||_F / ||A||_F\n");
    printf("========================================================\n");
    for( i=0; i < ntest; i++ ){
        N   = size[i];
        lda = N;
        n2  = lda*N;
        flops = FLOPS_SGETRI( (float)N ) / 1000000;
        
        ldda = ((N+31)/32)*32;

        /* Initialize the matrix */
        lapackf77_slarnv( &ione, ISEED, &n2, h_A );
        A_norm = lapackf77_slange( "f", &N, &N, h_A, &lda, rwork );

        /* Factor the matrix. Both MAGMA and LAPACK will use this factor. */
        magma_ssetmatrix( N, N, h_A, lda, d_A, ldda );
        magma_sgetrf_gpu( N, N, d_A, ldda, ipiv, &info );
        magma_sgetmatrix( N, N, d_A, ldda, h_A, lda );
        
        // check for exact singularity
        //h_A[ 10 + 10*lda ] = MAGMA_S_MAKE( 0.0, 0.0 );
        //magma_ssetmatrix( N, N, h_A, lda, d_A, ldda );

        /* ====================================================================
           Performs operation using MAGMA
           =================================================================== */
        start = get_current_time();
        magma_sgetri_gpu( N,    d_A, ldda, ipiv, dwork, ldwork, &info );
        end = get_current_time();
        if (info != 0)
            printf("magma_sgetri_gpu returned error %d\n", (int) info);

        gpu_perf = flops / GetTimerValue(start, end);
        
        magma_sgetmatrix( N, N, d_A, ldda, h_R, lda );
         
        /* =====================================================================
           Performs operation using LAPACK
           =================================================================== */
        start = get_current_time();
        lapackf77_sgetri( &N,     h_A, &lda, ipiv, work, &lwork, &info );
        end = get_current_time();
        if (info != 0)
            printf("lapackf77_sgetri returned error %d\n", (int) info);
        
        cpu_perf = flops / GetTimerValue(start, end);
        
        /* =====================================================================
           Check the result compared to LAPACK
           =================================================================== */
        blasf77_saxpy( &n2, &c_neg_one, h_A, &ione, h_R, &ione );
        R_norm = lapackf77_slange( "f", &N, &N, h_R, &lda, rwork );
        
        printf( "%5d    %6.2f         %6.2f        %e\n",
                (int) N, cpu_perf, gpu_perf, R_norm / A_norm );
        
        if (argc != 1)
            break;
    }

    /* Memory clean up */
    TESTING_FREE(     ipiv  );
    TESTING_FREE(     work  );
    TESTING_FREE(     h_A   );
    TESTING_HOSTFREE( h_R   );
    TESTING_DEVFREE(  d_A   );
    TESTING_DEVFREE(  dwork );

    /* Shutdown */
    TESTING_CUDA_FINALIZE();
}
