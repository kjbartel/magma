/*
 *  -- MAGMA (version 1.2.0) --
 *     Univ. of Tennessee, Knoxville
 *     Univ. of California, Berkeley
 *     Univ. of Colorado, Denver
 *     May 2012
 *
 * @generated s Tue May 15 18:18:16 2012
 *
 **/
// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cublas.h>

// includes, project
#include "flops.h"
#include "magma.h"
#include "magma_lapack.h"
#include "testings.h"

// Flops formula
#define PRECISION_s
#if defined(PRECISION_z) || defined(PRECISION_c)
#define FLOPS(n) ( 6. * FMULS_POTRF(n) + 2. * FADDS_POTRF(n) )
#else
#define FLOPS(n) (      FMULS_POTRF(n) +      FADDS_POTRF(n) )
#endif

/* ////////////////////////////////////////////////////////////////////////////
   -- Testing spotrf
*/
int main( int argc, char** argv)
{
    TESTING_CUDA_INIT();

    magma_timestr_t       start, end;
    float           flops, gpu_perf, cpu_perf;
    float *h_A, *h_R;
    magma_int_t      N=0, n2, lda;
    magma_int_t      size[13] = {1024,2048,3072,4032,5184,6048,7200,8064,8928,10240,20000,30000,40000};
    magma_int_t      size_n = 13;

    magma_int_t  i, info, flag = 0;
    const char  *uplo     = MagmaLowerStr;
    float c_neg_one = MAGMA_S_NEG_ONE;
    magma_int_t  ione     = 1;
    magma_int_t  ISEED[4] = {0,0,0,1};
    float       work[1];  // not referenced for lange norm 'f'
    float       matnorm;

    if (argc != 1){
        for(i = 1; i<argc; i++){
            if (strcmp("-N", argv[i])==0) {
                flag = 1;
                size[0] = size[size_n-1] = atoi(argv[++i]);
             } if (strcmp("-UPLO", argv[i])==0)
              if (strcmp("U", argv[++i])) uplo = MagmaUpperStr;
        }
    }
    N = size[size_n-1];
    if (N<=0) {
      printf( "N=%d\n",N );
      exit(1);
    } else {
        printf("\nUsage: \n");
        if (strcmp(MagmaLowerStr,uplo) )
        printf("  testing_spotrf -UPLO L -N %d:%d\n\n", size[0],size[size_n-1]);
        else
        printf("  testing_spotrf -UPLO U -N %d:%d\n\n", size[0],size[size_n-1]);
    }

    /* Allocate host memory for the matrix */
    n2 = N * N;
    TESTING_MALLOC(    h_A, float, n2);
    TESTING_HOSTALLOC( h_R, float, n2);

    printf("\n\n");
    printf("  N    CPU GFlop/s    GPU GFlop/s    ||R_magma - R_lapack||_F / ||R_lapack||_F\n");
    printf("========================================================\n");
    for(i=0; i<size_n; i++){
        N     = size[i];
        lda   = N;
        n2    = lda*N;
        flops = FLOPS( (float)N ) / 1000000;

        /* ====================================================================
           Initialize the matrix
           =================================================================== */
        lapackf77_slarnv( &ione, ISEED, &n2, h_A );
        /* Symmetrize and increase the diagonal */
        {
            magma_int_t i, j;
            for(i=0; i<N; i++) {
                MAGMA_S_SET2REAL( h_A[i*lda+i], ( MAGMA_S_REAL(h_A[i*lda+i]) + 1.*N ) );
                for(j=0; j<i; j++)
                    h_A[i*lda+j] = (h_A[j*lda+i]);
            }
        }
        lapackf77_slacpy( MagmaUpperLowerStr, &N, &N, h_A, &lda, h_R, &lda );

        /* ====================================================================
           Performs operation using MAGMA
           =================================================================== */
        //magma_spotrf(uplo[0], N, h_R, lda, &info);
        //lapackf77_slacpy( MagmaUpperLowerStr, &N, &N, h_A, &lda, h_R, &lda );

        start = get_current_time();
        magma_spotrf(uplo[0], N, h_R, lda, &info);
        end = get_current_time();
        if (info != 0)
            printf("Argument %d of magma_spotrf had an illegal value.\n", -info);

        gpu_perf = flops / GetTimerValue(start, end);

        /* =====================================================================
           Performs operation using LAPACK
           =================================================================== */
        start = get_current_time();
        lapackf77_spotrf(uplo, &N, h_A, &lda, &info);
        end = get_current_time();
        if (info != 0)
            printf("Argument %d of lapack_spotrf had an illegal value.\n", -info);

        cpu_perf = flops / GetTimerValue(start, end);

        /* =====================================================================
           Check the result compared to LAPACK
           =================================================================== */
        matnorm = lapackf77_slange("f", &N, &N, h_A, &N, work);
        blasf77_saxpy(&n2, &c_neg_one, h_A, &ione, h_R, &ione);
        printf("%5d    %6.2f         %6.2f        %e\n",
               size[i], cpu_perf, gpu_perf,
               lapackf77_slange("f", &N, &N, h_R, &N, work) / matnorm );

        if (flag == 1)
            break;
    }

    /* Memory clean up */
    TESTING_FREE( h_A );
    TESTING_HOSTFREE( h_R );

    TESTING_CUDA_FINALIZE();
}
