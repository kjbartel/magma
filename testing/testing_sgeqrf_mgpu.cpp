/*
    -- MAGMA (version 1.2.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       June 2012

       @generated s Thu Jun 28 12:31:45 2012

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

// Flops formula
#define PRECISION_s
#if defined(PRECISION_z) || defined(PRECISION_c)
#define FLOPS(m, n) ( 6.*FMULS_GEQRF(m, n) + 2.*FADDS_GEQRF(m, n) )
#else
#define FLOPS(m, n) (    FMULS_GEQRF(m, n) +    FADDS_GEQRF(m, n) )
#endif

#define MultiGPUs

/* ////////////////////////////////////////////////////////////////////////////
   -- Testing sgeqrf
*/

int main( int argc, char** argv)
{
    TESTING_CUDA_INIT();
    cudaSetDevice(0);

    magma_timestr_t       start, end;
    float           flops, gpu_perf, cpu_perf;
    float           matnorm, work[1];
    float  c_neg_one = MAGMA_S_NEG_ONE;
    float *h_A, *h_R, *tau, *hwork, tmp[1];
    float *d_lA[4];

    /* Matrix size */
    magma_int_t M = 0, N = 0, n2, n_local[4], lda, ldda, lhwork;
    magma_int_t size[10] = {1024,2048,3072,4032,5184,6016,7040,8064,9088,9984};

    magma_int_t i, k, nk, info, min_mn;
    int max_num_gpus, num_gpus = 1;
    
    magma_int_t ione     = 1;
    magma_int_t ISEED[4] = {0,0,0,1};

    if (argc != 1){
        for(i = 1; i<argc; i++){
            if (strcmp("-N", argv[i])==0)
                N = atoi(argv[++i]);
            else if (strcmp("-M", argv[i])==0)
                M = atoi(argv[++i]);
            else if (strcmp("-NGPU", argv[i])==0)
              num_gpus = atoi(argv[++i]);
        }
        if ( M == 0 ) {
            M = N;
        }
        if ( N == 0 ) {
            N = M;
        }
        if (M>0 && N>0)
          printf("  testing_sgeqrf_gpu -M %d -N %d -NGPU %d\n\n", (int) M, (int) N, (int) num_gpus);
        else
            {
                printf("\nUsage: \n");
                printf("  testing_sgeqrf_gpu -M %d -N %d -NGPU %d\n\n", 
                       1024, 1024, 1);
                exit(1);
            }
    }
    else {
        printf("\nUsage: \n");
        printf("  testing_sgeqrf_gpu -M %d -N %d -NGPU %d\n\n", 1024, 1024, 1);
        M = N = size[9];
    }
    
    ldda   = ((M+31)/32)*32;
    n2     = M * N;
    min_mn = min(M, N);

    magma_int_t nb  = magma_get_sgeqrf_nb(M);

    cudaGetDeviceCount(&max_num_gpus);
    if (num_gpus > max_num_gpus){
      printf("More GPUs requested than available. Have to change it.\n");
      num_gpus = max_num_gpus;
    }
    printf("Number of GPUs to be used = %d\n", (int) num_gpus);

    /* Allocate host memory for the matrix */
    TESTING_MALLOC(    tau, float, min_mn );
    TESTING_MALLOC(    h_A, float, n2     );
    TESTING_HOSTALLOC( h_R, float, n2     );

    for(i=0; i<num_gpus; i++){      
      n_local[i] = ((N/nb)/num_gpus)*nb;
      if (i < (N/nb)%num_gpus)
        n_local[i] += nb;
      else if (i == (N/nb)%num_gpus)
        n_local[i] += N%nb;
      
      #ifdef  MultiGPUs
         cudaSetDevice(i);
      #endif
      TESTING_DEVALLOC(  d_lA[i], float, ldda*n_local[i] );
      printf("device %2d n_local = %4d\n", (int) i, (int) n_local[i]);  
    }
    cudaSetDevice(0);

    lhwork = -1;
    lapackf77_sgeqrf(&M, &N, h_A, &M, tau, tmp, &lhwork, &info);
    lhwork = (magma_int_t)MAGMA_S_REAL( tmp[0] );

    TESTING_MALLOC( hwork, float, lhwork );

    printf("  M     N   CPU GFlop/s   GPU GFlop/s    ||R||_F / ||A||_F\n");
    printf("==========================================================\n");
    for(i=0; i<10; i++){
        if (argc == 1){
            M = N = size[i];
        }
        min_mn= min(M, N);
        lda   = M;
        n2    = lda*N;
        ldda  = ((M+31)/32)*32;
        flops = FLOPS( (float)M, (float)N ) / 1000000;

        /* Initialize the matrix */
        lapackf77_slarnv( &ione, ISEED, &n2, h_A );
        lapackf77_slacpy( MagmaUpperLowerStr, &M, &N, h_A, &lda, h_R, &lda );

        /* =====================================================================
           Performs operation using LAPACK
           =================================================================== */
        start = get_current_time();
        lapackf77_sgeqrf(&M, &N, h_A, &M, tau, hwork, &lhwork, &info);
        end = get_current_time();
        if (info < 0)
            printf("Argument %d of lapack_sgeqrf had an illegal value.\n", (int) -info);

        cpu_perf = flops / GetTimerValue(start, end);

        /* ====================================================================
           Performs operation using MAGMA
           =================================================================== */
        magmablas_ssetmatrix_1D_bcyclic(M, N, h_R, lda, d_lA, ldda, num_gpus, nb);
        
        start = get_current_time();
        magma_sgeqrf2_mgpu( num_gpus, M, N, d_lA, ldda, tau, &info);
        end = get_current_time();

        if (info < 0)
          printf("Argument %d of magma_sgeqrf2 had an illegal value.\n", (int) -info);
        
        gpu_perf = flops / GetTimerValue(start, end);
        
        /* =====================================================================
           Check the result compared to LAPACK
           =================================================================== */
        magmablas_sgetmatrix_1D_bcyclic(M, N, d_lA, ldda, h_R, lda, num_gpus, nb);
        
        matnorm = lapackf77_slange("f", &M, &N, h_A, &M, work);
        blasf77_saxpy(&n2, &c_neg_one, h_A, &ione, h_R, &ione);
        
        printf("%5d %5d  %6.2f         %6.2f        %e\n",
               (int) M, (int) N, cpu_perf, gpu_perf,
               lapackf77_slange("f", &M, &N, h_R, &M, work) / matnorm);
        
        if (argc != 1)
          break;
    }
    
    /* Memory clean up */
    TESTING_FREE( tau );
    TESTING_FREE( h_A );
    TESTING_FREE( hwork );
    TESTING_HOSTFREE( h_R );

    for(i=0; i<num_gpus; i++){
      #ifdef  MultiGPUs
         cudaSetDevice(i);
      #endif
      TESTING_DEVFREE( d_lA[i] );
    }

    /* Shutdown */
    TESTING_CUDA_FINALIZE();
    return EXIT_SUCCESS;
}
