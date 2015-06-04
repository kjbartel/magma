/*
 *  -- MAGMA (version 1.2.1) --
 *     Univ. of Tennessee, Knoxville
 *     Univ. of California, Berkeley
 *     Univ. of Colorado, Denver
 *     June 2012
 *
 * @generated c Thu Jun 28 12:31:43 2012
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
#define PRECISION_c
#if defined(PRECISION_z) || defined(PRECISION_c)
#define FLOPS(m, n) ( 6. * FMULS_GETRF(m, n) + 2. * FADDS_GETRF(m, n) )
#else
#define FLOPS(m, n) (      FMULS_GETRF(m, n) +      FADDS_GETRF(m, n) )
#endif

/* ========= definition of multiple GPU code ========= */
extern "C" magma_int_t
magma_cgetrf_mgpu(magma_int_t num_gpus, 
                     magma_int_t m, magma_int_t n,
                 cuFloatComplex **d_lA, magma_int_t ldda,
                         magma_int_t *ipiv, magma_int_t *info);
/* =================================================== */



float get_LU_error(magma_int_t M, magma_int_t N, 
                    cuFloatComplex *A,  magma_int_t lda, 
                    cuFloatComplex *LU, magma_int_t *IPIV)
{
    magma_int_t min_mn = min(M,N);
    magma_int_t ione   = 1;
    magma_int_t i, j;
    cuFloatComplex alpha = MAGMA_C_ONE;
    cuFloatComplex beta  = MAGMA_C_ZERO;
    cuFloatComplex *L, *U;
    float work[1], matnorm, residual;
                       
    TESTING_MALLOC( L, cuFloatComplex, M*min_mn);
    TESTING_MALLOC( U, cuFloatComplex, min_mn*N);
    memset( L, 0, M*min_mn*sizeof(cuFloatComplex) );
    memset( U, 0, min_mn*N*sizeof(cuFloatComplex) );

    lapackf77_claswp( &N, A, &lda, &ione, &min_mn, IPIV, &ione);
    lapackf77_clacpy( MagmaLowerStr, &M, &min_mn, LU, &lda, L, &M      );
    lapackf77_clacpy( MagmaUpperStr, &min_mn, &N, LU, &lda, U, &min_mn );

    for(j=0; j<min_mn; j++)
        L[j+j*M] = MAGMA_C_MAKE( 1., 0. );
    
    matnorm = lapackf77_clange("f", &M, &N, A, &lda, work);

    blasf77_cgemm("N", "N", &M, &N, &min_mn,
                  &alpha, L, &M, U, &min_mn, &beta, LU, &lda);

    for( j = 0; j < N; j++ ) {
        for( i = 0; i < M; i++ ) {
            LU[i+j*lda] = MAGMA_C_SUB( LU[i+j*lda], A[i+j*lda] );
        }
    }
    residual = lapackf77_clange("f", &M, &N, LU, &lda, work);

    TESTING_FREE(L);
    TESTING_FREE(U);

    return residual / (matnorm * N);
}

/* ////////////////////////////////////////////////////////////////////////////
   -- Testing cgetrf
*/
int main( int argc, char** argv)
{
    TESTING_CUDA_INIT();

    magma_timestr_t  start, end;
    float           flops, gpu_perf, cpu_perf, error;
    cuFloatComplex *h_A, *h_R;
    cuFloatComplex *d_lA[4];
    magma_int_t     *ipiv;

    /* Matrix size */
    magma_int_t M = 0, N = 0, flag = 0, n2, lda, ldda, num_gpus, num_gpus0 = 1, n_local;
    magma_int_t size[10] = {960,1920,3072,4032,4992,5952,7104,8064,9024,9984};
    magma_int_t n_size = 10;

    magma_int_t i, k, info, min_mn, nb0, nb, nk, maxn, ret, ldn_local;
    magma_int_t ione     = 1;
    magma_int_t ISEED[4] = {0,0,0,1};

    if (argc != 1){
        for(i = 1; i<argc; i++){
            if (strcmp("-N", argv[i])==0) {
                    N = atoi(argv[++i]);
                    flag = 1;
            } else if (strcmp("-M", argv[i])==0) {
                    M = atoi(argv[++i]);
                    flag = 1;
            } else if (strcmp("-NGPU", argv[i])==0)
                num_gpus0 = atoi(argv[++i]);
        }
    }
    if( flag != 0 ) {
        if (M>0 && N>0 && num_gpus0>0)
            printf("  testing_cgetrf_mgpu -M %d -N %d -NGPU %d\n\n", (int) M, (int) N, (int) num_gpus0);
        else {
            printf("\nError: \n");
            printf("  (m, n, num_gpus)=(%d, %d, %d) must be positive.\n\n", (int) M, (int) N, (int) num_gpus0);
            exit(1);
        }
    } else {
        M = N = size[n_size-1];
        printf("\nDefault: \n");
        printf("  testing_cgetrf_mgpu -M %d -N %d -NGPU %d\n\n", (int) M, (int) N, (int) num_gpus0);
    }

    ldda   = ((M+31)/32)*32;
    maxn   = ((N+31)/32)*32;
    n2     = M * N;
    min_mn = min(M, N);
    nb     = magma_get_cgetrf_nb(M);
        num_gpus = num_gpus0;

    /* Allocate host memory for the matrix */
    TESTING_MALLOC(ipiv, magma_int_t, min_mn);
    TESTING_MALLOC(    h_A, cuFloatComplex, n2     );
    TESTING_HOSTALLOC( h_R, cuFloatComplex, n2     );
        /* allocate device memory, assuming fixed nb and num_gpus */
    for(i=0; i<num_gpus; i++){
          n_local = ((N/nb)/num_gpus)*nb;
          if (i < (N/nb)%num_gpus)
            n_local += nb;
          else if (i == (N/nb)%num_gpus)
            n_local += N%nb;
          ldn_local = ((n_local+31)/32)*32;
      cudaSetDevice(i);
      //TESTING_DEVALLOC( d_lA[i], cuFloatComplex, ldda*n_local );
      TESTING_DEVALLOC( d_lA[i], cuFloatComplex, ldda*ldn_local );
    }
    cudaSetDevice(0);
    nb0 = nb;

    printf("  M     N   CPU GFlop/s    GPU GFlop/s   ||PA-LU||/(||A||*N)\n");
    printf("============================================================\n");
    for(i=0; i<n_size; i++){
        if (flag == 0){
              M = N = size[i];
        }
            min_mn= min(M, N);
            lda   = M;
            n2    = lda*N;
            ldda  = ((M+31)/32)*32;
            flops = FLOPS( (float)M, (float)N ) / 1000000;

        /* Initialize the matrix */
        lapackf77_clarnv( &ione, ISEED, &n2, h_A );
        lapackf77_clacpy( MagmaUpperLowerStr, &M, &N, h_A, &lda, h_R, &lda );

       /* =====================================================================
           Performs operation using LAPACK
           =================================================================== */
        start = get_current_time();
        lapackf77_cgetrf(&M, &N, h_A, &lda, ipiv, &info);
        end = get_current_time();
        if (info < 0) {
                        printf("Argument %d of cgetrf had an illegal value.\n", (int) -info);
                        break;
                } else if (info != 0 ) {
                        printf("cgetrf returned info=%d.\n", (int) info);
                        break;
                }
        cpu_perf = flops / GetTimerValue(start, end);
                lapackf77_clacpy( MagmaUpperLowerStr, &M, &N, h_R, &lda, h_A, &lda );

        /* ====================================================================
           Performs operation using MAGMA
           =================================================================== */
                /* == distributing the matrix == */
                //cudaSetDevice(0);
        //cublasSetMatrix( M, N, sizeof(cuFloatComplex), h_R, lda, d_lA[0], ldda);
        nb = magma_get_cgetrf_nb(M);
#ifdef TESTING_CGETRF_MGPU_CHECK
                if( nb != nb0 ) {
                  printf( " different nb used for memory-allocation (%d vs. %d)\n", (int) nb, (int) nb0 );
                }
#endif
        if( num_gpus0 > N/nb ) {
                  num_gpus = N/nb;
                  if( N%nb != 0 ) num_gpus ++;
                  printf( " * too many GPUs for the matrix size, using %d GPUs\n", (int) num_gpus );
                } else {
                  num_gpus = num_gpus0;
                }

                for(int j=0; j<N; j+=nb){
                  k = (j/nb)%num_gpus;
                  cudaSetDevice(k);
                  nk = min(nb, N-j);
                  magma_csetmatrix( M, nk,
                                    h_R+j*lda,                       lda,
                                    d_lA[k]+j/(nb*num_gpus)*nb*ldda, ldda );
                }
                cudaSetDevice(0);

                /* == calling MAGMA with multiple GPUs == */
        start = get_current_time();
        magma_cgetrf_mgpu( num_gpus, M, N, d_lA, ldda, ipiv, &info);
        end = get_current_time();
        gpu_perf = flops / GetTimerValue(start, end);
        if (info < 0) {
            printf("Argument %d of magma_cgetrf_mgpu had an illegal value.\n", (int) -info);
                        break;
                } else if (info != 0 ) {
            printf("magma_cgetrf_mgpu returned info=%d.\n", (int) info);
                        break;
                }
                /* == download the matrix from GPUs == */
                //cudaSetDevice(0);
                //cublasGetMatrix( M, N, sizeof(cuFloatComplex), d_lA[0], ldda, h_R, M);
        for(int j=0; j<N; j+=nb){
                  k = (j/nb)%num_gpus;
                  cudaSetDevice(k);
                  nk = min(nb, N-j);
                  magma_cgetmatrix( M, nk,
                                    d_lA[k]+j/(nb*num_gpus)*nb*ldda, ldda,
                                    h_R+j*lda,                       lda );
                }
                cudaSetDevice(0);

        /* =====================================================================
           Check the factorization
           =================================================================== */
        error = get_LU_error(M, N, h_A, lda, h_R, ipiv);
        
        printf("%5d %5d  %6.2f         %6.2f         %e\n",
               (int) M, (int) N, cpu_perf, gpu_perf, error);

        if (flag != 0)
            break;
    }

    /* Memory clean up */
    TESTING_FREE( ipiv );
    TESTING_FREE( h_A );
    TESTING_HOSTFREE( h_R );
    for(i=0; i<num_gpus0; i++){
                TESTING_DEVFREE( d_lA[i] );
        }

    /* Shutdown */
    TESTING_CUDA_FINALIZE();
}
