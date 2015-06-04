/*
 *  -- MAGMA (version 1.3.0) --
 *     Univ. of Tennessee, Knoxville
 *     Univ. of California, Berkeley
 *     Univ. of Colorado, Denver
 *     November 2012
 *
 * @precisions normal z -> c d s
 *
 **/
/* includes, system */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <cuda_runtime_api.h>
#include <cublas.h>

/* includes, project */
#include "flops.h"
#include "magma.h"
#include "magma_lapack.h"
#include "testings.h"

#define PRECISION_z
/* Flops formula */
#if defined(PRECISION_z) || defined(PRECISION_c)
#define FLOPS(n) ( 6. * FMULS_POTRF(n) + 2. * FADDS_POTRF(n) )
#else
#define FLOPS(n) (      FMULS_POTRF(n) +      FADDS_POTRF(n) )
#endif

/* definitions for multi-GPU code */
extern "C" magma_int_t
magma_zpotrf_mgpu(int num_gpus, char uplo, magma_int_t n,
                      cuDoubleComplex **d_lA, magma_int_t ldda, magma_int_t *info);
extern "C" magma_int_t
magma_zpotrf1_mgpu(int num_gpus, char uplo, magma_int_t n,
                      cuDoubleComplex **d_lA, magma_int_t ldda, magma_int_t *info);


/* ////////////////////////////////////////////////////////////////////////////
   -- Testing zpotrf_mgpu
*/
int main( int argc, char** argv) 
{
    TESTING_CUDA_INIT();
    magma_setdevice(0);

    magma_timestr_t  start, end;
    double      flops, gpu_perf, cpu_perf;
    cuDoubleComplex *h_A, *h_R;
    cuDoubleComplex *d_lA[4];
    magma_int_t N = 0, n2, mb, nb, nk, lda, ldda, n_local, ldn_local;
    //magma_int_t size[10] = {1000,2000,3000,4000,5000,6000,7000,8000,9000,10000};
    magma_int_t size[10] = {1024,2048,3072,4032,5184,6016,7040,8064,9088,10112};
    magma_int_t n_sizes = 10, flag = 0;
    
    magma_int_t i, j, k, info, num_gpus0 = 1, num_gpus;
    const char *uplo     = MagmaLowerStr;
    cuDoubleComplex c_neg_one = MAGMA_Z_NEG_ONE;
    magma_int_t ione     = 1;
    magma_int_t ISEED[4] = {0,0,0,1};
    double      work[1], matnorm;
    
    N = size[n_sizes-1];
    if (argc != 1){
        for(i = 1; i<argc; i++){        
            if (strcmp("-N", argv[i])==0) {
                flag = 1;
                N = atoi(argv[++i]);
                size[0] = size[n_sizes-1] = N;
            }
            if (strcmp("-NGPU", argv[i])==0)
                num_gpus0 = atoi(argv[++i]);
            if (strcmp("-UPLO",argv[i])==0) {
                if (strcmp("L",argv[++i])==0) uplo = MagmaLowerStr;
                else                          uplo = MagmaUpperStr;
            }
        }
        if (strcmp(uplo,MagmaLowerStr)==0)
        printf("\n  testing_zpotrf_mgpu -N %d -NGPU %d -UPLO L\n\n", (int) N, (int) num_gpus0 );
        else
        printf("\n  testing_zpotrf_mgpu -N %d -NGPU %d -UPLO U\n\n", (int) N, (int) num_gpus0 );
    } else {
        printf("\nDefault: \n");
        printf("  testing_zpotrf_mgpu -N %d:%d -NGPU %d -UPLO L\n\n", (int) size[0], (int) size[n_sizes-1], (int) num_gpus0 );
    }
    if( N <= 0 || num_gpus0 <= 0 )  {
        printf( " invalid input N=%d NGPU=%d\n", (int) N, (int) num_gpus0 );
        exit(1);
    }

    /* looking for max. ldda */
    ldda = 0;
    n2   = 0;
    for(i=0; i<n_sizes; i++){
        N     = size[i];
        nb = magma_get_zpotrf_nb(N);
        mb = nb;
        if( num_gpus0 > N/nb ) {
            num_gpus = N/nb;
            if( N%nb != 0 ) num_gpus ++;
        } else {
            num_gpus = num_gpus0;
        }
        n_local = nb*(1+N/(nb*num_gpus)) * mb*((N+mb-1)/mb);
        if( n_local > ldda ) ldda = n_local;
        if( n2 < N*N ) n2 = N*N;
        if (flag != 0) break;
    }

    /* Allocate host memory for the matrix */
    TESTING_HOSTALLOC( h_A, cuDoubleComplex, n2);
    TESTING_HOSTALLOC( h_R, cuDoubleComplex, n2);
    /* allocate local matrix on GPU */
    for(i=0; i<num_gpus0; i++){
        magma_setdevice(i);
        TESTING_DEVALLOC( d_lA[i], cuDoubleComplex, ldda );
    }
    magma_setdevice(0);

    printf("  N    CPU GFlop/s    GPU GFlop/s    ||R||_F / ||A||_F\n");
    printf("========================================================\n");
    for(i=0; i<n_sizes; i++){
        N     = size[i];
        lda   = N; 
        n2    = lda*N;
        flops = FLOPS( (double)N ) / 1000000;
        
        /* Initialize the matrix */
        lapackf77_zlarnv( &ione, ISEED, &n2, h_A );
        /* Symmetrize and increase the diagonal */
        {
            magma_int_t i, j;
            for(i=0; i<N; i++) {
                MAGMA_Z_SET2REAL( h_A[i*lda+i], ( MAGMA_Z_REAL(h_A[i*lda+i]) + 1.*N ) );
                for(j=0; j<i; j++)
                   h_A[i*lda+j] = cuConj(h_A[j*lda+i]);
            }
        }
        lapackf77_zlacpy( MagmaUpperLowerStr, &N, &N, h_A, &lda, h_R, &lda );

        /* ====================================================================
           Performs operation using MAGMA 
           =================================================================== */

        nb = magma_get_zpotrf_nb(N);
        if( num_gpus0 > N/nb ) {
            num_gpus = N/nb;
            if( N%nb != 0 ) num_gpus ++;
            printf( " * too many GPUs for the matrix size, using %d GPUs\n", (int) num_gpus );
        } else {
            num_gpus = num_gpus0;
        }

        /* distribute matrix to gpus */
        if( lapackf77_lsame(uplo, "U") ) {
            /* going through each block-column */
            ldda  = ((N+mb-1)/mb)*mb;
            for(j=0; j<N; j+=nb){
              k = (j/nb)%num_gpus;
              magma_setdevice(k);
              nk = min(nb, N-j);
              magma_zsetmatrix( N, nk,
                                h_A+j*lda,                       lda,
                                d_lA[k]+j/(nb*num_gpus)*nb*ldda, ldda );
            }
        } else {
            /* going through each block-row */
            ldda = (1+N/(nb*num_gpus))*nb;
            for(j=0; j<N; j+=nb){
              k = (j/nb)%num_gpus;
              magma_setdevice(k);
              nk = min(nb, N-j);
              magma_zsetmatrix( nk, N,
                                h_A+j,                      lda,
                                d_lA[k]+j/(nb*num_gpus)*nb, ldda );
            }
        }
        magma_setdevice(0);

        /* call magma_zpotrf_mgpu */
        start = get_current_time();
        magma_zpotrf_mgpu(num_gpus, uplo[0], N, d_lA, ldda, &info);
        end = get_current_time();
        if (info < 0) {
            printf("Argument %d of magma_zpotrf_mgpu had an illegal value.\n", (int) -info);
            break;
        } else if (info != 0) {
            printf("magma_zpotrf_mgpu returned info=%d\n", (int) info );
            break;
        }
        gpu_perf = flops / GetTimerValue(start, end);
        
        /* gather matrix from gpus */
        if( lapackf77_lsame(uplo, "U") ) {
            for(j=0; j<N; j+=nb){
                k = (j/nb)%num_gpus;
                magma_setdevice(k);
                nk = min(nb, N-j);
                magma_zgetmatrix( N, nk,
                                  d_lA[k]+j/(nb*num_gpus)*nb*ldda, ldda,
                                  h_R+j*lda,                       lda );
            }
        } else {
            for(j=0; j<N; j+=nb){
              k = (j/nb)%num_gpus;
              magma_setdevice(k);
              nk = min(nb, N-j);
              magma_zgetmatrix( nk, N,
                                d_lA[k]+j/(nb*num_gpus)*nb, ldda,
                                h_R+j,                      lda );
            }
        }
        magma_setdevice(0);

        /* =====================================================================
           Performs operation using LAPACK 
           =================================================================== */
        start = get_current_time();
        lapackf77_zpotrf(uplo, &N, h_A, &lda, &info);
        end = get_current_time();
        if (info < 0) {
              printf("Argument %d of zpotrf had an illegal value.\n", (int) -info);
              break;
        } else if (info != 0) {
              printf("lapackf77_zpotrf returned info=%d\n", (int) info );
              break;
        }
        cpu_perf = flops / GetTimerValue(start, end);
      
        /* =====================================================================
           Check the result compared to LAPACK
           =================================================================== */
        matnorm = lapackf77_zlange("f", &N, &N, h_A, &lda, work);
        blasf77_zaxpy(&n2, &c_neg_one, h_A, &ione, h_R, &ione);
        printf("%5d    %6.2f         %6.2f        %e\n", 
               (int) size[i], cpu_perf, gpu_perf,
               lapackf77_zlange("f", &N, &N, h_R, &lda, work) / matnorm);
        
        if (flag != 0) break;
    }

    /* Memory clean up */
    TESTING_HOSTFREE( h_A );
    TESTING_HOSTFREE( h_R );
    for(i=0; i<num_gpus; i++){
      magma_setdevice(i);
      TESTING_DEVFREE( d_lA[i] );
    }

    /* Shutdown */
    TESTING_CUDA_FINALIZE();
}
