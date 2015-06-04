/*
    -- MAGMA (version 1.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2011

    @author Raffaele Solca

    @generated c Sun Nov 13 20:48:55 2011

*/

// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cublas.h>

// includes, project
#include "magma.h"
#include "magma_lapack.h"
#include "testings.h"

#define absv(v1) ((v1)>0? (v1): -(v1))

/* ////////////////////////////////////////////////////////////////////////////
   -- Testing chegvd
*/
int main( int argc, char** argv) 
{
    TESTING_CUDA_INIT();

    cuFloatComplex *h_A, *h_R, *h_B, *h_S, *h_work;
    float *rwork, *w1, *w2;
    magma_int_t *iwork;
    float gpu_time, cpu_time;

    magma_timestr_t start, end;

    /* Matrix size */
    magma_int_t N=0, n2;
    magma_int_t size[4] = {1024,2048,4100,6001};

    magma_int_t i, itype, info;
    magma_int_t ione = 1, izero = 0;
    magma_int_t five = 5;

    cuFloatComplex czero = MAGMA_C_ZERO;
    cuFloatComplex zone  = MAGMA_C_ONE;
    cuFloatComplex mzone  = MAGMA_C_NEG_ONE;

    float done = 1.;
    float mdone = -1.;
    float dten = 10.;
    magma_int_t ISEED[4] = {0,0,0,1};

    //const char *uplo = MagmaLowerStr;
    char *uplo = (char*)MagmaLowerStr;
    //char *uplo = (char*)MagmaUpperStr;
    char *jobz = (char*)MagmaVectorsStr;
    itype = 1;

    magma_int_t checkres;
    float result[4];

    int flagN = 0;

    if (argc != 1){
        for(i = 1; i<argc; i++){
            if (strcmp("-N", argv[i])==0){
                N = atoi(argv[++i]);
                if (N>0){
                   printf("  testing_chegvd -N %d\n\n", N);
                   flagN=1;
                }
                else {
                   printf("\nUsage: \n");
                   printf("  testing_chegvd -N %d\n\n", N);
                   exit(1);
                }
            }
            if (strcmp("-itype", argv[i])==0){
                itype = atoi(argv[++i]);
                if (itype>0 && itype <= 3){
                   printf("  testing_chegvd -itype %d\n\n", itype);
                }
                else {
                   printf("\nUsage: \n");
                   printf("  testing_chegvd -itype %d\n\n", itype);
                   exit(1);
                }
            }
            if (strcmp("-L", argv[i])==0){
              uplo = (char*)MagmaLowerStr;
              printf("  testing_chegvd -L");
            }
            if (strcmp("-U", argv[i])==0){
              uplo = (char*)MagmaUpperStr;
              printf("  testing_chegvd -U");              
            }
          
        }
      
    } else {
        printf("\nUsage: \n");
        printf("  testing_chegvd -L/U -N %d -itype %d\n\n", 1024, 1);
    }

    if(!flagN)
        N = size[3];

    checkres  = getenv("MAGMA_TESTINGS_CHECK") != NULL;
    n2  = N * N;

    /* Allocate host memory for the matrix */
    TESTING_MALLOC(   h_A, cuFloatComplex, n2);
    TESTING_MALLOC(   h_B, cuFloatComplex, n2);
    TESTING_MALLOC(    w1, float         ,  N);
    TESTING_MALLOC(    w2, float         ,  N);
    TESTING_HOSTALLOC(h_R, cuFloatComplex, n2);
    TESTING_HOSTALLOC(h_S, cuFloatComplex, n2);

    magma_int_t nb = magma_get_chetrd_nb(N);
    magma_int_t lwork = 2*N*nb + N*N;
    magma_int_t lrwork = 1 + 5*N +2*N*N;
    magma_int_t liwork = 3 + 5*N;

    TESTING_HOSTALLOC(h_work, cuFloatComplex,  lwork);
    TESTING_MALLOC(    rwork,          float, lrwork);
    TESTING_MALLOC(    iwork,     magma_int_t, liwork);
    
    printf("\n\n");
    printf("  N     CPU Time(s)    GPU Time(s) \n");
    printf("===================================\n");
    for(i=0; i<4; i++){
        if (!flagN){
            N = size[i];
            n2 = N*N;
        }

        /* Initialize the matrix */
        lapackf77_clarnv( &ione, ISEED, &n2, h_A );
        //lapackf77_clatms( &N, &N, "U", ISEED, "P", w1, &five, &dten,
        //                 &done, &N, &N, uplo, h_B, &N, h_work, &info);
        //lapackf77_claset( "A", &N, &N, &czero, &zone, h_B, &N);
        lapackf77_clarnv( &ione, ISEED, &n2, h_B );
        /* increase the diagonal */
        {
          magma_int_t i, j;
          for(i=0; i<N; i++) {
            MAGMA_C_SET2REAL( h_B[i*N+i], ( MAGMA_C_GET_X(h_B[i*N+i]) + 1.*N ) );
            MAGMA_C_SET2REAL( h_A[i*N+i], MAGMA_C_GET_X(h_A[i*N+i]) );
          }
        }
        lapackf77_clacpy( MagmaUpperLowerStr, &N, &N, h_A, &N, h_R, &N );
        lapackf77_clacpy( MagmaUpperLowerStr, &N, &N, h_B, &N, h_S, &N );

        magma_chegvd(itype, jobz[0], uplo[0],
                     N, h_R, N, h_S, N, w1,
                     h_work, lwork, 
                     rwork, lrwork, 
                     iwork, liwork, 
                     &info);
        
        lapackf77_clacpy( MagmaUpperLowerStr, &N, &N, h_A, &N, h_R, &N );
        lapackf77_clacpy( MagmaUpperLowerStr, &N, &N, h_B, &N, h_S, &N );


        /* ====================================================================
           Performs operation using MAGMA
           =================================================================== */
        start = get_current_time();
        magma_chegvd(itype, jobz[0], uplo[0],
                     N, h_R, N, h_S, N, w1,
                     h_work, lwork,
                     rwork, lrwork,
                     iwork, liwork,
                     &info);
        end = get_current_time();

        gpu_time = GetTimerValue(start,end)/1000.;

        if ( checkres ) {
          /* =====================================================================
             Check the results following the LAPACK's [zc]hegvd routine.
             A x = lambda B x is solved
             and the following 3 tests computed:
             (1)    | A Z - B Z D | / ( |A||Z| N )  (itype = 1)
                    | A B Z - Z D | / ( |A||Z| N )  (itype = 2)
                    | B A Z - Z D | / ( |A||Z| N )  (itype = 3)
             (2)    | I - V V' B | / ( N )           (itype = 1,2)
                    | B - V V' | / ( |B| N )         (itype = 3)
             (3)    | S(with V) - S(w/o V) | / | S |
             =================================================================== */
          float temp1, temp2;
          cuFloatComplex *tau;

          if (itype == 1 || itype == 2){
            lapackf77_claset( "A", &N, &N, &czero, &zone, h_S, &N);
            blasf77_cgemm("N", "C", &N, &N, &N, &zone, h_R, &N, h_R, &N, &czero, h_work, &N);
            blasf77_chemm("R", uplo, &N, &N, &mzone, h_B, &N, h_work, &N, &zone, h_S, &N);
            result[1]= lapackf77_clange("1", &N, &N, h_S, &N, rwork) / N;
          }
          else if (itype == 3){
            lapackf77_clacpy( MagmaUpperLowerStr, &N, &N, h_B, &N, h_S, &N);
            blasf77_cherk(uplo, "N", &N, &N, &mdone, h_R, &N, &done, h_S, &N); 
            result[1]= lapackf77_clanhe("1",uplo, &N, h_S, &N, rwork) / N / lapackf77_clanhe("1",uplo, &N, h_B, &N, rwork);
          }

          result[0] = 1.;
          result[0] /= lapackf77_clanhe("1",uplo, &N, h_A, &N, rwork);
          result[0] /= lapackf77_clange("1",&N , &N, h_R, &N, rwork);

          if (itype == 1){
            blasf77_chemm("L", uplo, &N, &N, &zone, h_A, &N, h_R, &N, &czero, h_work, &N);
            for(int i=0; i<N; ++i)
              blasf77_csscal(&N, &w1[i], &h_R[i*N], &ione);
            blasf77_chemm("L", uplo, &N, &N, &mzone, h_B, &N, h_R, &N, &zone, h_work, &N);
            result[0] *= lapackf77_clange("1", &N, &N, h_work, &N, rwork)/N;
          }
          else if (itype == 2){
            blasf77_chemm("L", uplo, &N, &N, &zone, h_B, &N, h_R, &N, &czero, h_work, &N);
            for(int i=0; i<N; ++i)
              blasf77_csscal(&N, &w1[i], &h_R[i*N], &ione);
            blasf77_chemm("L", uplo, &N, &N, &zone, h_A, &N, h_work, &N, &mzone, h_R, &N);
            result[0] *= lapackf77_clange("1", &N, &N, h_R, &N, rwork)/N;
          }
          else if (itype == 3){
            blasf77_chemm("L", uplo, &N, &N, &zone, h_A, &N, h_R, &N, &czero, h_work, &N);
            for(int i=0; i<N; ++i)
              blasf77_csscal(&N, &w1[i], &h_R[i*N], &ione);
            blasf77_chemm("L", uplo, &N, &N, &zone, h_B, &N, h_work, &N, &mzone, h_R, &N);
            result[0] *= lapackf77_clange("1", &N, &N, h_R, &N, rwork)/N;
          }

/*          lapackf77_chet21(&ione, uplo, &N, &izero,
                           h_A, &N,
                           w1, w1,
                           h_R, &N,
                           h_R, &N,
                           tau, h_work, rwork, &result[0]);
*/          
          lapackf77_clacpy( MagmaUpperLowerStr, &N, &N, h_A, &N, h_R, &N );
          lapackf77_clacpy( MagmaUpperLowerStr, &N, &N, h_B, &N, h_S, &N );
 
          magma_chegvd(itype, 'N', uplo[0],
                       N, h_R, N, h_S, N, w2,
                       h_work, lwork,
                       rwork, lrwork,
                       iwork, liwork,
                       &info);

          temp1 = temp2 = 0;
          for(int j=0; j<N; j++){
            temp1 = max(temp1, absv(w1[j]));
            temp1 = max(temp1, absv(w2[j]));
            temp2 = max(temp2, absv(w1[j]-w2[j]));
          }
          result[2] = temp2 / temp1;
        }

        /* =====================================================================
           Performs operation using LAPACK
           =================================================================== */
        start = get_current_time();
        lapackf77_chegvd(&itype, jobz, uplo,
                         &N, h_A, &N, h_B, &N, w2,
                         h_work, &lwork,
                         rwork, &lrwork,
                         iwork, &liwork,
                         &info);
        end = get_current_time();
        if (info < 0)
          printf("Argument %d of chegvd had an illegal value.\n", -info);

        cpu_time = GetTimerValue(start,end)/1000.;


        /* =====================================================================
           Print execution time
           =================================================================== */
        printf("%5d     %6.2f         %6.2f\n",
               N, cpu_time, gpu_time);
        if ( checkres ){
          printf("Testing the eigenvalues and eigenvectors for correctness:\n");
          if(itype==1)
             printf("(1)    | A Z - B Z D | / (|A| |Z| N) = %e\n", result[0]);
          else if(itype==2)
             printf("(1)    | A B Z - Z D | / (|A| |Z| N) = %e\n", result[0]);
          else if(itype==3)
             printf("(1)    | B A Z - Z D | / (|A| |Z| N) = %e\n", result[0]);
          if(itype==1 || itype ==2)
             printf("(2)    | I -   Z Z' B | /  N      = %e\n", result[1]);
          else
             printf("(2)    | B -  Z Z' | / (|B| N)      = %e\n", result[1]);
          printf("(3)    | D(w/ Z)-D(w/o Z)|/ |D| = %e\n\n", result[2]);
        }

        if (flagN)
            break;
    }
 
    /* Memory clean up */
    TESTING_FREE(       h_A);
    TESTING_FREE(       h_B);
    TESTING_FREE(        w1);
    TESTING_FREE(        w2);
    TESTING_FREE(     rwork);
    TESTING_FREE(     iwork);
    TESTING_HOSTFREE(h_work);
    TESTING_HOSTFREE(   h_R);
    TESTING_HOSTFREE(   h_S);

    /* Shutdown */
    TESTING_CUDA_FINALIZE();
}
