/*
    -- MAGMA (version 1.2.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2010

    @author Raffaele Solca

    @precisions normal z -> c

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
   -- Testing zhegvdx
*/
int main( int argc, char** argv)
{

    TESTING_CUDA_INIT();

    cuDoubleComplex *h_A, *h_R, *h_B, *h_S, *h_work;
    double *rwork, *w1, *w2;
    magma_int_t *iwork;
    double gpu_time, cpu_time;

    magma_timestr_t start, end;

    /* Matrix size */
    magma_int_t N=0, n2;
    magma_int_t size[4] = {1024,2048,4100,6001};

    magma_int_t i, itype, info;
    magma_int_t ione = 1, izero = 0;
    magma_int_t five = 5;

    cuDoubleComplex c_zero    = MAGMA_Z_ZERO;
    cuDoubleComplex c_one     = MAGMA_Z_ONE;
    cuDoubleComplex c_neg_one = MAGMA_Z_NEG_ONE;

    double d_one     =  1.;
    double d_neg_one = -1.;
    double d_ten     = 10.;
    magma_int_t ISEED[4] = {0,0,0,1};

    magma_int_t il,iu,m1,m2;
    double vl,vu;

    double fraction_ev = 0;

    //const char *uplo = MagmaLowerStr;
    char *uplo = (char*)MagmaLowerStr;
    //char *uplo = (char*)MagmaUpperStr;
    char *jobz = (char*)MagmaVectorsStr;
    itype = 1;

    magma_int_t checkres;
    double result[2];

    int flagN = 0;

    if (argc != 1){
        for(i = 1; i<argc; i++){
            if (strcmp("-N", argv[i])==0){
                N = atoi(argv[++i]);
                if (N>0){
                   printf("  testing_zhegvdx -N %d\n\n", N);
                   flagN=1;
                }
                else {
                   printf("\nUsage: \n");
                   printf("  testing_zhegvdx -N %d\n\n", N);
                   exit(1);
                }
            }
            if (strcmp("-itype", argv[i])==0){
                itype = atoi(argv[++i]);
                if (itype>0 && itype <= 3){
                   printf("  testing_zhegvdx -itype %d\n\n", itype);
                }
                else {
                   printf("\nUsage: \n");
                   printf("  testing_zhegvdx -itype %d\n\n", itype);
                   exit(1);
                }
            }
            if (strcmp("-FE", argv[i])==0){
                fraction_ev = atof(argv[++i]);
                if (fraction_ev > 0 && fraction_ev <= 1){
                    printf("  testing_zhegvdx -FE %f\n\n", fraction_ev);
                }
                else {
                    fraction_ev = 0;
                }
            }
            if (strcmp("-L", argv[i])==0){
              uplo = (char*)MagmaLowerStr;
              printf("  testing_zhegvdx -L");
            }
            if (strcmp("-U", argv[i])==0){
              uplo = (char*)MagmaUpperStr;
              printf("  testing_zhegvdx -U");
            }

        }

    } else {
        printf("\nUsage: \n");
        printf("  testing_zhegvdx -L/U -N %d -itype %d\n\n", 1024, 1);
    }

    if(!flagN)
        N = size[3];

    checkres  = getenv("MAGMA_TESTINGS_CHECK") != NULL;

    n2  = N * N;

    /* Allocate host memory for the matrix */
    TESTING_MALLOC(   h_A, cuDoubleComplex, n2);
    TESTING_MALLOC(   h_B, cuDoubleComplex, n2);
    TESTING_MALLOC(    w1, double         ,  N);
    TESTING_MALLOC(    w2, double         ,  N);
    TESTING_HOSTALLOC(h_R, cuDoubleComplex, n2);
    TESTING_HOSTALLOC(h_S, cuDoubleComplex, n2);

    magma_int_t nb = magma_get_zhetrd_nb(N);
    magma_int_t lwork = 2*N*nb + N*N;
    magma_int_t lrwork = 1 + 5*N +2*N*N;
    magma_int_t liwork = 3 + 5*N;

    TESTING_HOSTALLOC(h_work, cuDoubleComplex,  lwork);
    TESTING_MALLOC(    rwork,          double, lrwork);
    TESTING_MALLOC(    iwork,     magma_int_t, liwork);

    printf("\n\n");
    printf("  N     M     GPU Time(s) \n");
    printf("==========================\n");
    for(i=0; i<4; i++){
        if (!flagN){
            N = size[i];
            n2 = N*N;
        }
        if (fraction_ev == 0){
            il = N / 10;
            iu = N / 5+il;
        }
        else {
            il = 1;
            iu = (int)(fraction_ev*N);
            if (iu < 1) iu = 1;
        }

        /* Initialize the matrix */
        lapackf77_zlarnv( &ione, ISEED, &n2, h_A );
        //lapackf77_zlatms( &N, &N, "U", ISEED, "P", w1, &five, &d_ten,
        //                 &d_one, &N, &N, uplo, h_B, &N, h_work, &info);
        //lapackf77_zlaset( "A", &N, &N, &c_zero, &c_one, h_B, &N);
        lapackf77_zlarnv( &ione, ISEED, &n2, h_B );
        /* increase the diagonal */
        {
          magma_int_t i, j;
          for(i=0; i<N; i++) {
            MAGMA_Z_SET2REAL( h_B[i*N+i], ( MAGMA_Z_REAL(h_B[i*N+i]) + 1.*N ) );
          }
        }
        lapackf77_zlacpy( MagmaUpperLowerStr, &N, &N, h_A, &N, h_R, &N );
        lapackf77_zlacpy( MagmaUpperLowerStr, &N, &N, h_B, &N, h_S, &N );

        magma_zhegvdx(itype, jobz[0], 'I', uplo[0],
                     N, h_R, N, h_S, N, vl, vu, il, iu, &m1, w1,
                     h_work, lwork,
                     rwork, lrwork,
                     iwork, liwork,
                     &info);

        lapackf77_zlacpy( MagmaUpperLowerStr, &N, &N, h_A, &N, h_R, &N );
        lapackf77_zlacpy( MagmaUpperLowerStr, &N, &N, h_B, &N, h_S, &N );


        /* ====================================================================
           Performs operation using MAGMA
           =================================================================== */
        start = get_current_time();
        magma_zhegvdx(itype, jobz[0], 'I', uplo[0],
                     N, h_R, N, h_S, N, vl, vu, il, iu, &m1, w1,
                     h_work, lwork,
                     rwork, lrwork,
                     iwork, liwork,
                     &info);
        end = get_current_time();

        gpu_time = GetTimerValue(start,end)/1000.;

        if ( checkres ) {
          /* =====================================================================
             Check the results following the LAPACK's [zc]hegvdx routine.
             A x = lambda B x is solved
             and the following 3 tests computed:
             (1)    | A Z - B Z D | / ( |A||Z| N )  (itype = 1)
                    | A B Z - Z D | / ( |A||Z| N )  (itype = 2)
                    | B A Z - Z D | / ( |A||Z| N )  (itype = 3)
             (2)    | S(with V) - S(w/o V) | / | S |
             =================================================================== */
          double temp1, temp2;
          cuDoubleComplex *tau;

          result[0] = 1.;
          result[0] /= lapackf77_zlanhe("1",uplo, &N, h_A, &N, rwork);
          result[0] /= lapackf77_zlange("1",&N , &m1, h_R, &N, rwork);

          if (itype == 1){
            blasf77_zhemm("L", uplo, &N, &m1, &c_one, h_A, &N, h_R, &N, &c_zero, h_work, &N);
            for(int i=0; i<m1; ++i)
              blasf77_zdscal(&N, &w1[i], &h_R[i*N], &ione);
            blasf77_zhemm("L", uplo, &N, &m1, &c_neg_one, h_B, &N, h_R, &N, &c_one, h_work, &N);
            result[0] *= lapackf77_zlange("1", &N, &m1, h_work, &N, rwork)/N;
          }
          else if (itype == 2){
            blasf77_zhemm("L", uplo, &N, &m1, &c_one, h_B, &N, h_R, &N, &c_zero, h_work, &N);
            for(int i=0; i<m1; ++i)
              blasf77_zdscal(&N, &w1[i], &h_R[i*N], &ione);
            blasf77_zhemm("L", uplo, &N, &m1, &c_one, h_A, &N, h_work, &N, &c_neg_one, h_R, &N);
            result[0] *= lapackf77_zlange("1", &N, &m1, h_R, &N, rwork)/N;
          }
          else if (itype == 3){
            blasf77_zhemm("L", uplo, &N, &m1, &c_one, h_A, &N, h_R, &N, &c_zero, h_work, &N);
            for(int i=0; i<m1; ++i)
              blasf77_zdscal(&N, &w1[i], &h_R[i*N], &ione);
            blasf77_zhemm("L", uplo, &N, &m1, &c_one, h_B, &N, h_work, &N, &c_neg_one, h_R, &N);
            result[0] *= lapackf77_zlange("1", &N, &m1, h_R, &N, rwork)/N;
          }


          lapackf77_zlacpy( MagmaUpperLowerStr, &N, &N, h_A, &N, h_R, &N );
          lapackf77_zlacpy( MagmaUpperLowerStr, &N, &N, h_B, &N, h_S, &N );

          magma_zhegvdx(itype, 'N', 'I', uplo[0],
                       N, h_R, N, h_S, N, vl, vu, il, iu, &m2, w2,
                       h_work, lwork,
                       rwork, lrwork,
                       iwork, liwork,
                       &info);

          temp1 = temp2 = 0;
          for(int j=0; j<m2; j++){
            temp1 = max(temp1, absv(w1[j]));
            temp1 = max(temp1, absv(w2[j]));
            temp2 = max(temp2, absv(w1[j]-w2[j]));
          }
          result[1] = temp2 / temp1;
        }


        /* =====================================================================
           Print execution time
           =================================================================== */
        printf("%5d %5d     %6.2f\n",
               N, m1, gpu_time);
        if ( checkres ){
          printf("Testing the eigenvalues and eigenvectors for correctness:\n");
          if(itype==1)
             printf("(1)    | A Z - B Z D | / (|A| |Z| N) = %e\n", result[0]);
          else if(itype==2)
             printf("(1)    | A B Z - Z D | / (|A| |Z| N) = %e\n", result[0]);
          else if(itype==3)
             printf("(1)    | B A Z - Z D | / (|A| |Z| N) = %e\n", result[0]);

          printf("(2)    | D(w/ Z)-D(w/o Z)|/ |D| = %e\n\n", result[1]);
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
