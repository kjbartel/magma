/*
    -- MAGMA (version 1.4.0-beta2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2010

    @author Raffaele Solca
    @author Azzam Haidar

    @generated c Fri Jun 28 19:34:05 2013

*/

// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <cuda_runtime_api.h>
#include <cublas.h>

// includes, project
#include "magma.h"
#include "magma_lapack.h"
#include "testings.h"
#include "magma_cbulge.h"

#define PRECISION_c

#define absv(v1) ((v1)>0? (v1): -(v1))
/* ////////////////////////////////////////////////////////////////////////////
   -- Testing chegvdx
*/
int main( int argc, char** argv)
{

    TESTING_INIT();

    real_Double_t gpu_time;

    magmaFloatComplex *h_A, *h_R, *h_B, *h_S, *h_work;

#if defined(PRECISION_z) || defined(PRECISION_c)
    float *rwork;
    magma_int_t lrwork;
#endif

    /* Matrix size */
    float *w1, *w2, result[2];
    magma_int_t *iwork;
    magma_int_t N, n2, info, lwork, liwork;
    magmaFloatComplex c_zero    = MAGMA_C_ZERO;
    magmaFloatComplex c_one     = MAGMA_C_ONE;
    magmaFloatComplex c_neg_one = MAGMA_C_NEG_ONE;
    magma_int_t ione     = 1;
    magma_int_t ISEED[4] = {0,0,0,1};;

    magma_timestr_t start, end;

    magma_opts opts;
    parse_opts( argc, argv, &opts );

    char jobz = opts.jobz;
    int checkres = opts.check;

    char range = 'A';
    char uplo = opts.uplo;
    magma_int_t itype = opts.itype;

    float f = opts.fraction;

    if (f != 1)
        range='I';

    if ( checkres && jobz == MagmaNoVec ) {
        fprintf( stderr, "checking results requires vectors; setting jobz=V (option -JV)\n" );
        jobz = MagmaVec;
    }

    printf("using: itype = %d, jobz = %c,range = %c, uplo = %c, checkres = %d, fraction = %6.4f\n", itype, jobz, range, uplo, checkres, f);

    printf("  N     M     GPU Time(s) \n");
    printf("==========================\n");
    for( int i = 0; i < opts.ntest; ++i ) {
        for( int iter = 0; iter < opts.niter; ++iter ) {
            N = opts.nsize[i];
            n2     = N*N;
#if defined(PRECISION_z) || defined(PRECISION_c)
            lwork  = magma_cbulge_get_lq2(N) + 2*N + N*N;
            lrwork = 1 + 5*N +2*N*N;
#else
            lwork  = magma_cbulge_get_lq2(N) + 1 + 6*N + 2*N*N;
#endif
            liwork = 3 + 5*N;

            /* Allocate host memory for the matrix */
            TESTING_MALLOC(   h_A, magmaFloatComplex, n2);
            TESTING_MALLOC(   h_B, magmaFloatComplex, n2);
            TESTING_MALLOC(    w1, float         ,  N);
            TESTING_MALLOC(    w2, float         ,  N);
            TESTING_HOSTALLOC(h_R, magmaFloatComplex, n2);
            TESTING_HOSTALLOC(h_S, magmaFloatComplex, n2);
            TESTING_HOSTALLOC(h_work, magmaFloatComplex,  lwork);
#if defined(PRECISION_z) || defined(PRECISION_c)
            TESTING_HOSTALLOC( rwork,          float, lrwork);
#endif
            TESTING_MALLOC(    iwork,     magma_int_t, liwork);

            /* Initialize the matrix */
            lapackf77_clarnv( &ione, ISEED, &n2, h_A );
            lapackf77_clarnv( &ione, ISEED, &n2, h_B );
            /* increase the diagonal */
            for(int i=0; i<N; i++) {
                MAGMA_C_SET2REAL( h_B[i*N+i], ( MAGMA_C_REAL(h_B[i*N+i]) + 1.*N ) );
                MAGMA_C_SET2REAL( h_A[i*N+i], MAGMA_C_REAL(h_A[i*N+i]) );
            }

            magma_int_t m1 = 0;
            float vl = 0;
            float vu = 0;
            magma_int_t il = 0;
            magma_int_t iu = 0;

            if (range == 'I'){
                il = 1;
                iu = f*N;
            }

            // ==================================================================
            // Warmup using MAGMA
            // ==================================================================
            if(opts.warmup){
                lapackf77_clacpy( MagmaUpperLowerStr, &N, &N, h_A, &N, h_R, &N );
                lapackf77_clacpy( MagmaUpperLowerStr, &N, &N, h_B, &N, h_S, &N );

                magma_chegvdx_2stage(itype, jobz, range, uplo,
                                     N, h_R, N, h_S, N, vl, vu, il, iu, &m1, w1,
                                     h_work, lwork,
#if defined(PRECISION_z) || defined(PRECISION_c)
                                     rwork, lrwork,
#endif
                                     iwork, liwork,
                                     &info);
            }
            // ===================================================================
            // Performs operation using MAGMA
            // ===================================================================
            lapackf77_clacpy( MagmaUpperLowerStr, &N, &N, h_A, &N, h_R, &N );
            lapackf77_clacpy( MagmaUpperLowerStr, &N, &N, h_B, &N, h_S, &N );

            start = get_current_time();
            magma_chegvdx_2stage(itype, jobz, range, uplo,
                                 N, h_R, N, h_S, N, vl, vu, il, iu, &m1, w1,
                                 h_work, lwork,
#if defined(PRECISION_z) || defined(PRECISION_c)
                                 rwork, lrwork,
#endif
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
#if defined(PRECISION_d) || defined(PRECISION_s)
                float *rwork = h_work + N*N;
#endif
                float temp1, temp2;

                result[0] = 1.;
                result[0] /= lapackf77_clanhe("1",&uplo, &N, h_A, &N, rwork);
                result[0] /= lapackf77_clange("1",&N , &m1, h_R, &N, rwork);

                if (itype == 1){
                    blasf77_chemm("L", &uplo, &N, &m1, &c_one, h_A, &N, h_R, &N, &c_zero, h_work, &N);
                    for(int i=0; i<m1; ++i)
                        blasf77_csscal(&N, &w1[i], &h_R[i*N], &ione);
                    blasf77_chemm("L", &uplo, &N, &m1, &c_neg_one, h_B, &N, h_R, &N, &c_one, h_work, &N);
                    result[0] *= lapackf77_clange("1", &N, &m1, h_work, &N, rwork)/N;
                }
                else if (itype == 2){
                    blasf77_chemm("L", &uplo, &N, &m1, &c_one, h_B, &N, h_R, &N, &c_zero, h_work, &N);
                    for(int i=0; i<m1; ++i)
                        blasf77_csscal(&N, &w1[i], &h_R[i*N], &ione);
                    blasf77_chemm("L", &uplo, &N, &m1, &c_one, h_A, &N, h_work, &N, &c_neg_one, h_R, &N);
                    result[0] *= lapackf77_clange("1", &N, &m1, h_R, &N, rwork)/N;
                }
                else if (itype == 3){
                    blasf77_chemm("L", &uplo, &N, &m1, &c_one, h_A, &N, h_R, &N, &c_zero, h_work, &N);
                    for(int i=0; i<m1; ++i)
                        blasf77_csscal(&N, &w1[i], &h_R[i*N], &ione);
                    blasf77_chemm("L", &uplo, &N, &m1, &c_one, h_B, &N, h_work, &N, &c_neg_one, h_R, &N);
                    result[0] *= lapackf77_clange("1", &N, &m1, h_R, &N, rwork)/N;
                }

                lapackf77_clacpy( MagmaUpperLowerStr, &N, &N, h_A, &N, h_R, &N );
                lapackf77_clacpy( MagmaUpperLowerStr, &N, &N, h_B, &N, h_S, &N );

                magma_int_t m2 = 0;
                magma_chegvdx(itype, 'N', range, uplo,
                              N, h_R, N, h_S, N, vl, vu, il, iu, &m2, w2,
                              h_work, lwork,
#if defined(PRECISION_z) || defined(PRECISION_c)
                              rwork, lrwork,
#endif
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
                   (int) N, (int) m1, gpu_time);
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

            TESTING_FREE(       h_A);
            TESTING_FREE(       h_B);
            TESTING_FREE(        w1);
            TESTING_FREE(        w2);
#if defined(PRECISION_z) || defined(PRECISION_c)
            TESTING_HOSTFREE( rwork);
#endif
            TESTING_FREE(     iwork);
            TESTING_HOSTFREE(h_work);
            TESTING_HOSTFREE(   h_R);
            TESTING_HOSTFREE(   h_S);
        }
        if ( opts.niter > 1 ) {
            printf( "\n" );
        }
    }

    /* Shutdown */
    TESTING_FINALIZE();

    return 0;
}
