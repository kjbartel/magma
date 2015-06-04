/*
    -- MAGMA (version 1.3.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2012

       @generated s Wed Nov 14 22:54:26 2012

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
#define PRECISION_s

/* ////////////////////////////////////////////////////////////////////////////
   -- Testing sgesvd
*/
int main( int argc, char** argv)
{
    TESTING_CUDA_INIT();

    real_Double_t   gpu_time, cpu_time;
    float *h_A, *h_R, *U, *VT, *h_work;
    float *S1, *S2;
#if defined(PRECISION_z) || defined(PRECISION_c)
    float *rwork;
#endif

    /* Matrix size */
    magma_int_t M=0, N=0, n2, min_mn;
    const int MAXTESTS = 10;
    magma_int_t msize[MAXTESTS] = { 1024, 2048, 3072, 4032, 5184, 6016, 7040, 8064, 9088, 10112 };
    magma_int_t nsize[MAXTESTS] = { 1024, 2048, 3072, 4032, 5184, 6016, 7040, 8064, 9088, 10112 };

    magma_int_t info;
    magma_int_t ione     = 1;
    magma_int_t ISEED[4] = {0,0,0,1};
    
    const char* jobu = "S";
    const char* jobv = "S";

    int checkres = getenv("MAGMA_TESTINGS_CHECK") != NULL;
    int lapack   = getenv("MAGMA_RUN_LAPACK")     != NULL;
    int test_all = false;
    int workspace = 1;
    
    // process command line arguments
    printf( "\nUsage: %s -N <m,n> -U[ASON] -V[ASON] -all -c -l -w[123]\n"
            "  -N can be repeated up to %d times. If only m is given, then m=n.\n"
            "  -c or setting $MAGMA_TESTINGS_CHECK checks result.\n"
            "  -l or setting $MAGMA_RUN_LAPACK runs LAPACK and checks singular values.\n"
            "  -U* and -V* set jobu and jobv.\n"
            "  -all tests all 15 combinations of jobu and jobv.\n"
            "  -w* sets workspace size, from default min (1) to max (3).\n\n",
            argv[0], MAXTESTS );
    
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
        else if ( strcmp("-UA", argv[i]) == 0 )
            jobu = "A";
        else if ( strcmp("-US", argv[i]) == 0 )
            jobu = "S";
        else if ( strcmp("-UO", argv[i]) == 0 )
            jobu = "O";
        else if ( strcmp("-UN", argv[i]) == 0 )
            jobu = "N";
        
        else if ( strcmp("-VA", argv[i]) == 0 )
            jobv = "A";
        else if ( strcmp("-VS", argv[i]) == 0 )
            jobv = "S";
        else if ( strcmp("-VO", argv[i]) == 0 )
            jobv = "O";
        else if ( strcmp("-VN", argv[i]) == 0 )
            jobv = "N";
        
        else if ( strcmp("-all", argv[i]) == 0 )
            test_all = true;
        else if ( strcmp("-c", argv[i]) == 0 )
            checkres = true;
        else if ( strcmp("-l", argv[i]) == 0 )
            lapack = true;
        else if ( strcmp("-w1", argv[i]) == 0 )
            workspace = 1;
        else if ( strcmp("-w2", argv[i]) == 0 )
            workspace = 2;
        else if ( strcmp("-w3", argv[i]) == 0 )
            workspace = 3;
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

    n2  = M * N;
    min_mn = min(M, N);

    /* Allocate host memory for the matrix */
    TESTING_MALLOC(h_A, float,  n2);
    TESTING_MALLOC( VT, float, N*N);
    TESTING_MALLOC(  U, float, M*M);
    TESTING_MALLOC( S1, float,       min_mn);
    TESTING_MALLOC( S2, float,       min_mn);

#if defined(PRECISION_z) || defined(PRECISION_c)
    TESTING_MALLOC(rwork, float,   5*min_mn);
#endif
    TESTING_HOSTALLOC(h_R, float, n2);

    magma_int_t nb = magma_get_sgesvd_nb(N);
    magma_int_t lwork;

    switch( workspace ) {
        default:
#if defined(PRECISION_z) || defined(PRECISION_c)
        case 1: lwork = (M+N)*nb + 2*min_mn;                   break;  // minimum
        case 2: lwork = (M+N)*nb + 2*min_mn +   min_mn*min_mn; break;  // optimal for some paths
        case 3: lwork = (M+N)*nb + 2*min_mn + 2*min_mn*min_mn; break;  // optimal for all paths
#else
        case 1: lwork = (M+N)*nb + 3*min_mn;                   break;  // minimum
        case 2: lwork = (M+N)*nb + 3*min_mn +   min_mn*min_mn; break;  // optimal for some paths
        case 3: lwork = (M+N)*nb + 3*min_mn + 2*min_mn*min_mn; break;  // optimal for all paths
#endif
    }

    TESTING_HOSTALLOC(h_work, float, lwork);
    
    const char* jobs[] = { "None", "Some", "Over", "All" };

    printf("-1.00 indicates non-applicable test that was skipped. See code for norm formulas.\n");
    printf("jobu jobv     M     N  CPU time (sec)  GPU time (sec)  |S1-S2|/.  |A-USV'|/. |I-UU'|/M  |I-VV'|/N  S (0=okay)\n");
    printf("===============================================================================================================\n");
    for( int i = 0; i < ntest; ++i ) {
        for( int ijobu = 0; ijobu < 4; ++ijobu ) {
        for( int ijobv = 0; ijobv < 4; ++ijobv ) {
            if ( test_all ) {
                jobu = jobs[ ijobu ];
                jobv = jobs[ ijobv ];
                if ( jobu[0] == 'O' && jobv[0] == 'O' ) {
                    // illegal combination; skip
                    continue;
                }
            }
            else if ( ijobu > 0 || ijobv > 0 ) {
                // if not testing all, run only once, with ijobu = ijobv = 0
                continue;
            }
            
            M = msize[i];
            N = nsize[i];
            n2 = M*N;
            min_mn = min(M, N);
    
            /* Initialize the matrix */
            lapackf77_slarnv( &ione, ISEED, &n2, h_A );
            lapackf77_slacpy( MagmaUpperLowerStr, &M, &N, h_A, &M, h_R, &M );
    
            /* ====================================================================
               Performs operation using MAGMA
               =================================================================== */
            gpu_time = magma_wtime();
            #if defined(PRECISION_z) || defined(PRECISION_c)
            magma_sgesvd( jobu[0], jobv[0], M, N,
                          h_R, M, S1, U, M,
                          VT, N, h_work, lwork, rwork, &info );
            #else
            magma_sgesvd( jobu[0], jobv[0], M, N,
                          h_R, M, S1, U, M,
                          VT, N, h_work, lwork, &info );
            #endif
            gpu_time = magma_wtime() - gpu_time;
            if (info != 0)
                printf("magma_sgesvd returned error %d.\n", (int) info);
            
            float eps = lapackf77_slamch( "E" );
            float result[4] = { -1/eps, -1/eps, -1/eps, -1/eps };
            if ( checkres ) {
                /* =====================================================================
                   Check the results following the LAPACK's [zcds]drvbd routine.
                   A is factored as A = U diag(S) VT and the following 4 tests computed:
                   (1)    | A - U diag(S) VT | / ( |A| max(M,N) )
                   (2)    | I - U'U | / ( M )
                   (3)    | I - VT VT' | / ( N )
                   (4)    S contains MNMIN nonnegative values in decreasing order.
                          (Return 0 if true, 1/ULP if false.)
                   =================================================================== */
                magma_int_t izero = 0;
                float *E;
                float *h_work_err;
                magma_int_t lwork_err = max(5*min_mn, (3*min_mn + max(M,N)))*128;
                TESTING_MALLOC(h_work_err, float, lwork_err);
                
                // get size and location of U and V^T depending on jobu and jobv
                // U2=NULL and VT2=NULL if they were not computed (e.g., jobu=N)
                magma_int_t M2  = (jobu[0] == 'A' ? M : min_mn);
                magma_int_t N2  = (jobv[0] == 'A' ? N : min_mn);
                magma_int_t ldu = M;
                magma_int_t ldv = (jobv[0] == 'O' ? M : N);
                float *U2  = NULL;
                float *VT2 = NULL;
                if ( jobu[0] == 'S' || jobu[0] == 'A' ) {
                    U2 = U;
                } else if ( jobu[0] == 'O' ) {
                    U2 = h_R;
                }
                if ( jobv[0] == 'S' || jobv[0] == 'A' ) {
                    VT2 = VT;
                } else if ( jobv[0] == 'O' ) {
                    VT2 = h_R;
                }
                
                #if defined(PRECISION_z) || defined(PRECISION_c)
                if ( U2 != NULL && VT2 != NULL ) {
                    lapackf77_sbdt01(&M, &N, &izero, h_A, &M,
                                     U2, &ldu, S1, E, VT2, &ldv, h_work_err, rwork, &result[0]);
                }
                if ( U2 != NULL ) {
                    lapackf77_sort01("Columns", &M, &M2, U2,  &ldu, h_work_err, &lwork_err, rwork, &result[1]);
                }
                if ( VT2 != NULL ) {
                    lapackf77_sort01(   "Rows", &N2, &N, VT2, &ldv, h_work_err, &lwork_err, rwork, &result[2]);
                }
                #else
                if ( U2 != NULL && VT2 != NULL ) {
                    lapackf77_sbdt01(&M, &N, &izero, h_A, &M,
                                     U2, &ldu, S1, E, VT2, &ldv, h_work_err, &result[0]);
                }
                if ( U2 != NULL ) {
                    lapackf77_sort01("Columns", &M, &M2, U2,  &ldu, h_work_err, &lwork_err, &result[1]);
                }
                if ( VT2 != NULL ) {
                    // this step may be really slow for large N
                    lapackf77_sort01(   "Rows", &N2, &N, VT2, &ldv, h_work_err, &lwork_err, &result[2]);
                }
                #endif
                
                result[3] = 0.;
                for(int j=0; j < min_mn-1; j++){
                    if ( S1[j] < S1[j+1] )
                        result[3] = 1.;
                    if ( S1[j] < 0. )
                        result[3] = 1.;
                }
                if (min_mn > 1 && S1[min_mn-1] < 0.)
                    result[3] = 1.;
                
                result[0] *= eps;
                result[1] *= eps;
                result[2] *= eps;
                
                TESTING_FREE( h_work_err );
            }
    
            /* =====================================================================
               Performs operation using LAPACK
               =================================================================== */
            if ( lapack ) {
                cpu_time = magma_wtime();
                #if defined(PRECISION_z) || defined(PRECISION_c)
                lapackf77_sgesvd( jobu, jobv, &M, &N,
                                  h_A, &M, S2, U, &M,
                                  VT, &N, h_work, &lwork, rwork, &info);
                #else
                lapackf77_sgesvd( jobu, jobv, &M, &N,
                                  h_A, &M, S2, U, &M,
                                  VT, &N, h_work, &lwork, &info);
                #endif
                cpu_time = magma_wtime() - cpu_time;
                if (info != 0)
                    printf("lapackf77_sgesvd returned error %d.\n", (int) info);
                
                /* =====================================================================
                   Check the result compared to LAPACK
                   =================================================================== */
                float work[1], error = 1., mone = -1;
                magma_int_t one = 1;
        
                error = lapackf77_slange("f", &min_mn, &one, S1, &min_mn, work);
                blasf77_saxpy(&min_mn, &mone, S1, &one, S2, &one);
                error = lapackf77_slange("f", &min_mn, &one, S2, &min_mn, work) / error;
                
                printf("   %c    %c %5d %5d  %7.2f         %7.2f         %8.2e",
                       jobu[0], jobv[0], (int) M, (int) N, cpu_time, gpu_time, error );
            }
            else {
                printf("   %c    %c %5d %5d    ---           %7.2f         ---",
                       jobu[0], jobv[0], (int) M, (int) N, gpu_time );
            }
            if ( checkres ) {
                printf("  %#9.3g  %#9.3g  %#9.3g   %1.0f\n",
                       result[0], result[1], result[2], result[3] );
            }
            else {
                printf("\n");
            }
        }}
        if ( test_all ) {
            printf("\n");
        }
    }

    /* Memory clean up */
    TESTING_FREE(       h_A);
    TESTING_FREE(        VT);
    TESTING_FREE(        S1);
    TESTING_FREE(        S2);
#if defined(PRECISION_z) || defined(PRECISION_c)
    TESTING_FREE(     rwork);
#endif
    TESTING_FREE(         U);
    TESTING_HOSTFREE(h_work);
    TESTING_HOSTFREE(   h_R);

    TESTING_CUDA_FINALIZE();
    return 0;
}
