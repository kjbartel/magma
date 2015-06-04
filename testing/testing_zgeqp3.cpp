/*
    -- MAGMA (version 1.2.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       June 2012

       @precisions normal z -> c

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
#define PRECISION_z
#if defined(PRECISION_z) || defined(PRECISION_c)
#define FLOPS(m, n) ( 6.*FMULS_GEQRF(m, n) + 2.*FADDS_GEQRF(m, n) )
#else
#define FLOPS(m, n) (    FMULS_GEQRF(m, n) +    FADDS_GEQRF(m, n) )
#endif

/* ////////////////////////////////////////////////////////////////////////////
   -- Testing zgeqp3
*/
int main( int argc, char** argv) 
{
    TESTING_CUDA_INIT();

    magma_timestr_t       start, end;
    magma_int_t      checkres;
    double           flops, gpu_perf, cpu_perf, gpu_time, cpu_time;
    double           matnorm, work[1];
    cuDoubleComplex  c_neg_one = MAGMA_Z_NEG_ONE;
    cuDoubleComplex *h_A, *h_R, *tau, *h_work, tmp[1];

    double *rwork;
    magma_int_t *jpvt;

    /* Matrix size */
    magma_int_t M = 0, N = 0, n2, lda, lwork;
    magma_int_t size[10] = {1024,2048,3072,4032,5184,6016,7040,8064,9088,9984};

    magma_int_t i, j, info, min_mn, nb;
    magma_int_t ione     = 1;
    magma_int_t ISEED[4] = {0,0,0,1};

    if (argc != 1){
        for(i = 1; i<argc; i++){
            if (strcmp("-N", argv[i])==0)
                N = atoi(argv[++i]);
            else if (strcmp("-M", argv[i])==0)
                M = atoi(argv[++i]);
        }
        if ( M == 0 ) {
            M = N;
        }
        if ( N == 0 ) {
            N = M;
        }
        if (N>0 && M>0)
            printf("  testing_zgeqp3 -M %d -N %d\n\n", M, N);
        else
            {
                printf("\nUsage: \n");
                printf("  testing_zgeqp3 -M %d -N %d\n\n", M, N);
                exit(1);
            }
    }
    else {
        printf("\nUsage: \n");
        printf("  testing_zgeqp3 -M %d -N %d\n\n", 1024, 1024);
        M = N = size[9];
    }

    checkres = getenv("MAGMA_TESTINGS_CHECK") != NULL;
    checkres = 1;

    n2  = M * N;
    min_mn = min(M, N);
    nb = magma_get_zgeqp3_nb(min_mn);

    TESTING_MALLOC(    jpvt, magma_int_t,     N );

    TESTING_MALLOC(    tau,  cuDoubleComplex, min_mn);
    TESTING_MALLOC(    rwork, double, 2*N);
    TESTING_MALLOC(    h_A,  cuDoubleComplex, n2 );
    TESTING_HOSTALLOC( h_R,  cuDoubleComplex, n2 );
    //TESTING_MALLOC( h_R,  cuDoubleComplex, n2   );

    lwork = ( N+1 )*nb;
    if ( checkres )
        lwork = max(lwork, M * N + N);
    TESTING_MALLOC( h_work, cuDoubleComplex, lwork );

    printf("  M     N   CPU GFlop/s (sec)   GPU GFlop/s (sec)  ||A*P - Q*R||_F  \n");
    printf("==================================================================\n");
    for(i=0; i<10; i++){
        if (argc == 1){
            M = N = size[i];
        }
        min_mn= min(M, N);
        lda   = M;
        n2    = lda*N;
        flops = FLOPS( (double)M, (double)N ) / 1000000;

        /* Initialize the matrix */
        lapackf77_zlarnv( &ione, ISEED, &n2, h_A );
        lapackf77_zlacpy( MagmaUpperLowerStr, &M, &N, h_A, &lda, h_R, &lda );

        /* =====================================================================
           Performs operation using LAPACK
           =================================================================== */
        for (j = 0; j < N; j++)
            jpvt[j] = 0;

        start = get_current_time();
        lapackf77_zgeqp3(&M, &N, h_R, &lda, jpvt, tau, h_work, &lwork, rwork, &info);
        end = get_current_time();
        if (info < 0)
            printf("Argument %d of lapack_zgeqp3 had an illegal value.\n", -info);

        cpu_perf = flops / GetTimerValue(start, end);
        cpu_time = GetTimerValue(start, end) * 1e-3;

        /* ====================================================================
           Performs operation using MAGMA
           =================================================================== */
        lapackf77_zlacpy( MagmaUpperLowerStr, &M, &N, h_A, &lda, h_R, &lda );
        for (j = 0; j < N; j++) 
            jpvt[j] = 0 ;

        start = get_current_time();
        magma_zgeqp3(M, N, h_R, lda, jpvt, tau, h_work, lwork, rwork, &info);
        end = get_current_time();
        if (info < 0)
            printf("Argument %d of magma_zgeqp3 had an illegal value.\n", -info);
        
        gpu_perf = flops / GetTimerValue(start, end);
        gpu_time = GetTimerValue(start, end) * 1e-3;

        /* =====================================================================
           Check the result 
           =================================================================== */
        if ( checkres ) 
            {
                double result[3], ulp;

                magma_int_t minmn = min(M, N);
                ulp = lapackf77_dlamch( "P" );
                
                // Compute norm( A*P - Q*R )
                result[0] = lapackf77_zqpt01(&M, &N, &minmn, h_A, h_R, &lda, 
                                             tau, jpvt, h_work, &lwork);
                result[0] *= ulp;

                printf("%5d %5d  %6.2f (%6.2f)    %6.2f (%6.2f)     %e\n",
                       M, N, cpu_perf, cpu_time, gpu_perf, gpu_time, result[0]);
            }
        else
            printf("%5d %5d  %6.2f (%6.2f)    %6.2f (%6.2f)\n",
                   M, N, cpu_perf, cpu_time, gpu_perf, gpu_time);

        if (argc != 1)
            break;
    }

    /* Memory clean up */
    TESTING_FREE( jpvt );
    TESTING_FREE( tau );
    TESTING_FREE( rwork );
    TESTING_FREE( h_A );
    TESTING_HOSTFREE( h_R );
    //TESTING_FREE( h_R );
    TESTING_FREE( h_work );

    /* Shutdown */
    TESTING_CUDA_FINALIZE();
}
