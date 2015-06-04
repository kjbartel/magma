/*
    -- MAGMA (version 1.2.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       May 2012

       @generated c Tue May 15 18:18:22 2012

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
#include "flops.h"
#include "magma.h"
#include "magma_lapack.h"
#include "testings.h"

// Flops formula
#define PRECISION_c
#if defined(PRECISION_z) || defined(PRECISION_c)
#define FLOPS_GEQRF(m, n      ) ( 6.*FMULS_GEQRF(m, n      ) + 2.*FADDS_GEQRF(m, n      ) )
#define FLOPS_GEQRS(m, n, nrhs) ( 6.*FMULS_GEQRS(m, n, nrhs) + 2.*FADDS_GEQRS(m, n, nrhs) )
#else
#define FLOPS_GEQRF(m, n      ) (    FMULS_GEQRF(m, n      ) +    FADDS_GEQRF(m, n      ) )
#define FLOPS_GEQRS(m, n, nrhs) (    FMULS_GEQRS(m, n, nrhs) +    FADDS_GEQRS(m, n, nrhs) )
#endif

/* ////////////////////////////////////////////////////////////////////////////
   -- Testing cgeqrs
*/
int main( int argc, char** argv)
{
    TESTING_CUDA_INIT();
   
    magma_timestr_t       start, end;
    float           flops, gpu_perf, cpu_perf;
    float           matnorm, work[1];
    cuFloatComplex  c_one     = MAGMA_C_ONE;
    cuFloatComplex  c_neg_one = MAGMA_C_NEG_ONE;
    cuFloatComplex *h_A, *h_A2, *h_B, *h_X, *h_R, *tau, *hwork, tmp[1];
    cuFloatComplex *d_A, *d_B;

    /* Matrix size */
    magma_int_t M = 0, N = 0, n2;
    magma_int_t lda, ldb, ldda, lddb, lworkgpu, lhwork;
    magma_int_t size[10] = {1024,2048,3072,4032,5184,6016,7040,8064,9088,10112};

    magma_int_t i, info, min_mn, nb, l1, l2;
    magma_int_t ione     = 1;
    magma_int_t nrhs     = 3;
    magma_int_t ISEED[4] = {0,0,0,1};

    if (argc != 1){
        for(i = 1; i<argc; i++){
            if (strcmp("-N", argv[i])==0)
                N = atoi(argv[++i]);
            else if (strcmp("-M", argv[i])==0)
                M = atoi(argv[++i]);
            else if (strcmp("-nrhs", argv[i])==0)
                nrhs = atoi(argv[++i]);
        }
        if (N>0 && M>0 && M >= N)
            printf("  testing_cgeqrs_gpu -nrhs %d -M %d -N %d\n\n", nrhs, M, N);
        else
            {
                printf("\nUsage: \n");
                printf("  testing_cgeqrs_gpu -nrhs %d  -M %d  -N %d\n\n", nrhs, M, N);
                printf("  M has to be >= N, exit.\n");
                exit(1);
            }
    }
    else {
        printf("\nUsage: \n");
        printf("  testing_cgeqrs_gpu -nrhs %d  -M %d  -N %d\n\n", nrhs, 1024, 1024);
        M = N = size[9];
    }

    ldda   = ((M+31)/32)*32;
    lddb   = ldda;
    n2     = M * N;
    min_mn = min(M, N);
    nb     = magma_get_cgeqrf_nb(M);
    lda = ldb = M;
    lworkgpu = (M-N + nb)*(nrhs+2*nb);

    /* Allocate host memory for the matrix */
    TESTING_MALLOC( tau,  cuFloatComplex, min_mn   );
    TESTING_MALLOC( h_A,  cuFloatComplex, lda*N    );
    TESTING_MALLOC( h_A2, cuFloatComplex, lda*N    );
    TESTING_MALLOC( h_B,  cuFloatComplex, ldb*nrhs );
    TESTING_MALLOC( h_X,  cuFloatComplex, ldb*nrhs );
    TESTING_MALLOC( h_R,  cuFloatComplex, ldb*nrhs );

    TESTING_DEVALLOC( d_A, cuFloatComplex, ldda*N      );
    TESTING_DEVALLOC( d_B, cuFloatComplex, lddb*nrhs   );

    /*
     * Get size for host workspace
     */
    lhwork = -1;
    lapackf77_cgeqrf(&M, &N, h_A, &M, tau, tmp, &lhwork, &info);
    l1 = (magma_int_t)MAGMA_C_REAL( tmp[0] );
    lhwork = -1;
    lapackf77_cunmqr( MagmaLeftStr, MagmaConjTransStr, 
                      &M, &nrhs, &min_mn, h_A, &lda, tau,
                      h_X, &ldb, tmp, &lhwork, &info);
    l2 = (magma_int_t)MAGMA_C_REAL( tmp[0] );
    lhwork = max( max( l1, l2 ), lworkgpu );

    TESTING_MALLOC( hwork, cuFloatComplex, lhwork );

    printf("\n");
    printf("                                         ||b-Ax|| / (N||A||)\n");
    printf("  M     N    CPU GFlop/s   GPU GFlop/s      CPU      GPU    \n");
    printf("============================================================\n");
    for(i=0; i<10; i++){
        if (argc == 1){
            M = N = size[i];
        }
        min_mn= min(M, N);
        ldb = lda = M;
        n2    = lda*N;
        ldda  = ((M+31)/32)*32;
        flops = (FLOPS_GEQRF( (float)M, (float)N ) 
                 + FLOPS_GEQRS( (float)M, (float)N, (float)nrhs )) / 1000000;

        /* Initialize the matrices */
        lapackf77_clarnv( &ione, ISEED, &n2, h_A );
        lapackf77_clacpy( MagmaUpperLowerStr, &M, &N, h_A, &lda, h_A2, &lda );

        n2 = M*nrhs;
        lapackf77_clarnv( &ione, ISEED, &n2, h_B );
        lapackf77_clacpy( MagmaUpperLowerStr, &M, &nrhs, h_B, &ldb, h_R, &ldb );

        /* ====================================================================
           Performs operation using MAGMA
           =================================================================== */
        magma_csetmatrix( M, N,    h_A, lda, d_A, ldda );
        magma_csetmatrix( M, nrhs, h_B, ldb, d_B, lddb );

        start = get_current_time();
        magma_cgels3_gpu( MagmaNoTrans, M, N, nrhs, d_A, ldda,
                          d_B, lddb, hwork, lworkgpu, &info);
        end = get_current_time();
        if (info < 0)
            printf("Argument %d of magma_cgels had an illegal value.\n", -info);
        
        gpu_perf = flops / GetTimerValue(start, end);

        // Get the solution in h_X
        magma_cgetmatrix( N, nrhs, d_B, lddb, h_X, ldb );

        // compute the residual
        blasf77_cgemm( MagmaNoTransStr, MagmaNoTransStr, &M, &nrhs, &N, 
                       &c_neg_one, h_A, &lda, 
                                   h_X, &ldb, 
                       &c_one,     h_R, &ldb);
        matnorm = lapackf77_clange("f", &M, &N, h_A, &lda, work);

        /* =====================================================================
           Performs operation using LAPACK
           =================================================================== */
        lapackf77_clacpy( MagmaUpperLowerStr, &M, &nrhs, h_B, &ldb, h_X, &ldb );

        start = get_current_time();
        lapackf77_cgels( MagmaNoTransStr, &M, &N, &nrhs,
                         h_A, &lda, h_X, &ldb, hwork, &lhwork, &info);
        end = get_current_time();
        cpu_perf = flops / GetTimerValue(start, end);
        if (info < 0)
          printf("Argument %d of lapackf77_cgels had an illegal value.\n", -info);

        blasf77_cgemm( MagmaNoTransStr, MagmaNoTransStr, &M, &nrhs, &N, 
                       &c_neg_one, h_A2, &lda, 
                                   h_X,  &ldb, 
                       &c_one,     h_B,  &ldb);

        printf("%5d %5d   %6.1f       %6.1f       %7.2e   %7.2e\n",
               M, N, cpu_perf, gpu_perf,
               lapackf77_clange("f", &M, &nrhs, h_B, &M, work)/(min_mn*matnorm),
               lapackf77_clange("f", &M, &nrhs, h_R, &M, work)/(min_mn*matnorm) );

        if (argc != 1)
            break;
    }

    /* Memory clean up */
    TESTING_FREE( tau );
    TESTING_FREE( h_A );
    TESTING_FREE( h_A2 );
    TESTING_FREE( h_B );
    TESTING_FREE( h_X );
    TESTING_FREE( h_R );
    TESTING_FREE( hwork );
    TESTING_DEVFREE( d_A );
    TESTING_DEVFREE( d_B );

    /* Shutdown */
    TESTING_CUDA_FINALIZE();
}
