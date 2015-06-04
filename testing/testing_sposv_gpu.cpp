/*
 *  -- MAGMA (version 1.1) --
 *     Univ. of Tennessee, Knoxville
 *     Univ. of California, Berkeley
 *     Univ. of Colorado, Denver
 *     November 2011
 *
 * @generated s Sun Nov 13 20:48:51 2011
 *
 **/
// includes, system
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cublas.h>
#include <unistd.h>

// includes, project
#include "flops.h"
#include "magma.h"
#include "magma_lapack.h"
#include "testings.h"

#define PRECISION_s
// Flops formula
#if defined(PRECISION_z) || defined(PRECISION_c)
#define FLOPS_POTRF(n      ) ( 6.*FMULS_POTRF(n      ) + 2.*FADDS_POTRF(n      ) )
#define FLOPS_POTRS(n, nrhs) ( 6.*FMULS_POTRS(n, nrhs) + 2.*FADDS_POTRS(n, nrhs) )
#else
#define FLOPS_POTRF(n      ) (    FMULS_POTRF(n      ) +    FADDS_POTRF(n      ) )
#define FLOPS_POTRS(n, nrhs) (    FMULS_POTRS(n, nrhs) +    FADDS_POTRS(n, nrhs) )
#endif

/* ////////////////////////////////////////////////////////////////////////////
   -- Testing sposv_gpu
*/
int main( int argc, char** argv)
{
    TESTING_CUDA_INIT();

    magma_timestr_t start, end;
    float          flops, gpu_perf;
    float          Rnorm, Anorm, Bnorm, *work;
    float zone  = MAGMA_S_ONE;
    float mzone = MAGMA_S_NEG_ONE;
    float *h_A, *h_B, *h_X;
    float *d_A, *d_B;
    const char  *uplo     = MagmaLowerStr;
    magma_int_t lda, ldb, N;
    magma_int_t ldda, lddb;
    magma_int_t i, info, szeA, szeB;
    magma_int_t ione     = 1;
    magma_int_t NRHS     = 100;
    magma_int_t ISEED[4] = {0,0,0,1};
    const int MAXTESTS   = 10;
    magma_int_t size[MAXTESTS] = {1024,2048,3072,4032,5184,6016,7040,8064,9088,10112};

    // process command line arguments
    printf( "\nUsage:\n" );
    printf( "  %s -N <matrix size> -R <right hand sides>\n", argv[0] );
    printf( "  -N can be repeated up to %d times\n", MAXTESTS );
    int ntest = 0;
    int ch;
    while( (ch = getopt( argc, argv, "N:R:" )) != -1 ) {
        switch( ch ) {
            case 'N':
                if ( ntest == MAXTESTS ) {
                    printf( "error: -N exceeded maximum %d tests\n", MAXTESTS );
                    exit(1);
                }
                else {
                    size[ntest] = atoi( optarg );
                    if ( size[ntest] <= 0 ) {
                        printf( "error: -N value %d <= 0\n", size[ntest] );
                        exit(1);
                    }
                    ntest++;
                }
                break;
            case 'R':
                NRHS = atoi( optarg );
                break;
            case '?':
            default:
                exit(1);
        }
    }
    argc -= optind;
    argv += optind;
    if ( ntest == 0 ) {
        ntest = MAXTESTS;
    }
    
    // allocate maximum amount of memory required
    N = 0;
    for( i = 0; i < ntest; ++i ) {
        N = max( N, size[i] );
    }
    lda = ldb = N;
    lddb = ldda = ((N+31)/32)*32;
    
    TESTING_MALLOC( h_A, float, lda*N    );
    TESTING_MALLOC( h_B, float, ldb*NRHS );
    TESTING_MALLOC( h_X, float, ldb*NRHS );
    TESTING_MALLOC( work, float,         N        );

    TESTING_DEVALLOC( d_A, float, ldda*N    );
    TESTING_DEVALLOC( d_B, float, lddb*NRHS );

    printf("\n\n");
    printf("  N     NRHS       GPU GFlop/s      || b-Ax || / ||A||*||B||\n");
    printf("========================================================\n");
    
    for( i = 0; i < ntest; ++i ) {
        N   = size[i];
        lda = ldb = N;
        ldda = ((N+31)/32)*32;
        lddb = ldda;
        flops = ( FLOPS_POTRF( (float)N ) +
                  FLOPS_POTRS( (float)N, (float)NRHS ) ) / 1e6;

        /* ====================================================================
           Initialize the matrix
           =================================================================== */
        szeA = lda*N;
        szeB = ldb*NRHS;
        lapackf77_slarnv( &ione, ISEED, &szeA, h_A );
        lapackf77_slarnv( &ione, ISEED, &szeB, h_B );
        /* Symmetrize and increase the diagonal */
        {
            magma_int_t i, j;
            for(i=0; i<N; i++) {
                MAGMA_S_SET2REAL( h_A[i*lda+i], ( MAGMA_S_GET_X(h_A[i*lda+i]) + 1.*N ) );
                for(j=0; j<i; j++)
                    h_A[i*lda+j] = (h_A[j*lda+i]);
            }
        }
        
        cublasSetMatrix( N, N,    sizeof( float ), h_A, N, d_A, ldda );
        cublasSetMatrix( N, NRHS, sizeof( float ), h_B, N, d_B, lddb );

        /* ====================================================================
           Performs operation using MAGMA
           =================================================================== */
        start = get_current_time();
        magma_sposv_gpu( uplo[0], N, NRHS, d_A, ldda, d_B, lddb, &info );
        end = get_current_time();
        if (info != 0)
            printf("Argument %d of magma_spotrf had an illegal value.\n", -info);

        gpu_perf = flops / GetTimerValue(start, end);

        /* =====================================================================
           Residual
           =================================================================== */
        cublasGetMatrix( N, NRHS, sizeof( float ), d_B, lddb, h_X, ldb );
        
        Anorm = lapackf77_slange("I", &N, &N,    h_A, &lda, work);
        Bnorm = lapackf77_slange("I", &N, &NRHS, h_B, &ldb, work);

        blasf77_sgemm( MagmaNoTransStr, MagmaNoTransStr, &N, &NRHS, &N,
                       &zone,  h_A, &lda,
                               h_X, &ldb,
                       &mzone, h_B, &ldb );
        
        Rnorm = lapackf77_slange("I", &N, &NRHS, h_B, &ldb, work);

        printf("%5d  %4d             %6.2f        %e\n",
               N, NRHS, gpu_perf, Rnorm/(Anorm*Bnorm) );
    }

    /* Memory clean up */
    TESTING_FREE( h_A );
    TESTING_FREE( h_B );
    TESTING_FREE( h_X );
    TESTING_FREE( work );

    TESTING_DEVFREE( d_A );
    TESTING_DEVFREE( d_B );

    TESTING_CUDA_FINALIZE();
}
