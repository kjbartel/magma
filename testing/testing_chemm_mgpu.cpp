/*
    -- MAGMA (version 1.2.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       June 2012

       @generated c Thu Jun 28 12:31:52 2012
       
       @author Mark Gates
*/

// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <cuda_runtime_api.h>
#include <cublas.h>
#include <assert.h>
#include <stdarg.h>  // va_start

// includes, project
#include "flops.h"
#include "magma.h"
#include "magma_lapack.h"
#include "testings.h"


// --------------------
void ensure( bool condition, const char* msg, ... )
{
    if ( not condition ) {
        va_list va;
        va_start( va, msg );
        vprintf( msg, va );
        exit(1);
    }
}


/* ////////////////////////////////////////////////////////////////////////////
   -- Testing magma_chemm_mgpu
*/
int main( int argc, char** argv)
{
    TESTING_CUDA_INIT();

    cuFloatComplex c_neg_one = MAGMA_C_NEG_ONE;
    cuFloatComplex c_one     = MAGMA_C_ONE;
    float d_one = 1.0;
    
    real_Double_t    gflops, gpu_perf=0., cpu_perf=0., gpu_time=0., cpu_time=0.;
    real_Double_t    gpu_perf2, gpu_time2;
    float           error=0., error2=0., work[1];
    cuFloatComplex *hA, *hX, *hB, *hR;
    cuFloatComplex *dA[MagmaMaxGPUs], *dX[MagmaMaxGPUs], *dB[MagmaMaxGPUs];
    cuFloatComplex *dA2;
    
    /* Matrix size */
    magma_int_t m, size, lda, ldda;
    const int MAXTESTS = 10;
    // sizes are 1024*N - 32
    magma_int_t msize[MAXTESTS] =
        { 992, 2016, 3040, 4064, 5088, 6112, 7136, 8160, 9184, 10208 };
    magma_int_t n       = 64;
    magma_int_t nb      = 64;
    int nstream = 2;
    int count   = 3;
    int ngpu    = magma_num_gpus();
    
    magma_int_t info;
    magma_int_t ione     = 1;
    magma_int_t iseed[4] = {0,0,0,1};
        
    printf( "Usage: %s -M m -N n -nb nb -nstream nstream -ngpu ngpu -count count -c\n"
            "    -M can be repeated %d times.\n"
            "    -ngpu or setting $MAGMA_NUM_GPUS sets number of GPUs to use.\n"
            "    -c or setting $MAGMA_TESTINGS_CHECK runs LAPACK and checks result.\n",
            argv[0], MAXTESTS );

    int checkres = (getenv("MAGMA_TESTINGS_CHECK") != NULL);

    int ntest = 0;
    int mmax = 0;
    for( int i = 1; i < argc; i++ ) {
        if ( strcmp("-M", argv[i]) == 0 and i+1 < argc ) {
            ensure( ntest < MAXTESTS, "error: -M repeated more than maximum %d tests\n", MAXTESTS );
            msize[ntest] = atoi( argv[++i] );
            ensure( msize[ntest] > 0, "error: -M %s is invalid; must be > 0.\n", argv[i] );
            mmax = max( mmax, msize[ntest] );
            ntest++;
        }
        else if ( strcmp("-N", argv[i]) == 0 and i+1 < argc ) {
            n = atoi( argv[++i] );
            ensure( n > 0, "error: -N %s is invalid; must be > 0.\n", argv[i] );
        }
        else if ( strcmp("-nb", argv[i]) == 0 and i+1 < argc ) {
            nb = atoi( argv[++i] );
            ensure( nb > 0, "error: -nb %s is invalid; must be > 0.\n", argv[i] );
        }
        else if ( strcmp("-count", argv[i]) == 0 and i+1 < argc ) {
            count = atoi( argv[++i] );
            ensure( count > 0, "error: -count %s is invalid; must be > 0.\n", argv[i] );
        }
        else if ( strcmp("-nstream", argv[i]) == 0 and i+1 < argc ) {
            nstream = atoi( argv[++i] );
            ensure( nstream > 0 and nstream <= 20,
                    "error: -nstream %s is invalid; must be > 0 and <= 20.\n", argv[i] );
        }
        else if ( strcmp("-ngpu", argv[i]) == 0 and i+1 < argc ) {
            ngpu = atoi( argv[++i] );
            ensure( ngpu > 0, "error: -ngpu %s is invalid; must be > 0.\n", argv[i] );
        }
        else if ( strcmp("-c", argv[i]) == 0 ) {
            checkres = true;
        }
        else {
            printf( "invalid argument: %s\n", argv[i] );
            exit(1);
        }
    }
    if ( ntest == 0 ) {
        ntest = MAXTESTS;
        mmax = msize[ntest-1];
    }
    assert( mmax > 0 and n > 0 );
    
    // allocate memory for largest problem
    m = mmax;
    lda  = m;
    ldda = ((m + 31)/32)*32;

    TESTING_MALLOC( hA, cuFloatComplex, lda*m );
    TESTING_MALLOC( hX, cuFloatComplex, lda*n );
    TESTING_MALLOC( hB, cuFloatComplex, lda*n );
    TESTING_MALLOC( hR, cuFloatComplex, lda*n*2 );
    
    cudaStream_t streams[MagmaMaxGPUs][20];    
    for( int d = 0; d < ngpu; ++d ) {
        magma_int_t mlocal = ((m / nb) / ngpu + 1) * nb;
        cudaSetDevice( d );
        TESTING_DEVALLOC( dA[d], cuFloatComplex, ldda*mlocal );
        TESTING_DEVALLOC( dX[d], cuFloatComplex, ldda*n      );
        TESTING_DEVALLOC( dB[d], cuFloatComplex, ldda*n      );
        for( int i = 0; i < nstream; ++i ) {
            cudaStreamCreate( &streams[d][i] );
        }
    }
    if ( checkres ) {
    cudaSetDevice( 0 );
    TESTING_DEVALLOC( dA2, cuFloatComplex, ldda*m );
    }
    
    printf("\n");
    printf( "nb %d, ngpu %d, nstream %d\n", (int) nb, ngpu, nstream );
    printf("    m     n    CPU GFlop/s (sec)   GPU GFlop/s (sec)   CUBLAS hemm (sec)   ||R|| / ||A||*||X||\n");
    printf("==============================================================================================\n");
    for( int i = 0; i < ntest; ++i ) {
    for( int j = 0; j < count; ++j ) {
        m = msize[i];
        assert( m > 0 and n > 0 );
        
        lda  = m;
        ldda = ((m + 31)/32)*32;
        gflops = FLOPS_CHEMM( MagmaLeft, (float)m, (float)n ) / 1e9;

        size = lda*m;
        lapackf77_clarnv( &ione, iseed, &size, hA );
        // make diagonal real
        for( int i = 0; i < m; ++i ) {
            hA[i + i*lda] = MAGMA_C_MAKE( MAGMA_C_REAL( hA[i+i*lda] ), 0. );
        }
        size = lda*n;
        lapackf77_clarnv( &ione, iseed, &size, hX );
        lapackf77_clarnv( &ione, iseed, &size, hB );
        lapackf77_clacpy( "Full", &m, &n, hB, &lda, hR, &lda );
        
        /* ====================================================================
           Performs operation using MAGMA
           =================================================================== */
        magmablas_csetmatrix_1D_bcyclic( m, m, hA, lda, dA, ldda, ngpu, nb );
        for( int d = 0; d < ngpu; ++d ) {
            cudaSetDevice( d );
            magma_csetmatrix( m, n, hX, lda, dX[d], ldda );
            //magma_csetmatrix( m, n, hB, lda, dB[d], ldda );
        }
        
        cudaDeviceSynchronize();
        gpu_time = magma_wtime();
        magmablas_chemm_mgpu(
            MagmaLeft, MagmaLower, m, n,
            c_neg_one, dA, ldda, 0,
                       dX, ldda,
            c_one,     dB, ldda, hR, lda,
            ngpu, nb, streams, nstream );
        cudaDeviceSynchronize();
        gpu_time = magma_wtime() - gpu_time;
                
        gpu_perf = gflops / gpu_time;
        
        /* ====================================================================
           Performs operation using CUBLAS
           =================================================================== */
        if ( checkres ) {
            magma_setdevice( 0 );
            magma_csetmatrix( m, m, hA, lda, dA2, ldda );
            magma_csetmatrix( m, n, hX, lda, dX[0], ldda );
            magma_csetmatrix( m, n, hB, lda, dB[0], ldda );
            
            real_Double_t gpu_time2, gpu_perf2;
            cudaDeviceSynchronize();
            gpu_time2 = magma_wtime();
            magma_chemm(
                MagmaLeft, MagmaLower, m, n,
                c_neg_one, dA2,   ldda,
                           dX[0], ldda,
                c_one,     dB[0], ldda );
            cudaDeviceSynchronize();
            gpu_time2 = magma_wtime() - gpu_time2;
            
            gpu_perf2 = gflops / gpu_time2;
        }
        
        /* =====================================================================
           Performs operation using LAPACK
           =================================================================== */
        if ( checkres ) {
            // store ||A||*||X||
            error  = lapackf77_clange("fro", &m, &m, hA, &lda, work );
            error *= lapackf77_clange("fro", &m, &n, hX, &lda, work );
            
            //printf( "A =" ); magma_cprint( m, m, hA, lda );
            //printf( "X =" ); magma_cprint( m, n, hX, lda );
            //printf( "B =" ); magma_cprint( m, n, hB, lda );
            
            cpu_time = magma_wtime();
            blasf77_chemm( "Left", "Lower", &m, &n,
                            &c_neg_one, hA, &lda,
                                        hX, &lda,
                            &c_one,     hB, &lda );
            cpu_time = magma_wtime() - cpu_time;
            
            //printf( "B  =" ); magma_cprint( m, n, hB, lda );
            //printf( "dB =" ); magma_cprint( m, n, hR, lda );
            
            cpu_perf = gflops / cpu_time;
    
            // compute relative error ||R||/||A||*||X||, where R := B_magma - B_lapack = R - B
            size = lda*n;
            blasf77_caxpy( &size, &c_neg_one, hB, &ione, hR, &ione );
            error = lapackf77_clange("fro", &m, &n, hR, &lda, work) / error;
            
            //printf( "R ="  ); magma_cprint( m, n, hR, lda );
            
            printf( "%5d %5d   %7.1f (%7.4f)   %7.1f (%7.4f)   %7.1f (%7.4f)   %8.2e\n",
                    (int) m, (int) n, cpu_perf, cpu_time, gpu_perf, gpu_time, gpu_perf2, gpu_time2, error );
        }
        else {
            printf( "%5d %5d     ---   (  ---  )   %7.1f (%7.4f)   %7.1f (%7.4f)   ---\n",
                    (int) m, (int) n, /*cpu_perf, cpu_time,*/ gpu_perf, gpu_time, gpu_perf2, gpu_time2 /*, error*/ );
        }
    }}
    
    /* Memory clean up */
    for( int d = 0; d < ngpu; ++d ) {
        cudaSetDevice( d );
        TESTING_DEVFREE( dA[d] );
        TESTING_DEVFREE( dX[d] );
        //TESTING_DEVFREE( dB[d] );
    }
    
    TESTING_FREE( hA );
    TESTING_FREE( hX );
    TESTING_FREE( hB );
    TESTING_FREE( hR );
    
    /* Shutdown */
    TESTING_CUDA_FINALIZE();
    return 0;
}
