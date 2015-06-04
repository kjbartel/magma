/*
    -- MAGMA (version 1.2.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       June 2012

       @precisions normal z -> s d c
       
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

// version from zhetrd_mgpu.cpp
extern "C" void
magma_zher2k_mgpu(
    int num_gpus, char uplo, char trans, int nb, int n, int k,
    cuDoubleComplex alpha, cuDoubleComplex **db, int lddb, 
    double beta,           cuDoubleComplex **dc, int lddc, int offset,
    int num_streams, cudaStream_t stream[][10]);


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
   -- Testing magma_zher2k_mgpu
*/
int main( int argc, char** argv)
{
    TESTING_CUDA_INIT();

    cuDoubleComplex c_neg_one = MAGMA_Z_NEG_ONE;
    cuDoubleComplex c_one     = MAGMA_Z_ONE;
    double d_one = 1.0;
    
    real_Double_t    gflops, gpu_perf=0., cpu_perf=0., gpu_time=0., cpu_time=0.;
    real_Double_t    gpu_perf2, gpu_time2;
    double           error=0., error2=0., work[1];
    cuDoubleComplex *hA, *hR, *hR2, *hV, *hW;
    cuDoubleComplex *dV[MagmaMaxGPUs], *dW[MagmaMaxGPUs], *dA[MagmaMaxGPUs];

    /* Matrix size */
    magma_int_t n, size, lda, ldda;
    const int MAXTESTS = 10;
    // sizes are 1024*N - 32
    magma_int_t nsize[MAXTESTS] = { 992, 2016, 3040, 4064, 5088, 6112, 7136, 8160, 9184, 10208 };
    magma_int_t k       = 64;
    magma_int_t nb      = 64;
    magma_int_t nstream = 2;
    magma_int_t count   = 3;
    magma_int_t ngpu    = magma_num_gpus();
    
    magma_int_t info;
    magma_int_t ione     = 1;
    magma_int_t iseed[4] = {0,0,0,1};
        
    printf( "Usage: %s -N n -K k -nb nb -nstream nstream -ngpu ngpu -count count -c\n"
            "    -N can be repeated %d times.\n"
            "    -ngpu or setting $MAGMA_NUM_GPUS sets number of GPUs to use.\n"
            "    -c or setting $MAGMA_TESTINGS_CHECK runs LAPACK and checks result.\n",
            argv[0], MAXTESTS );

    int checkres = (getenv("MAGMA_TESTINGS_CHECK") != NULL);

    int ntest = 0;
    int nmax = 0;
    for( int i = 1; i < argc; i++ ) {
        if ( strcmp("-N", argv[i]) == 0 and i+1 < argc ) {
            ensure( ntest < MAXTESTS, "error: -N repeated more than maximum %d tests\n", MAXTESTS );
            nsize[ntest] = atoi( argv[++i] );
            ensure( nsize[ntest] > 0, "error: -N %s is invalid; must be > 0.\n", argv[i] );
            nmax = max( nmax, nsize[ntest] );
            ntest++;
        }
        else if ( strcmp("-K", argv[i]) == 0 and i+1 < argc ) {
            k = atoi( argv[++i] );
            ensure( k > 0, "error: -K %s is invalid; must be > 0.\n", argv[i] );
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
            ensure( nstream > 0 and nstream <= 10,
                    "error: -nstream %s is invalid; must be > 0 and <= 10.\n", argv[i] );
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
        nmax = nsize[ntest-1];
    }
    assert( nmax > 0 );
    
    // allocate memory for largest problem
    n = nmax;
    lda    = n;
    ldda   = ((n + 31)/32)*32;

    TESTING_MALLOC( hA,  cuDoubleComplex, lda*n );
    TESTING_MALLOC( hR,  cuDoubleComplex, lda*n );
    TESTING_MALLOC( hR2, cuDoubleComplex, lda*n );
    TESTING_MALLOC( hV,  cuDoubleComplex, lda*k*2 );
    //TESTING_MALLOC( hW,  cuDoubleComplex, lda*k );
    
    cudaStream_t streams[MagmaMaxGPUs][10];
    
    for( int d = 0; d < ngpu; ++d ) {
        magma_int_t nlocal = ((n / k) / ngpu + 1) * k;
        cudaSetDevice( d );
        TESTING_DEVALLOC( dA[d], cuDoubleComplex, ldda*nlocal );
        TESTING_DEVALLOC( dV[d], cuDoubleComplex, ldda*k*2      );
        //TESTING_DEVALLOC( dW[d], cuDoubleComplex, ldda*k      );
        for( int i = 0; i < nstream; ++i ) {
            cudaStreamCreate( &streams[d][i] );
        }
    }
    
    printf("\n");
    printf( "k %d, nb %d, ngpu %d, nstream %d\n", (int) k, (int) nb, (int) ngpu, (int) nstream );
    printf("    n    CPU GFlop/s (sec)   GPU GFlop/s (sec)   Ichi's code (sec)   ||R|| / ||A||\n");
    printf("=========================================================================\n");
    for( int i = 0; i < ntest; ++i ) {
    for( int j = 0; j < count; ++j ) {
        n = nsize[i];
        assert( n > 0 and k > 0 );
        
        lda  = n;
        ldda = ((n + 31)/32)*32;
        gflops = FLOPS_ZHER2K( (double)k, (double)n ) / 1e9;

        size = lda*n;
        lapackf77_zlarnv( &ione, iseed, &size, hA );
        size = lda*k*2;
        lapackf77_zlarnv( &ione, iseed, &size, hV );
        hW = hV + lda*k;
        //lapackf77_zlarnv( &ione, iseed, &size, hW );
        
        /* ====================================================================
           Performs operation using MAGMA
           =================================================================== */
        magmablas_zsetmatrix_1D_bcyclic( n, n, hA, lda, dA, ldda, ngpu, nb );
        for( int d = 0; d < ngpu; ++d ) {
            cudaSetDevice( d );
            dW[d] = dV[d] + ldda*k;
            magma_zsetmatrix( n, k, hV, lda, dV[d], ldda );
            magma_zsetmatrix( n, k, hW, lda, dW[d], ldda );
        }
        
        cudaDeviceSynchronize();
        gpu_time = magma_wtime();
        magmablas_zher2k_mgpu2(
            MagmaLower, MagmaNoTrans, n, k,
            c_neg_one, dV, ldda,
                       dW, ldda,
            d_one,     dA, ldda, 0,
            ngpu, nb, streams, nstream );
        cudaDeviceSynchronize();
        gpu_time = magma_wtime() - gpu_time;
        
        // Get dA back to the CPU to compare with the CPU result.
        magmablas_zgetmatrix_1D_bcyclic( n, n, dA, ldda, hR, lda, ngpu, nb );
                
        gpu_perf = gflops / gpu_time;
        
        /* ====================================================================
           Performs operation using MAGMA (Ichi's code)
           =================================================================== */
#if 0
        magmablas_zsetmatrix_1D_bcyclic( n, n, hA, lda, dA, ldda, ngpu, nb );
        for( int d = 0; d < ngpu; ++d ) {
            cudaSetDevice( d );
            magma_zsetmatrix( n, k, hV, lda, dV[d], ldda );
            magma_zsetmatrix( n, k, hW, lda, dW[d], ldda );
        }
        
        cudaDeviceSynchronize();
        gpu_time2 = magma_wtime();
        magma_zher2k_mgpu(
            ngpu, MagmaLower, MagmaNoTrans, nb, n, k,
            c_neg_one, dV, ldda,
                     //dW, ldda,
            d_one,     dA, ldda, 0,
            nstream, streams );
        cudaDeviceSynchronize();
        gpu_time2 = magma_wtime() - gpu_time2;
        
        // Get dA back to the CPU to compare with the CPU result.
        magmablas_zgetmatrix_1D_bcyclic( n, n, dA, ldda, hR2, lda, ngpu, nb );
                
        gpu_perf2 = gflops / gpu_time2;
#endif

        /* =====================================================================
           Performs operation using LAPACK
           =================================================================== */
        if ( checkres ) {
            // store ||A||
            error = lapackf77_zlange("f", &n, &n, hA, &lda, work );
            error2 = error;
            
            //printf( "A=" ); magma_zprint( n, n, hA, lda );
            //printf( "V=" ); magma_zprint( n, k, hV, lda );
            //printf( "W=" ); magma_zprint( n, k, hW, lda );
            
            cpu_time = magma_wtime();
            blasf77_zher2k( "Lower", "NoTrans", &n, &k,
                            &c_neg_one, hV, &lda,
                                        hW, &lda,
                            &d_one,     hA, &lda );
            cpu_time = magma_wtime() - cpu_time;
            
            //printf( "Ahat ="   );  magma_zprint( n, n, hA,  lda );
            //printf( "dAhat ="  );  magma_zprint( n, n, hR,  lda );
            //printf( "dAhat2 =" );  magma_zprint( n, n, hR2, lda );
            
            cpu_perf = gflops / cpu_time;
    
            // compute relative error ||R||/||A||, where R := A_magma - A_lapack = R - A
            size = lda*n;
            blasf77_zaxpy( &size, &c_neg_one, hA, &ione, hR, &ione );
            error = lapackf77_zlanhe("fro", "Lower", &n, hR, &lda, work) / error;
            
            blasf77_zaxpy( &size, &c_neg_one, hA, &ione, hR2, &ione );
            error2 = lapackf77_zlanhe("fro", "Lower", &n, hR2, &lda, work) / error2;
            
            printf( "%5d   %7.1f (%7.4f)   %7.1f (%7.4f)   %7.1f (%7.4f)   %8.2e   %8.2e\n",
                    (int) n, cpu_perf, cpu_time, gpu_perf, gpu_time, gpu_perf2, gpu_time2, error, error2 );
        }
        else {
            printf( "%5d     ---   (  ---  )   %7.1f (%7.4f)   %7.1f (%7.4f)     ---        ---\n",
                    (int) n, /*cpu_perf, cpu_time,*/ gpu_perf, gpu_time, gpu_perf2, gpu_time2 /*, error, error2*/ );
        }
    }}
    
    /* Memory clean up */
    for( int d = 0; d < ngpu; ++d ) {
        cudaSetDevice( d );
        TESTING_DEVFREE( dA[d] );
        TESTING_DEVFREE( dV[d] );
        //TESTING_DEVFREE( dW[d] );
    }
    
    TESTING_FREE( hA );
    TESTING_FREE( hR );
    TESTING_FREE( hV );
    //TESTING_FREE( hW );
    
    /* Shutdown */
    TESTING_CUDA_FINALIZE();
    return 0;
}
