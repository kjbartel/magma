/*
 *  -- MAGMA (version 1.3.0) --
 *     Univ. of Tennessee, Knoxville
 *     Univ. of California, Berkeley
 *     Univ. of Colorado, Denver
 *     November 2012
 *
 *  @generated c Wed Nov 14 22:54:09 2012
 *
 **/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <cuda_runtime_api.h>
#include <cublas.h>
#include <cblas.h>

#include "flops.h"
#include "magma.h"
#include "magmablas.h"
#include "magma_lapack.h"
#include "testings.h"

#define PRECISION_c
#if defined(PRECISION_z) || defined(PRECISION_c)
#define FLOPS(n) ( 6. * FMULS_SYMV(n) + 2. * FADDS_SYMV(n))
#else
#define FLOPS(n) (      FMULS_SYMV(n) +      FADDS_SYMV(n))
#endif

#if (GPUSHMEM == 200)
extern "C" magma_int_t
magmablas_chemv2( char uplo, magma_int_t n,
                      cuFloatComplex alpha,
                      cuFloatComplex *A, magma_int_t lda,
                      cuFloatComplex *X, magma_int_t incx,
                      cuFloatComplex beta,
                      cuFloatComplex *Y, magma_int_t incy,
                      cuFloatComplex *work, magma_int_t lwork);
#endif

int main(int argc, char **argv)
{        
    TESTING_CUDA_INIT();

    magma_timestr_t  start, end;
    float      flops, magma_perf, cuda_perf, error, work[1];
    magma_int_t ione     = 1;
    magma_int_t ISEED[4] = {0,0,0,1};
    cuFloatComplex c_neg_one = MAGMA_C_NEG_ONE;

    FILE        *fp ; 
    magma_int_t N, m, i, lda, LDA;
    magma_int_t matsize;
    magma_int_t vecsize;
    magma_int_t istart = 64;
    magma_int_t incx = 1;
    char        uplo = MagmaLower;
    cuFloatComplex alpha = MAGMA_C_MAKE(1., 0.); // MAGMA_C_MAKE(  1.5, -2.3 );
    cuFloatComplex beta  = MAGMA_C_MAKE(0., 0.); // MAGMA_C_MAKE( -0.6,  0.8 );
    cuFloatComplex *A, *X, *Y, *Ycublas, *Ymagma;
    cuFloatComplex *dA, *dX, *dY;
    int nb = 64;

     cuFloatComplex *C_work;
     cuFloatComplex *dC_work;
    
    fp = fopen ("results_chemv.txt", "w") ;
    if( fp == NULL ){ printf("Couldn't open output file\n"); exit(1);}

    printf("HEMV cuFloatComplex Precision\n\n"
           "Usage\n\t\t testing_chemv U|L N\n\n");

#ifdef PRECISION_z
    N = 8*1024+64;
#else
    N = 15*1024+64;
#endif 

    if( argc > 1 ) {
      uplo = argv[1][0];
    }
    if( argc > 2 ) {
      istart = N = atoi( argv[2] );
    }
    LDA = ((N+31)/32)*32;
    matsize = N*LDA;
    vecsize = N*incx;

    TESTING_MALLOC( A, cuFloatComplex, matsize );
    TESTING_MALLOC( X, cuFloatComplex, vecsize );
    TESTING_MALLOC( Y, cuFloatComplex, vecsize );
    TESTING_MALLOC( Ycublas, cuFloatComplex, vecsize );
    TESTING_MALLOC( Ymagma,  cuFloatComplex, vecsize );

    TESTING_DEVALLOC( dA, cuFloatComplex, matsize );
    TESTING_DEVALLOC( dX, cuFloatComplex, vecsize );
    TESTING_DEVALLOC( dY, cuFloatComplex, vecsize );


    int blocks    = N / nb + (N % nb != 0);
    int workspace = LDA * (blocks + 1);
    TESTING_MALLOC(    C_work, cuFloatComplex, workspace );
    TESTING_DEVALLOC( dC_work, cuFloatComplex, workspace );

    /* Initialize the matrix */
    lapackf77_clarnv( &ione, ISEED, &matsize, A );
    /* Make A hermitian */
    { 
        magma_int_t i, j;
        for(i=0; i<N; i++) {
            A[i*LDA+i] = MAGMA_C_MAKE( MAGMA_C_REAL(A[i*LDA+i]), 0. );
            for(j=0; j<i; j++)
                A[i*LDA+j] = cuConjf(A[j*LDA+i]);
        }
    }
        
    printf( "   n   CUBLAS,Gflop/s   MAGMABLAS,Gflop/s      \"error\"\n" 
            "==============================================================\n");
    fprintf(fp, "   n   CUBLAS,Gflop/s   MAGMABLAS,Gflop/s      \"error\"\n" 
            "==============================================================\n");
    
    for( i = istart; i<N+1; i = (int)((i+1)*1.1) )
    {
        m = i;
        lda = ((m+31)/32)*32;
        flops = FLOPS( (float)m ) / 1e6;

        printf(      "%5d ", (int) m );
        fprintf( fp, "%5d ", (int) m );

        vecsize = m * incx;
        lapackf77_clarnv( &ione, ISEED, &vecsize, X );
        lapackf77_clarnv( &ione, ISEED, &vecsize, Y );

        /* =====================================================================
           Performs operation using CUDA-BLAS
           =================================================================== */
        magma_csetmatrix( m, m, A, LDA, dA, lda );
        magma_csetvector( m, X, incx, dX, incx );
        magma_csetvector( m, Y, incx, dY, incx );


            blocks    = m / nb + (m % nb != 0);
            magma_csetmatrix( lda, blocks, C_work, LDA, dC_work, lda );
        start = get_current_time();
        cublasChemv( uplo, m, alpha, dA, lda, dX, incx, beta, dY, incx );
        end = get_current_time();

        magma_cgetvector( m, dY, incx, Ycublas, incx );
        
        
        cuda_perf = flops / GetTimerValue(start,end);
        printf(     "%11.2f", cuda_perf );
        fprintf(fp, "%11.2f", cuda_perf );
        
        magma_csetvector( m, Y, incx, dY, incx );
        magmablas_chemv( uplo, m, alpha, dA, lda, dX, incx, beta, dY, incx );
        magma_csetvector( m, Y, incx, dY, incx );
        
        
        start = get_current_time();
#if (GPUSHMEM == 200)
        magmablas_chemv2( uplo, m, alpha, dA, lda, dX, incx, beta, dY, incx, dC_work, lda * (blocks + 1));
#else
        magmablas_chemv( uplo, m, alpha, dA, lda, dX, incx, beta, dY, incx );
#endif
        end = get_current_time();
        
        magma_cgetvector( m, dY, incx, Ymagma, incx );

        magma_perf = flops / GetTimerValue(start,end);
        printf(     "%11.2f", magma_perf );
        fprintf(fp, "%11.2f", magma_perf );

        /* =====================================================================
           Computing the Difference Cublas VS Magma
           =================================================================== */
        
        blasf77_caxpy( &m, &c_neg_one, Ymagma, &incx, Ycublas, &incx);
        error = lapackf77_clange( "M", &m, &ione, Ycublas, &m, work );
            
#if 0
        printf(      "\t\t %8.6e", error / m );
        fprintf( fp, "\t\t %8.6e", error / m );

        /*
         * Extra check with cblas vs magma
         */
        cblas_ccopy( m, Y, incx, Ycublas, incx );
        cblas_chemv( CblasColMajor, CblasLower, m, 
                     CBLAS_SADDR(alpha), A, LDA, X, incx, 
                     CBLAS_SADDR(beta), Ycublas, incx );
 
        blasf77_caxpy( &m, &c_neg_one, Ymagma, &incx, Ycublas, &incx);
        error = lapackf77_clange( "M", &m, &ione, Ycublas, &m, work );
#endif

        printf(      "\t\t %8.6e\n", error / m );
        fprintf( fp, "\t\t %8.6e\n", error / m );
    }
    
    fclose( fp ) ; 

    /* Free Memory */
    TESTING_FREE( A );
    TESTING_FREE( X );
    TESTING_FREE( Y );
    TESTING_FREE( Ycublas );
    TESTING_FREE( Ymagma );

    TESTING_DEVFREE( dA );
    TESTING_DEVFREE( dX );
    TESTING_DEVFREE( dY );
    TESTING_FREE( C_work );
    TESTING_DEVFREE( dC_work );

    /* Free device */
    TESTING_CUDA_FINALIZE();
    return 0;
}        
