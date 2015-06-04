/*
    -- MAGMA (version 1.3.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2012

       @generated d Wed Nov 14 22:54:11 2012
*/

// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <cuda_runtime_api.h>
#include <cublas.h>
#include <cblas.h>

// includes, project
#include "flops.h"
#include "magma.h"
#include "magma_lapack.h"
#include "testings.h"

/* ////////////////////////////////////////////////////////////////////////////
   -- Testing dlacpy
*/
#define PRECISION_d

int main( int argc, char** argv) 
{
    TESTING_CUDA_INIT();

    double c_zero = MAGMA_D_ZERO;
    double c_one  = MAGMA_D_ONE;
    
    magma_timestr_t  start, end;
    double *hA, *hB, *hR, *dA, *dB;
    double           gpu_time, gpu_perf;

    magma_int_t ione     = 1;
    magma_int_t ISEED[4] = {0, 0, 0, 1};
    
    // groups of tests are:
    // whole matrix, sub-matrix, around k*64 rows, around k*64 cols,
    // zero rows, one row, zero cols, one col
    magma_int_t TESTS_I1[] = {     0,  100,     63,   64,   64,   64,   65,     10,   10,   10,   10,   10,      4,   4,   4,      4,   4,   4,     64,  64,  64,     64,  64,  64 };
    magma_int_t TESTS_I2[] = {  1000,  500,    511,  511,  512,  513,  513,    900,  900,  900,  900,  900,      4,   4,   4,      5,   5,   5,    127, 128, 129,    255, 256, 257 };
    magma_int_t TESTS_J1[] = {     0,   50,     10,   10,   10,   10,   10,     63,   64,   64,   64,   65,     64,  64,  64,     64,  64,  64,      4,   4,   4,      4,   4,   4 };
    magma_int_t TESTS_J2[] = {  1000,  400,    900,  900,  900,  900,  900,    511,  511,  512,  513,  513,    127, 128, 129,    255, 256, 257,      4,   4,   4,      5,   5,   5 };
    int ntest = sizeof(TESTS_J2) / sizeof(magma_int_t);
    
    magma_int_t n   = 1000;
    magma_int_t lda = n;
    
    TESTING_MALLOC   ( hA, double, lda*n );
    TESTING_MALLOC   ( hB, double, lda*n );
    TESTING_MALLOC   ( hR, double, lda*n );
    TESTING_DEVALLOC ( dA, double, lda*n );
    TESTING_DEVALLOC ( dB, double, lda*n );
    
    // initialize matrices; entries are (i.j) for A and (800 + i.j) for B.
    double nf = n;
    for( int i = 0; i < n; ++i ) {
        for( int j = 0; j < n; ++j ) {
            hA[i + j*lda] = MAGMA_D_MAKE( i + j/nf,       0. );
            hB[i + j*lda] = MAGMA_D_MAKE( i + j/nf + 800, 0. );
        }
    }
    
    printf( "\nNote: ranges use Python notation,\n"
            "i.e., A[i:j] is A[ i, i+1, ..., j-1 ], excluding A[j].\n\n" );
    for( int t = 0; t < ntest; ++t ) {
        magma_dsetmatrix( n, n, hA, lda, dA, lda );
        magma_dsetmatrix( n, n, hB, lda, dB, lda );
        
        // copy submatrix
        int i1 = TESTS_I1[ t ];
        int i2 = TESTS_I2[ t ];
        int j1 = TESTS_J1[ t ];
        int j2 = TESTS_J2[ t ];
        magmablas_dlacpy( 'F', i2-i1, j2-j1,
                          &dA[i1 + j1*lda], lda,
                          &dB[i1 + j1*lda], lda );
        
        // verify result
        int bad_copies = 0;
        int overwrites = 0;
        magma_dgetmatrix( n, n, dB, lda, hR, lda );
        
        for( int j = 0; j < n; ++j ) {
            for( int i = 0; i < n; ++i ) {
                if ( i1 <= i && i < i2 && j1 <= j && j < j2 ) {
                    if ( ! MAGMA_D_EQUAL( hR[i + j*lda], hA[i + j*lda] )) {
                        bad_copies += 1;
                        printf( "Copy failed at B[%d,%d], expected %.4f, got %.4f\n",
                                i, j, MAGMA_D_REAL( hA[i + j*lda] ),
                                      MAGMA_D_REAL( hR[i + j*lda] ));
                    }
                }
                else {
                    if ( ! MAGMA_D_EQUAL( hR[i + j*lda], hB[i + j*lda] )) {
                        overwrites += 1;
                        printf( "Overwrote at B[%d,%d], expected %.4f, got %.4f\n",
                                i, j, MAGMA_D_REAL( hA[i + j*lda] ),
                                      MAGMA_D_REAL( hR[i + j*lda] ));
                    }
                }
            }
        }
        printf( "B(%4d:%4d, %4d:%4d) = A(%4d:%4d, %4d:%4d) ",
                i1, i2, j1, j2,
                i1, i2, j1, j2 );
        if ( bad_copies > 0 || overwrites > 0 ) {
            printf( "failed, %d bad copies, %d overwrites\n", bad_copies, overwrites );
        }
        else {
            printf( "passed\n" );
        }
    }
    
    TESTING_FREE( hA );
    TESTING_FREE( hB );
    TESTING_FREE( hR );
    TESTING_DEVFREE( dA );
    TESTING_DEVFREE( dB );
    
    // --------------------------------------------------
    // speed tests
    magma_int_t SIZE[] = {
        1024, 1280, 1536, 1792, 2048, 2304, 2560, 2816, 3072, 3328, 3584, 3840,
        4096, 4352, 4608, 4864, 5120, 5376, 5632, 5888, 6144, 6400, 6656, 6912,
        7168, 7424, 7680, 7936, 8192, 8448, 8704, 8960, 9216, 9472, 9728, 9984
    };
    magma_int_t nsize = sizeof(SIZE) / sizeof(magma_int_t);
    
    printf("\n  N      GPU MB/s (sec)\n");
    printf("========================================\n");
    for( int t = 0; t < nsize; ++t ) {
        n = SIZE[ t ];
        lda = n;
        TESTING_MALLOC   ( hA, double, lda*n );
        TESTING_MALLOC   ( hB, double, lda*n );
        TESTING_DEVALLOC ( dA, double, lda*n );
        TESTING_DEVALLOC ( dB, double, lda*n );
        
        // initialize matrices
        magma_int_t n2 = lda*n;
        lapackf77_dlarnv( &ione, ISEED, &n2, hA );
        lapackf77_dlaset( "F", &n, &n, &c_zero, &c_zero, hB, &lda );
        magma_dsetmatrix( n, n, hA, lda, dA, lda );
        magmablas_dlaset( 'F', n, n, /*c_zero,*/ dB, lda );
        
        start = get_current_time();
        magmablas_dlacpy( 'F', n, n, dA, lda, dB, lda );
        end = get_current_time();
        
        // verify copy
        magma_dgetmatrix( n, n, dB, lda, hB, lda );
        for( int j = 0; j < n; ++j ) {
            for( int i = 0; i < n; ++i ) {
                if ( ! MAGMA_D_EQUAL( hA[i + j*lda], hB[i + j*lda] )) {
                    printf( "Copy failed at B[%d,%d], expected %.4f, got %.4f\n",
                            i, j, MAGMA_D_REAL( hA[i + j*lda] ),
                                  MAGMA_D_REAL( hB[i + j*lda] ));
                    exit(1);
                }
            }
        }
        
        gpu_time = GetTimerValue( start, end ) * 1e-3;
        gpu_perf = n*n*sizeof(double) / 1024. / 1024. / gpu_time;
        printf( "%5d    %6.2f (%8.6f)\n", (int) n, gpu_perf, gpu_time );
        
        TESTING_FREE   ( hA );
        TESTING_FREE   ( hB );
        TESTING_DEVFREE( dA );
        TESTING_DEVFREE( dB );
    }
    
    /* Shutdown */
    TESTING_CUDA_FINALIZE();
    return EXIT_SUCCESS;
}
