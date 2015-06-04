/*
    -- MAGMA (version 1.2.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       June 2012

       @precisions normal z -> c d s
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

extern "C" void
magmablas_zsymmetrize( char uplo, int m, cuDoubleComplex *A, int lda );

int main( int argc, char** argv) 
{
    #define hA(i,j) (hA + (i) + (j)*lda)
    
    TESTING_CUDA_INIT();

    cuDoubleComplex c_zero = MAGMA_Z_ZERO;
    cuDoubleComplex c_one  = MAGMA_Z_ONE;
    
    cuDoubleComplex *hA, *hR, *dA;
    real_Double_t   gpu_time, gpu_perf;

    int ione     = 1;
    int ISEED[4] = {0, 0, 0, 1};
    
    int nsize[] = { 32, 64, 96, 128, 100, 200 };
    int ntest = sizeof(nsize) / sizeof(int);
    int n   = nsize[ntest-1];
    int lda = ((n + 31)/32)*32;
    
    TESTING_MALLOC   ( hA, cuDoubleComplex, lda*n );
    TESTING_MALLOC   ( hR, cuDoubleComplex, lda*n );
    TESTING_DEVALLOC ( dA, cuDoubleComplex, lda*n );
    
    for( int t = 0; t < ntest; ++t ) {
        n = nsize[t];
        lda = ((n + 31)/32)*32;
        
        // initialize matrices; entries are (i.j) for A
        double nf = 1000.;
        for( int i = 0; i < n; ++i ) {
            for( int j = 0; j < n; ++j ) {
                *hA(i,j) = MAGMA_Z_MAKE( i + j/nf, 0. );
            }
        }
        printf( "A%d = ", n );
        magma_zprint( n, n, hA, lda );
        
        magma_zsetmatrix( n, n, hA, lda, dA, lda );
        magmablas_zsymmetrize( MagmaLower, n, dA, lda );
        magma_zgetmatrix( n, n, dA, lda, hR, lda );
        printf( "L%d = ", n );
        magma_zprint( n, n, hR, lda );
        
        magma_zsetmatrix( n, n, hA, lda, dA, lda );
        magmablas_zsymmetrize( MagmaUpper, n, dA, lda );
        magma_zgetmatrix( n, n, dA, lda, hR, lda );
        printf( "U%d = ", n );
        magma_zprint( n, n, hR, lda );
    }
    
    TESTING_FREE( hA );
    TESTING_FREE( hR );
    TESTING_DEVFREE( dA );
    
    /* Shutdown */
    TESTING_CUDA_FINALIZE();
    return 0;
}
