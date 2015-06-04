/*
    -- MAGMA (version 1.2.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       June 2012

       @generated c Thu Jun 28 12:31:52 2012
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
magmablas_csymmetrize( char uplo, int m, cuFloatComplex *A, int lda );

int main( int argc, char** argv) 
{
    #define hA(i,j) (hA + (i) + (j)*lda)
    
    TESTING_CUDA_INIT();

    cuFloatComplex c_zero = MAGMA_C_ZERO;
    cuFloatComplex c_one  = MAGMA_C_ONE;
    
    cuFloatComplex *hA, *hR, *dA;
    real_Double_t   gpu_time, gpu_perf;

    int ione     = 1;
    int ISEED[4] = {0, 0, 0, 1};
    
    int nsize[] = { 32, 64, 96, 128, 100, 200 };
    int ntest = sizeof(nsize) / sizeof(int);
    int n   = nsize[ntest-1];
    int lda = ((n + 31)/32)*32;
    
    TESTING_MALLOC   ( hA, cuFloatComplex, lda*n );
    TESTING_MALLOC   ( hR, cuFloatComplex, lda*n );
    TESTING_DEVALLOC ( dA, cuFloatComplex, lda*n );
    
    for( int t = 0; t < ntest; ++t ) {
        n = nsize[t];
        lda = ((n + 31)/32)*32;
        
        // initialize matrices; entries are (i.j) for A
        float nf = 1000.;
        for( int i = 0; i < n; ++i ) {
            for( int j = 0; j < n; ++j ) {
                *hA(i,j) = MAGMA_C_MAKE( i + j/nf, 0. );
            }
        }
        printf( "A%d = ", n );
        magma_cprint( n, n, hA, lda );
        
        magma_csetmatrix( n, n, hA, lda, dA, lda );
        magmablas_csymmetrize( MagmaLower, n, dA, lda );
        magma_cgetmatrix( n, n, dA, lda, hR, lda );
        printf( "L%d = ", n );
        magma_cprint( n, n, hR, lda );
        
        magma_csetmatrix( n, n, hA, lda, dA, lda );
        magmablas_csymmetrize( MagmaUpper, n, dA, lda );
        magma_cgetmatrix( n, n, dA, lda, hR, lda );
        printf( "U%d = ", n );
        magma_cprint( n, n, hR, lda );
    }
    
    TESTING_FREE( hA );
    TESTING_FREE( hR );
    TESTING_DEVFREE( dA );
    
    /* Shutdown */
    TESTING_CUDA_FINALIZE();
    return 0;
}
