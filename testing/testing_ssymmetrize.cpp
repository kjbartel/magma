/*
    -- MAGMA (version 1.3.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2012

       @generated s Wed Nov 14 22:54:11 2012
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

int main( int argc, char** argv) 
{
    #define hA(i,j) (hA + (i) + (j)*lda)
    
    TESTING_CUDA_INIT();

    float c_zero = MAGMA_S_ZERO;
    float c_one  = MAGMA_S_ONE;
    
    float *hA, *hR, *dA;
    //real_Double_t   gpu_time, gpu_perf;

    //int ione     = 1;
    //int ISEED[4] = {0, 0, 0, 1};
    
    int nsize[] = { 32, 64, 96, 256, 100, 200, 512 };
    int ntest = sizeof(nsize) / sizeof(int);
    int n   = nsize[ntest-1];
    int lda = ((n + 31)/32)*32;
    int ntile, nb;
    
    TESTING_MALLOC   ( hA, float, lda*n );
    TESTING_MALLOC   ( hR, float, lda*n );
    TESTING_DEVALLOC ( dA, float, lda*n );
    
    for( int t = 0; t < ntest; ++t ) {
        n = nsize[t];
        lda = ((n + 31)/32)*32;
        
        // initialize matrices; entries are (i.j) for A
        float nf = 100.;
        for( int j = 0; j < n; ++j ) {
            // upper
            for( int i = 0; i < j; ++i ) {
                *hA(i,j) = MAGMA_S_MAKE( (i + j/nf)/nf, 0. );
            }
            // lower
            for( int i = j; i < n; ++i ) {
                *hA(i,j) = MAGMA_S_MAKE( i + j/nf, 0. );
            }
        }
        printf( "A%d = ", n );
        magma_sprint( n, n, hA, lda );
        
        magma_ssetmatrix( n, n, hA, lda, dA, lda );
        magmablas_ssymmetrize( MagmaLower, n, dA, lda );
        magma_sgetmatrix( n, n, dA, lda, hR, lda );
        printf( "L%d = ", n );
        magma_sprint( n, n, hR, lda );
        
        magma_ssetmatrix( n, n, hA, lda, dA, lda );
        magmablas_ssymmetrize( MagmaUpper, n, dA, lda );
        magma_sgetmatrix( n, n, dA, lda, hR, lda );
        printf( "U%d = ", n );
        magma_sprint( n, n, hR, lda );
        
        // -----
        //lapackf77_slaset( "u", &n, &n, &c_zero, &c_one, hA, &lda );
        
        nb = 64;
        ntile = n / nb;
        magma_ssetmatrix( n, n, hA, lda, dA, lda );
        magmablas_ssymmetrize_tiles( MagmaLower, nb, dA, lda, ntile, nb, nb );
        magma_sgetmatrix( n, n, dA, lda, hR, lda );
        printf( "L%d_%d = ", n, nb );
        magma_sprint( n, n, hR, lda );
        
        nb = 32;
        ntile = n / nb;
        magma_ssetmatrix( n, n, hA, lda, dA, lda );
        magmablas_ssymmetrize_tiles( MagmaLower, nb, dA, lda, ntile, nb, nb );
        magma_sgetmatrix( n, n, dA, lda, hR, lda );
        printf( "L%d_%d = ", n, nb );
        magma_sprint( n, n, hR, lda );
        
        ntile = (n - nb < 0 ? 0 : (n - nb) / (2*nb) + 1);
        magma_ssetmatrix( n, n, hA, lda, dA, lda );
        magmablas_ssymmetrize_tiles( MagmaLower, nb, dA, lda, ntile, 2*nb, nb );
        magma_sgetmatrix( n, n, dA, lda, hR, lda );
        printf( "L%d_%d_2m = ", n, nb );
        magma_sprint( n, n, hR, lda );
        
        nb = 25;
        ntile = n / nb;
        magma_ssetmatrix( n, n, hA, lda, dA, lda );
        magmablas_ssymmetrize_tiles( MagmaLower, nb, dA, lda, ntile, nb, nb );
        magma_sgetmatrix( n, n, dA, lda, hR, lda );
        printf( "L%d_%d = ", n, nb );
        magma_sprint( n, n, hR, lda );
        
        nb = 25;
        ntile = (n - nb < 0 ? 0 : (n - nb) / (3*nb) + 1);
        magma_ssetmatrix( n, n, hA, lda, dA, lda );
        magmablas_ssymmetrize_tiles( MagmaLower, nb, dA, lda, ntile, nb, 3*nb );
        magma_sgetmatrix( n, n, dA, lda, hR, lda );
        printf( "L%d_%d_3n = ", n, nb );
        magma_sprint( n, n, hR, lda );
        
        nb = 100;
        ntile = n / nb;
        magma_ssetmatrix( n, n, hA, lda, dA, lda );
        magmablas_ssymmetrize_tiles( MagmaLower, nb, dA, lda, ntile, nb, nb );
        magmablas_ssymmetrize( MagmaLower, n%nb, &dA[ ntile*nb*(1+lda) ], lda );  // last partial block
        magma_sgetmatrix( n, n, dA, lda, hR, lda );
        printf( "L%d_%d = ", n, nb );
        magma_sprint( n, n, hR, lda );
        
        // -----
        nb = 64;
        ntile = n / nb;
        magma_ssetmatrix( n, n, hA, lda, dA, lda );
        magmablas_ssymmetrize_tiles( MagmaUpper, nb, dA, lda, ntile, nb, nb );
        magma_sgetmatrix( n, n, dA, lda, hR, lda );
        printf( "U%d_%d = ", n, nb );
        magma_sprint( n, n, hR, lda );
    }
    
    TESTING_FREE( hA );
    TESTING_FREE( hR );
    TESTING_DEVFREE( dA );
    
    /* Shutdown */
    TESTING_CUDA_FINALIZE();
    return 0;
}
