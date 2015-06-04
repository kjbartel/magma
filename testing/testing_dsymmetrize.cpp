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

int main( int argc, char** argv) 
{
    #define hA(i,j) (hA + (i) + (j)*lda)
    
    TESTING_CUDA_INIT();

    double c_zero = MAGMA_D_ZERO;
    double c_one  = MAGMA_D_ONE;
    
    double *hA, *hR, *dA;
    //real_Double_t   gpu_time, gpu_perf;

    //int ione     = 1;
    //int ISEED[4] = {0, 0, 0, 1};
    
    int nsize[] = { 32, 64, 96, 256, 100, 200, 512 };
    int ntest = sizeof(nsize) / sizeof(int);
    int n   = nsize[ntest-1];
    int lda = ((n + 31)/32)*32;
    int ntile, nb;
    
    TESTING_MALLOC   ( hA, double, lda*n );
    TESTING_MALLOC   ( hR, double, lda*n );
    TESTING_DEVALLOC ( dA, double, lda*n );
    
    for( int t = 0; t < ntest; ++t ) {
        n = nsize[t];
        lda = ((n + 31)/32)*32;
        
        // initialize matrices; entries are (i.j) for A
        double nf = 100.;
        for( int j = 0; j < n; ++j ) {
            // upper
            for( int i = 0; i < j; ++i ) {
                *hA(i,j) = MAGMA_D_MAKE( (i + j/nf)/nf, 0. );
            }
            // lower
            for( int i = j; i < n; ++i ) {
                *hA(i,j) = MAGMA_D_MAKE( i + j/nf, 0. );
            }
        }
        printf( "A%d = ", n );
        magma_dprint( n, n, hA, lda );
        
        magma_dsetmatrix( n, n, hA, lda, dA, lda );
        magmablas_dsymmetrize( MagmaLower, n, dA, lda );
        magma_dgetmatrix( n, n, dA, lda, hR, lda );
        printf( "L%d = ", n );
        magma_dprint( n, n, hR, lda );
        
        magma_dsetmatrix( n, n, hA, lda, dA, lda );
        magmablas_dsymmetrize( MagmaUpper, n, dA, lda );
        magma_dgetmatrix( n, n, dA, lda, hR, lda );
        printf( "U%d = ", n );
        magma_dprint( n, n, hR, lda );
        
        // -----
        //lapackf77_dlaset( "u", &n, &n, &c_zero, &c_one, hA, &lda );
        
        nb = 64;
        ntile = n / nb;
        magma_dsetmatrix( n, n, hA, lda, dA, lda );
        magmablas_dsymmetrize_tiles( MagmaLower, nb, dA, lda, ntile, nb, nb );
        magma_dgetmatrix( n, n, dA, lda, hR, lda );
        printf( "L%d_%d = ", n, nb );
        magma_dprint( n, n, hR, lda );
        
        nb = 32;
        ntile = n / nb;
        magma_dsetmatrix( n, n, hA, lda, dA, lda );
        magmablas_dsymmetrize_tiles( MagmaLower, nb, dA, lda, ntile, nb, nb );
        magma_dgetmatrix( n, n, dA, lda, hR, lda );
        printf( "L%d_%d = ", n, nb );
        magma_dprint( n, n, hR, lda );
        
        ntile = (n - nb < 0 ? 0 : (n - nb) / (2*nb) + 1);
        magma_dsetmatrix( n, n, hA, lda, dA, lda );
        magmablas_dsymmetrize_tiles( MagmaLower, nb, dA, lda, ntile, 2*nb, nb );
        magma_dgetmatrix( n, n, dA, lda, hR, lda );
        printf( "L%d_%d_2m = ", n, nb );
        magma_dprint( n, n, hR, lda );
        
        nb = 25;
        ntile = n / nb;
        magma_dsetmatrix( n, n, hA, lda, dA, lda );
        magmablas_dsymmetrize_tiles( MagmaLower, nb, dA, lda, ntile, nb, nb );
        magma_dgetmatrix( n, n, dA, lda, hR, lda );
        printf( "L%d_%d = ", n, nb );
        magma_dprint( n, n, hR, lda );
        
        nb = 25;
        ntile = (n - nb < 0 ? 0 : (n - nb) / (3*nb) + 1);
        magma_dsetmatrix( n, n, hA, lda, dA, lda );
        magmablas_dsymmetrize_tiles( MagmaLower, nb, dA, lda, ntile, nb, 3*nb );
        magma_dgetmatrix( n, n, dA, lda, hR, lda );
        printf( "L%d_%d_3n = ", n, nb );
        magma_dprint( n, n, hR, lda );
        
        nb = 100;
        ntile = n / nb;
        magma_dsetmatrix( n, n, hA, lda, dA, lda );
        magmablas_dsymmetrize_tiles( MagmaLower, nb, dA, lda, ntile, nb, nb );
        magmablas_dsymmetrize( MagmaLower, n%nb, &dA[ ntile*nb*(1+lda) ], lda );  // last partial block
        magma_dgetmatrix( n, n, dA, lda, hR, lda );
        printf( "L%d_%d = ", n, nb );
        magma_dprint( n, n, hR, lda );
        
        // -----
        nb = 64;
        ntile = n / nb;
        magma_dsetmatrix( n, n, hA, lda, dA, lda );
        magmablas_dsymmetrize_tiles( MagmaUpper, nb, dA, lda, ntile, nb, nb );
        magma_dgetmatrix( n, n, dA, lda, hR, lda );
        printf( "U%d_%d = ", n, nb );
        magma_dprint( n, n, hR, lda );
    }
    
    TESTING_FREE( hA );
    TESTING_FREE( hR );
    TESTING_DEVFREE( dA );
    
    /* Shutdown */
    TESTING_CUDA_FINALIZE();
    return 0;
}
