/*
    -- MAGMA (version 1.5.0-beta1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date April 2014
  
       @generated from testing_zprint.cpp normal z -> s, Fri Apr 25 15:06:07 2014
       @author Mark Gates
*/
// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <cuda_runtime_api.h>
#include <cublas.h>

// includes, project
#include "flops.h"
#include "magma.h"
#include "magma_lapack.h"
#include "testings.h"


/* ////////////////////////////////////////////////////////////////////////////
   -- Testing sprint
*/
int main( int argc, char** argv)
{
    TESTING_INIT();

    float *hA, *dA;
    //magma_int_t ione     = 1;
    //magma_int_t ISEED[4] = {0,0,0,1};
    //magma_int_t size;
    magma_int_t M, N, lda, ldda, size;
    
    magma_opts opts;
    parse_opts( argc, argv, &opts );
    
    for( int itest = 0; itest < opts.ntest; ++itest ) {
        for( int iter = 0; iter < opts.niter; ++iter ) {
            M     = opts.msize[itest];
            N     = opts.nsize[itest];
            lda   = N;
            ldda  = ((M + 31)/32)*32;
            size  = lda*N;

            /* Allocate host memory for the matrix */
            TESTING_MALLOC_CPU( hA, float, lda *N );
            TESTING_MALLOC_DEV( dA, float, ldda*N );
        
            //lapackf77_slarnv( &ione, ISEED, &size, hA );
            for( int j = 0; j < N; ++j ) {
                for( int i = 0; i < M; ++i ) {
                    hA[i + j*lda] = MAGMA_S_MAKE( i + j*0.01, 0. );
                }
            }
            magma_ssetmatrix( M, N, hA, lda, dA, ldda );
            
            printf( "A=" );
            magma_sprint( M, N, hA, lda );
            
            printf( "dA=" );
            magma_sprint_gpu( M, N, dA, ldda );
            
            TESTING_FREE_CPU( hA );
            TESTING_FREE_DEV( dA );
        }
    }

    TESTING_FINALIZE();
    return 0;
}
