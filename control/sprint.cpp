/*
    -- MAGMA (version 1.2.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       June 2012

       @author Mark Gates
       @generated s Thu Jun 28 12:30:05 2012

*/
#include "common_magma.h"

#define A(i,j) (A + i + j*lda)

// -------------------------
// Prints a matrix that is on the CPU host.
extern "C"
void magma_sprint( magma_int_t m, magma_int_t n, float *A, magma_int_t lda )
{
    if ( magma_is_devptr( A ) == 1 ) {
        fprintf( stderr, "ERROR: sprint called with device pointer.\n" );
        exit(1);
    }
    
    float c_zero = MAGMA_S_ZERO;
    
    printf( "[\n" );
    for( int i = 0; i < m; ++i ) {
        for( int j = 0; j < n; ++j ) {
            if ( MAGMA_S_EQUAL( *A(i,j), c_zero )) {
                printf( "   0.    " );
            }
            else {
                printf( " %8.4f", MAGMA_S_REAL( *A(i,j) ));
            }
        }
        printf( "\n" );
    }
    printf( "];\n" );
}

// -------------------------
// Prints a matrix that is on the GPU device.
// Internally allocates memory on host, copies it to the host, prints it,
// and de-allocates host memory.
extern "C"
void magma_sprint_gpu( magma_int_t m, magma_int_t n, float *dA, magma_int_t ldda )
{
    if ( magma_is_devptr( dA ) == 0 ) {
        fprintf( stderr, "ERROR: sprint_gpu called with host pointer.\n" );
        exit(1);
    }
    
    int lda = m;
    float* A = (float*) malloc( lda*n*sizeof(float) );
    cublasGetMatrix( m, n, sizeof(float), dA, ldda, A, lda );
    
    magma_sprint( m, n, A, lda );
    
    free( A );
}
