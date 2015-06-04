/*
    -- MAGMA (version 1.2.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       May 2012

       @author Mark Gates
       @precisions normal z -> s d c

*/
#include "common_magma.h"

// -------------------------
// Returns:
//  1 if A is a device pointer (definitely),
//  0 if A is a host   pointer (definitely or inferred from error),
// -1 if unknown.
// On 2.0 cards with unified addressing, CUDA can tell if this is a device pointer.
// For malloc'd host pointers, cudaPointerGetAttributes returns error.
static int is_devptr( void* A )
{
    cudaError_t err;
    cudaDeviceProp prop;
    cudaPointerAttributes attr;
    int dev;
    err = cudaGetDevice( &dev );
    if ( ! err ) {
        err = cudaGetDeviceProperties( &prop, dev );
        if ( ! err and prop.unifiedAddressing ) {
            err = cudaPointerGetAttributes( &attr, A );
            if ( ! err ) {
                // definitely know type
                return (attr.memoryType == cudaMemoryTypeDevice);
            }
            else if ( err == cudaErrorInvalidValue ) {
                // infer as host pointer
                return 0;
            }
        }
    }
    // unknown, e.g., device doesn't support unified addressing
    return -1;
}


#define A(i,j) (A + i + j*lda)

// -------------------------
// Prints a matrix that is on the CPU host.
extern "C"
void magma_zprint( int m, int n, cuDoubleComplex *A, int lda )
{
    if ( is_devptr( A ) == 1 ) {
        fprintf( stderr, "ERROR: zprint called with device pointer.\n" );
        exit(1);
    }
    
    cuDoubleComplex c_zero = MAGMA_Z_ZERO;
    
    printf( "[\n" );
    for( int i = 0; i < m; ++i ) {
        for( int j = 0; j < n; ++j ) {
            if ( MAGMA_Z_EQUAL( *A(i,j), c_zero )) {
                printf( "   0.    " );
            }
            else {
                printf( " %8.4f", MAGMA_Z_REAL( *A(i,j) ));
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
void magma_zprint_gpu( int m, int n, cuDoubleComplex *dA, int ldda )
{
    if ( is_devptr( dA ) == 0 ) {
        fprintf( stderr, "ERROR: zprint_gpu called with host pointer.\n" );
        exit(1);
    }
    
    int lda = m;
    cuDoubleComplex* A = (cuDoubleComplex*) malloc( lda*n*sizeof(cuDoubleComplex) );
    cublasGetMatrix( m, n, sizeof(cuDoubleComplex), dA, ldda, A, lda );
    
    magma_zprint( m, n, A, lda );
    
    free( A );
}
