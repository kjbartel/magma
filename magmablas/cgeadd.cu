 /*
    -- MAGMA (version 1.3.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2012

       @generated c Wed Nov 14 22:53:45 2012
       @author Mark Gates
*/
#include "common_magma.h"
#include <assert.h>

#define NB 64

/*
    Matrix is m x n, and is divided into block rows, each NB x n.
    Each block has NB threads.
    Each thread adds one row, iterating across all columns.
    The bottom block of rows may be partially outside the matrix;
    if so, rows outside the matrix (i >= m) are disabled.
    
    TODO. Block in both directions, for large matrices.
    E.g., each block does 64x64 tile, instead of 64xN tile.
*/
__global__ void
cgeadd_kernel(
    int m, int n,
    cuFloatComplex alpha,
    const cuFloatComplex *dA, int ldda,
    cuFloatComplex       *dB, int lddb )
{
    // dA and dB iterate across row i
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    if ( i < m ) {
        dA += i;
        dB += i;
        const cuFloatComplex *dAend = dA + n*ldda;
        while( dA < dAend ) {
            *dB = alpha*(*dA) + (*dB);
            dA += ldda;
            dB += lddb;
        }
    }
}


extern "C" void
magmablas_cgeadd(
    magma_int_t m, magma_int_t n,
    cuFloatComplex alpha,
    const cuFloatComplex *dA, magma_int_t ldda,
    cuFloatComplex       *dB, magma_int_t lddb )
{
/*
    Purpose
    =======
    
    ZGEADD adds two matrices, B = alpha*A + B.
    
    Arguments
    =========
    
    M       (input) INTEGER
            The number of rows of the matrix dA.  M >= 0.
    
    N       (input) INTEGER
            The number of columns of the matrix dA.  N >= 0.
    
    ALPHA   (input) COMPLEX DOUBLE PRECISION
            The scalar alpha.
            
    dA      (input/output) COMPLEX DOUBLE PRECISION array, dimension (LDDA,N)
            The m by n matrix dA.
    
    LDDA    (input) INTEGER
            The leading dimension of the array dA.  LDDA >= max(1,M).
            
    dB      (input/output) COMPLEX DOUBLE PRECISION array, dimension (LDDB,N)
            The m by n matrix dB.
    
    LDDB    (input) INTEGER
            The leading dimension of the array dB.  LDDB >= max(1,M).
    
    =====================================================================   */

    //printf( "m %d, grid %d, threads %d\n", m, grid.x, threads.x );
    if ( m == 0 || n == 0 )
        return;
    
    assert( m > 0 );
    assert( n > 0 );
    assert( ldda >= m );
    assert( lddb >= m );
    
    dim3 threads( NB );
    dim3 grid( (m + NB - 1)/NB );
    
    cgeadd_kernel<<< grid, threads, 0, magma_stream >>>
        ( m, n, alpha, dA, ldda, dB, lddb );
}
