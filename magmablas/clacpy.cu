/*
    -- MAGMA (version 1.3.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2012

       @generated c Wed Nov 14 22:53:46 2012

*/
#include "common_magma.h"

/*
    Matrix is divided into 64 x n block rows.
    Each block has 64 threads.
    Each thread copies one row, iterating across all columns.
    The bottom block of rows may be partially outside the matrix;
    if so, rows outside the matrix (row >= m) are disabled.
    
    @author Mark Gates
*/
__global__ void 
clacpy_kernel( int m, int n,
               const cuFloatComplex *A, int lda,
               cuFloatComplex       *B, int ldb ) 
{ 
    int row = blockIdx.x*64 + threadIdx.x;
    if ( row < m ) {
        A += row;
        B += row;
        const cuFloatComplex *Aend = A + lda*n;
        while( A < Aend ) {
            *B = *A;
            A += lda;
            B += ldb;
        }
    }
}


extern "C" void 
magmablas_clacpy( char uplo, magma_int_t m, magma_int_t n,
                  const cuFloatComplex *A, magma_int_t lda,
                  cuFloatComplex       *B, magma_int_t ldb )
{
/*
    Note
  ========
  - UPLO Parameter is disabled
  - Do we want to provide a generic function to the user with all the options?

  Purpose
  =======

  CLACPY copies all or part of a two-dimensional matrix A to another
  matrix B.

  Arguments
  =========

  UPLO    (input) CHARACTER*1
          Specifies the part of the matrix A to be copied to B.
          = 'U':      Upper triangular part
          = 'L':      Lower triangular part
          Otherwise:  All of the matrix A

  M       (input) INTEGER
          The number of rows of the matrix A.  M >= 0.

  N       (input) INTEGER
          The number of columns of the matrix A.  N >= 0.

  A       (input) COMPLEX DOUBLE PRECISION array, dimension (LDA,N)
          The m by n matrix A.  If UPLO = 'U', only the upper triangle
          or trapezoid is accessed; if UPLO = 'L', only the lower
          triangle or trapezoid is accessed.

  LDA     (input) INTEGER
          The leading dimension of the array A.  LDA >= max(1,M).

  B       (output) COMPLEX DOUBLE PRECISION array, dimension (LDB,N)
          On exit, B = A in the locations specified by UPLO.

  LDB     (input) INTEGER
          The leading dimension of the array B.  LDB >= max(1,M).

  =====================================================================   */

    dim3 threads( 64 );
    dim3 grid( m/64 + (m%64 != 0) );
    
    //printf( "m %d, n %d, grid %d, threads %d\n", m, n, grid.x, threads.x );
    if ( m == 0 || n == 0 )
        return;
    
    if ( (uplo == 'U') || (uplo == 'u') ) {
        fprintf(stderr, "lacpy upper is not implemented\n");
    }
    else if ( (uplo == 'L') || (uplo == 'l') ) {
        fprintf(stderr, "lacpy lower is not implemented\n");
    }
    else {
        clacpy_kernel<<< grid, threads, 0, magma_stream >>> ( m, n, A, lda, B, ldb );
    }
}
