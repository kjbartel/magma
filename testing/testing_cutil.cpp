/*
 *  -- MAGMA (version 1.3.0) --
 *     Univ. of Tennessee, Knoxville
 *     Univ. of California, Berkeley
 *     Univ. of Colorado, Denver
 *     November 2012
 *
 * @generated c Wed Nov 14 22:54:12 2012
 *
 * Utilities for testing.
 * @author Mark Gates
 **/

#include "testings.h"

#define A(i,j)  A[i + j*lda]

// --------------------
// Make a matrix symmetric/Hermitian.
// Makes diagonal real.
// Sets Aji = conj( Aij ) for j < i, that is, copy lower triangle to upper triangle.
extern "C"
void magma_chermitian( magma_int_t N, cuFloatComplex* A, magma_int_t lda )
{
    magma_int_t i, j;
    for( i=0; i<N; ++i ) {
        A(i,i) = MAGMA_C_MAKE( MAGMA_C_REAL( A(i,i) ), 0. );
        for( j=0; j<i; ++j ) {
            A(j,i) = MAGMA_C_CNJG( A(i,j) );
        }
    }
}


// --------------------
// Make a matrix symmetric/Hermitian positive definite.
// Increases diagonal by N, and makes it real.
// Sets Aji = conj( Aij ) for j < i, that is, copy lower triangle to upper triangle.
extern "C"
void magma_chpd( magma_int_t N, cuFloatComplex* A, magma_int_t lda )
{
    magma_int_t i, j;
    for( i=0; i<N; ++i ) {
        A(i,i) = MAGMA_C_MAKE( MAGMA_C_REAL( A(i,i) ) + N, 0. );
        for( j=0; j<i; ++j ) {
            A(j,i) = MAGMA_C_CNJG( A(i,j) );
        }
    }
}
