/*
    -- MAGMA (version 1.5.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date September 2014

       @generated from magma_zutil.cpp normal z -> c, Wed Sep 17 15:08:40 2014

       @author Mark Gates

       Utilities for testing.
*/

#include "testings.h"

#define A(i,j)  A[i + j*lda]

// --------------------
// Make a matrix symmetric/Hermitian.
// Makes diagonal real.
// Sets Aji = conj( Aij ) for j < i, that is, copy lower triangle to upper triangle.
extern "C"
void magma_cmake_hermitian( magma_int_t N, magmaFloatComplex* A, magma_int_t lda )
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
void magma_cmake_hpd( magma_int_t N, magmaFloatComplex* A, magma_int_t lda )
{
    magma_int_t i, j;
    for( i=0; i<N; ++i ) {
        A(i,i) = MAGMA_C_MAKE( MAGMA_C_REAL( A(i,i) ) + N, 0. );
        for( j=0; j<i; ++j ) {
            A(j,i) = MAGMA_C_CNJG( A(i,j) );
        }
    }
}
