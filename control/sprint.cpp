/*
    -- MAGMA (version 1.5.0-beta3) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date July 2014

       @author Mark Gates
       @generated from zprint.cpp normal z -> s, Fri Jul 18 17:34:20 2014

*/
#include "common_magma.h"

#define PRECISION_s

/**
    Purpose
    -------

    magma_sprint prints a matrix that is located on the CPU host.
    The output is intended to be Matlab compatible, to be useful in debugging.

    Arguments
    ---------

    @param[in]
    m       INTEGER
            The number of rows of the matrix A.  M >= 0.

    @param[in]
    n       INTEGER
            The number of columns of the matrix A.  N >= 0.

    @param[in]
    A       REAL array, dimension (LDA,N), on the CPU host.
            The M-by-N matrix to be printed.

    @param[in]
    lda     INTEGER
            The leading dimension of the array A.  LDA >= max(1,M).

    @ingroup magma_saux2
    ********************************************************************/
extern "C"
void magma_sprint( magma_int_t m, magma_int_t n, const float *A, magma_int_t lda )
{
    #define A(i,j) (A + (i) + (j)*lda)
    
    magma_int_t info = 0;
    if ( m < 0 )
        info = -1;
    else if ( n < 0 )
        info = -2;
    else if ( magma_is_devptr( A ) == 1 )
        info = -3;
    else if ( lda < max(1,m) )
        info = -4;
    
    if (info != 0) {
        magma_xerbla( __func__, -(info) );
        return;  //info;
    }
    
    float c_zero = MAGMA_S_ZERO;
    
    if ( m == 1 ) {
        printf( "[ " );
    }
    else {
        printf( "[\n" );
    }
    for( int i = 0; i < m; ++i ) {
        for( int j = 0; j < n; ++j ) {
            if ( MAGMA_S_EQUAL( *A(i,j), c_zero )) {
#if defined(PRECISION_z) || defined(PRECISION_c)
                printf( "   0.              " );
#else
                printf( "   0.    " );
#endif
            }
            else {
#if defined(PRECISION_z) || defined(PRECISION_c)
                printf( " %8.4f+%8.4fi", MAGMA_S_REAL( *A(i,j) ), MAGMA_S_IMAG( *A(i,j) ));
#else
                printf( " %8.4f", MAGMA_S_REAL( *A(i,j) ));
#endif
            }
        }
        if ( m > 1 ) {
            printf( "\n" );
        }
        else {
            printf( " " );
        }
    }
    printf( "];\n" );
}


/**
    Purpose
    -------
    magma_sprint_gpu prints a matrix that is located on the GPU device.
    Internally, it allocates CPU memory and copies the matrix to the CPU.
    The output is intended to be Matlab compatible, to be useful in debugging.

    Arguments
    ---------
    @param[in]
    m       INTEGER
            The number of rows of the matrix A.  M >= 0.

    @param[in]
    n       INTEGER
            The number of columns of the matrix A.  N >= 0.

    @param[in]
    dA      REAL array, dimension (LDDA,N), on the GPU device.
            The M-by-N matrix to be printed.

    @param[in]
    ldda    INTEGER
            The leading dimension of the array A.  LDDA >= max(1,M).

    @ingroup magma_saux2
    ********************************************************************/
extern "C"
void magma_sprint_gpu( magma_int_t m, magma_int_t n, const float *dA, magma_int_t ldda )
{
    magma_int_t info = 0;
    if ( m < 0 )
        info = -1;
    else if ( n < 0 )
        info = -2;
    else if ( magma_is_devptr( dA ) == 0 )
        info = -3;
    else if ( ldda < max(1,m) )
        info = -4;
    
    if (info != 0) {
        magma_xerbla( __func__, -(info) );
        return;  //info;
    }
    
    magma_int_t lda = m;
    float* A;
    magma_smalloc_cpu( &A, lda*n );
    magma_sgetmatrix( m, n, dA, ldda, A, lda );
    
    magma_sprint( m, n, A, lda );
    
    magma_free_cpu( A );
}
