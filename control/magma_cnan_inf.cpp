/*
    -- MAGMA (version 1.5.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date September 2014

       @author Mark Gates
       @generated from magma_znan_inf.cpp normal z -> c, Wed Sep 17 15:08:33 2014

*/
#include <limits>

#include "common_magma.h"

#define COMPLEX


const magmaFloatComplex MAGMA_C_NAN
    = MAGMA_C_MAKE( std::numeric_limits<float>::quiet_NaN(),
                    std::numeric_limits<float>::quiet_NaN() );

const magmaFloatComplex MAGMA_C_INF
    = MAGMA_C_MAKE( std::numeric_limits<float>::infinity(),
                    std::numeric_limits<float>::infinity() );

/** @return true if either real(x) or imag(x) is NAN. */
inline bool magma_c_isnan( magmaFloatComplex x )
{
#ifdef COMPLEX
    return isnan( MAGMA_C_REAL( x )) ||
           isnan( MAGMA_C_IMAG( x ));
#else
    return isnan( x );
#endif
}


/** @return true if either real(x) or imag(x) is INF. */
inline bool magma_c_isinf( magmaFloatComplex x )
{
#ifdef COMPLEX
    return isinf( MAGMA_C_REAL( x )) ||
           isinf( MAGMA_C_IMAG( x ));
#else
    return isinf( x );
#endif
}


/**
    Purpose
    -------

    magma_cnan_inf checks a matrix that is located on the CPU host
    for NAN (not-a-number) and INF (infinity) values.
    
    NAN is created by 0/0 and similar.
    INF is created by x/0 and similar, where x != 0.

    Arguments
    ---------
    @param[in]
    uplo    magma_uplo_t
            Specifies what part of the matrix A to check.
      -     = MagmaUpper:  Upper triangular part of A
      -     = MagmaLower:  Lower triangular part of A
      -     = MagmaFull:   All of A

    @param[in]
    m       INTEGER
            The number of rows of the matrix A.  M >= 0.

    @param[in]
    n       INTEGER
            The number of columns of the matrix A.  N >= 0.

    @param[in]
    A       COMPLEX array, dimension (LDA,N), on the CPU host.
            The M-by-N matrix to be printed.

    @param[in]
    lda     INTEGER
            The leading dimension of the array A.  LDA >= max(1,M).

    @param[out]
    cnt_nan INTEGER*
            If non-NULL, on exit contains the number of NAN values in A.

    @param[out]
    cnt_inf INTEGER*
            If non-NULL, on exit contains the number of INF values in A.
            
    @return
      -     >= 0:  Returns number of NAN + number of INF values.
      -     <  0:  If it returns -i, the i-th argument had an illegal value,
                   or another error occured, such as memory allocation failed.

    @ingroup magma_caux2
    ********************************************************************/
extern "C"
magma_int_t magma_cnan_inf(
    magma_uplo_t uplo, magma_int_t m, magma_int_t n,
    const magmaFloatComplex *A, magma_int_t lda,
    magma_int_t *cnt_nan,
    magma_int_t *cnt_inf )
{
    #define A(i,j) (A + (i) + (j)*lda)
    
    magma_int_t info = 0;
    if ( uplo != MagmaLower && uplo != MagmaUpper && uplo != MagmaFull )
        info = -1;
    else if ( m < 0 )
        info = -2;
    else if ( n < 0 )
        info = -3;
    else if ( magma_is_devptr( A ) == 1 )
        info = -4;
    else if ( lda < max(1,m) )
        info = -5;
    
    if (info != 0) {
        magma_xerbla( __func__, -(info) );
        return info;
    }
    
    int c_nan = 0;
    int c_inf = 0;
    
    if ( uplo == MagmaLower ) {
        for( int j = 0; j < n; ++j ) {
            for( int i = j; i < m; ++i ) {  // i >= j
                if      ( magma_c_isnan( *A(i,j) )) { c_nan++; }
                else if ( magma_c_isinf( *A(i,j) )) { c_inf++; }
            }
        }
    }
    else if ( uplo == MagmaUpper ) {
        for( int j = 0; j < n; ++j ) {
            for( int i = 0; i < m && i <= j; ++i ) {  // i <= j
                if      ( magma_c_isnan( *A(i,j) )) { c_nan++; }
                else if ( magma_c_isinf( *A(i,j) )) { c_inf++; }
            }
        }
    }
    else if ( uplo == MagmaFull ) {
        for( int j = 0; j < n; ++j ) {
            for( int i = 0; i < m; ++i ) {
                if      ( magma_c_isnan( *A(i,j) )) { c_nan++; }
                else if ( magma_c_isinf( *A(i,j) )) { c_inf++; }
            }
        }
    }
    
    if ( cnt_nan != NULL ) { *cnt_nan = c_nan; }
    if ( cnt_inf != NULL ) { *cnt_inf = c_inf; }
    
    return (c_nan + c_inf);
}


/**
    Purpose
    -------

    magma_cnan_inf checks a matrix that is located on the CPU host
    for NAN (not-a-number) and INF (infinity) values.
    
    NAN is created by 0/0 and similar.
    INF is created by x/0 and similar, where x != 0.

    Arguments
    ---------
    @param[in]
    uplo    magma_uplo_t
            Specifies what part of the matrix A to check.
      -     = MagmaUpper:  Upper triangular part of A
      -     = MagmaLower:  Lower triangular part of A
      -     = MagmaFull:   All of A

    @param[in]
    m       INTEGER
            The number of rows of the matrix A.  M >= 0.

    @param[in]
    n       INTEGER
            The number of columns of the matrix A.  N >= 0.

    @param[in]
    dA      COMPLEX array, dimension (LDDA,N), on the GPU device.
            The M-by-N matrix to be printed.

    @param[in]
    ldda    INTEGER
            The leading dimension of the array A.  LDDA >= max(1,M).

    @param[out]
    cnt_nan INTEGER*
            If non-NULL, on exit contains the number of NAN values in A.

    @param[out]
    cnt_inf INTEGER*
            If non-NULL, on exit contains the number of INF values in A.
            
    @return
      -     >= 0:  Returns number of NAN + number of INF values.
      -     <  0:  If it returns -i, the i-th argument had an illegal value,
                   or another error occured, such as memory allocation failed.

    @ingroup magma_caux2
    ********************************************************************/
extern "C"
magma_int_t magma_cnan_inf_gpu(
    magma_uplo_t uplo, magma_int_t m, magma_int_t n,
    const magmaFloatComplex *dA, magma_int_t ldda,
    magma_int_t *cnt_nan,
    magma_int_t *cnt_inf )
{
    magma_int_t info = 0;
    if ( uplo != MagmaLower && uplo != MagmaUpper && uplo != MagmaFull )
        info = -1;
    else if ( m < 0 )
        info = -2;
    else if ( n < 0 )
        info = -3;
    else if ( magma_is_devptr( dA ) == 0 )
        info = -4;
    else if ( ldda < max(1,m) )
        info = -5;
    
    if (info != 0) {
        magma_xerbla( __func__, -(info) );
        return info;
    }
    
    magma_int_t lda = m;
    magmaFloatComplex* A;
    magma_cmalloc_cpu( &A, lda*n );
    magma_cgetmatrix( m, n, dA, ldda, A, lda );
    
    magma_int_t cnt = magma_cnan_inf( uplo, m, n, A, lda, cnt_nan, cnt_inf );
    
    magma_free_cpu( A );
    return cnt;
}
