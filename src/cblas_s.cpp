/*
    -- MAGMA (version 1.6.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date November 2014
 
       @author Mark Gates
       @generated from cblas_z.cpp normal z -> s, Sat Nov 15 19:54:03 2014

    Wrappers around a few CBLAS functions.
    
    Primarily, we use the standard Fortran BLAS interface in MAGMA. However,
    functions that return a value (as opposed to subroutines that are void)
    are not portable, as they depend on how Fortran returns values. The routines
    here provide a portable interface. These are not identical to CBLAS, in
    particular, [cz]dot[uc] return real numbers (as in Fortran BLAS) rather
    than return values via an argument.
    
    Only these BLAS-1 functions are provided:
    
    magma_cblas_sasum / dasum
    magma_cblas_snrm2 / dnrm2
    magma_cblas_sdot  / ddot
    magma_cblas_sdot  / ddot

*/
#include <cblas.h>

#include "common_magma.h"

#define REAL

// ========================================
// Level 1 BLAS

// --------------------
/** Returns the sum of absolute values of vector x; i.e., one norm.

    @param[in]
    n       Number of elements in vector x. n >= 0.

    @param[in]
    x       REAL array on CPU host.
            The n element vector x of dimension (1 + (n-1)*incx).

    @param[in]
    incx    Stride between consecutive elements of x. incx > 0.

    @ingroup magma_sblas1
*/
extern "C"
float magma_cblas_sasum(
    magma_int_t n,
    const float *x, magma_int_t incx )
{
    return cblas_sasum( n, x, incx );
}

// --------------------
/** Returns 2-norm of vector x. Avoids unnecesary over/underflow.

    @param[in]
    n       Number of elements in vector x. n >= 0.

    @param[in]
    x       REAL array on CPU host.
            The n element vector x of dimension (1 + (n-1)*incx).

    @param[in]
    incx    Stride between consecutive elements of x. incx > 0.

    @ingroup magma_sblas1
*/
extern "C"
float magma_cblas_snrm2(
    magma_int_t n,
    const float *x, magma_int_t incx )
{
    return cblas_snrm2( n, x, incx );
}

// --------------------
/** Returns dot product of vectors x and y; \f$ x^H y \f$.

    @param[in]
    n       Number of elements in vector x and y. n >= 0.

    @param[in]
    x       REAL array on CPU host.
            The n element vector x of dimension (1 + (n-1)*incx).

    @param[in]
    incx    Stride between consecutive elements of x. incx > 0.

    @param[in]
    y       REAL array on CPU host.
            The n element vector y of dimension (1 + (n-1)*incy).

    @param[in]
    incy    Stride between consecutive elements of dy. incy > 0.

    @ingroup magma_sblas1
*/
extern "C"
float magma_cblas_sdot(
    magma_int_t n,
    const float *x, magma_int_t incx,
    const float *y, magma_int_t incy )
{
    // after too many issues with MKL and other BLAS, just write our own dot product!
    float value = MAGMA_S_ZERO;
    magma_int_t i;
    if ( incx == 1 && incy == 1 ) {
        for( i=0; i < n; ++i ) {
            value += conj( x[i] ) * y[i];
        }
    }
    else {
        magma_int_t ix=0, iy=0;
        if ( incx < 0 ) { ix = (-n + 1)*incx + 1; }
        if ( incy < 0 ) { iy = (-n + 1)*incy + 1; }
        for( magma_int_t i=0; i < n; ++i ) {
            value += conj( x[ix] ) * y[iy];
            ix += incx;
            iy += incy;
        }
    }
    return value;
    //#ifdef COMPLEX
    //float value;
    //cblas_sdot_sub( n, x, incx, y, incy, &value );
    //return value;
    //#else
    //return cblas_sdot( n, x, incx, y, incy );
    //#endif
}

#ifdef COMPLEX
// --------------------
/** Returns dot product (unconjugated) of vectors x and y; \f$ x^T y \f$.

    @param[in]
    n       Number of elements in vector x and y. n >= 0.

    @param[in]
    x       REAL array on CPU host.
            The n element vector x of dimension (1 + (n-1)*incx).

    @param[in]
    incx    Stride between consecutive elements of x. incx > 0.

    @param[in]
    y       REAL array on CPU host.
            The n element vector y of dimension (1 + (n-1)*incy).

    @param[in]
    incy    Stride between consecutive elements of dy. incy > 0.

    @ingroup magma_sblas1
*/
extern "C"
float magma_cblas_sdot(
    magma_int_t n,
    const float *x, magma_int_t incx,
    const float *y, magma_int_t incy )
{
    // after too many issues with MKL and other BLAS, just write our own dot product!
    float value = MAGMA_S_ZERO;
    magma_int_t i;
    if ( incx == 1 && incy == 1 ) {
        for( i=0; i < n; ++i ) {
            value += x[i] * y[i];
        }
    }
    else {
        magma_int_t ix=0, iy=0;
        if ( incx < 0 ) { ix = (-n + 1)*incx + 1; }
        if ( incy < 0 ) { iy = (-n + 1)*incy + 1; }
        for( magma_int_t i=0; i < n; ++i ) {
            value += x[ix] * y[iy];
            ix += incx;
            iy += incy;
        }
    }
    return value;
    //#ifdef COMPLEX
    //float value;
    //cblas_sdot_sub( n, x, incx, y, incy, &value );
    //return value;
    //#else
    //return cblas_sdot( n, x, incx, y, incy );
    //#endif
}
#endif

#undef REAL
