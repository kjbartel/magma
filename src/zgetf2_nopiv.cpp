/*
    -- MAGMA (version 1.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2011

       @precisions normal z -> s d c

*/
#include "common_magma.h"

// === Define what BLAS to use ============================================
#define PRECISION_z
#if (defined(PRECISION_s) || defined(PRECISION_d))
  #define cublasZgemm magmablas_zgemm
  #define cublasZtrsm magmablas_ztrsm
#endif
// === End defining what BLAS to use =======================================

extern "C" magma_int_t
magma_zgetf2_nopiv(magma_int_t *m, magma_int_t *n, cuDoubleComplex *a, 
                   magma_int_t *lda, magma_int_t *info)
{
/*  -- MAGMA (version 1.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2011

    Purpose   
    =======   
    ZGETF2_NOPIV computes an LU factorization of a general m-by-n
    matrix A without pivoting.

    The factorization has the form   
       A = L * U   
    where L is lower triangular with unit diagonal elements (lower
    trapezoidal if m > n), and U is upper triangular (upper 
    trapezoidal if m < n).   

    This is the right-looking Level 2 BLAS version of the algorithm.   

    Arguments   
    =========   
    M       (input) INTEGER   
            The number of rows of the matrix A.  M >= 0.   

    N       (input) INTEGER   
            The number of columns of the matrix A.  N >= 0.   

    A       (input/output) COMPLEX*16 array, dimension (LDA,N)   
            On entry, the m by n matrix to be factored.   
            On exit, the factors L and U from the factorization   
            A = P*L*U; the unit diagonal elements of L are not stored.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= max(1,M).   

    INFO    (output) INTEGER   
            = 0: successful exit   
            < 0: if INFO = -k, the k-th argument had an illegal value   
            > 0: if INFO = k, U(k,k) is exactly zero. The factorization   
                 has been completed, but the factor U is exactly   
                 singular, and division by zero will occur if it is used   
                 to solve a system of equations.   
    =====================================================================   */
    
  cuDoubleComplex c_one = MAGMA_Z_ONE, c_zero = MAGMA_Z_ZERO;
    static magma_int_t c__1 = 1;
    
    magma_int_t a_dim1, a_offset, i__1, i__2, i__3;
    cuDoubleComplex z__1;
    static magma_int_t i__, j;
    static double sfmin;

    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;

    /* Function Body */
    *info = 0;
    if (*m < 0) {
        *info = -1;
    } else if (*n < 0) {
        *info = -2;
    } else if (*lda < max(1,*m)) {
        *info = -4;
    }
    if (*info != 0) {
        magma_xerbla( __func__, -(*info) );
        return MAGMA_ERR_ILLEGAL_VALUE;
    }

    /* Quick return if possible */
    if (*m == 0 || *n == 0)
        return MAGMA_SUCCESS;

    /* Compute machine safe minimum */
    sfmin = lapackf77_dlamch("S");

    i__1 = min(*m,*n);
    for (j = 1; j <= i__1; ++j) 
      {
        /* Test for singularity. */
        i__2 = j + j * a_dim1;
        if (!MAGMA_Z_EQUAL(a[i__2], c_zero)) {
    
          /* Compute elements J+1:M of J-th column. */
          if (j < *m) {
            if (MAGMA_Z_ABS(a[j + j * a_dim1]) >= sfmin) 
              {
                i__2 = *m - j;
                MAGMA_Z_ASSIGN(z__1, MAGMA_Z_DIV(c_one, a[j + j * a_dim1]));
                blasf77_zscal(&i__2, &z__1, &a[j + 1 + j * a_dim1], &c__1);
              } 
            else 
              {
                i__2 = *m - j;
                for (i__ = 1; i__ <= i__2; ++i__) {
                  i__3 = j + i__ + j * a_dim1;
                  MAGMA_Z_ASSIGN(a[i__3], MAGMA_Z_DIV(a[j + i__ + j * a_dim1], a[j + j*a_dim1]));
                }
              }
          }
          
        } else if (*info == 0)
          *info = j;
        
        if (j < min(*m,*n)) {  
          /* Update trailing submatrix. */
          i__2 = *m - j;
          i__3 = *n - j;
          z__1 = MAGMA_Z_NEG_ONE; 
          zgeru_(&i__2, &i__3, &z__1, &a[j + 1 + j * a_dim1], &c__1,
                 &a[j + (j+1) * a_dim1], lda, &a[j + 1 + (j+1) * a_dim1], lda);
        }
      }
    
    return MAGMA_SUCCESS;
} /* magma_zgetf2_nopiv */
