/*
    -- MAGMA (version 1.2.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       May 2012

       @author Hatem Ltaief
       @author Mathieu Faverge

       @generated s Tue May 15 18:17:51 2012

*/
#include "common_magma.h"

#define magma_sgemm magmablas_sgemm
//#define magma_strsm magmablas_strsm
//#define magma_strmm magmablas_strmm

extern "C" magma_int_t
magma_sgessm_gpu( char storev, magma_int_t m, magma_int_t n, magma_int_t k, magma_int_t ib, 
                  magma_int_t *ipiv, 
                  float *dL1, magma_int_t lddl1, 
                  float *dL,  magma_int_t lddl, 
                  float *dA,  magma_int_t ldda, 
                  magma_int_t *info)
{
/*  -- MAGMA (version 1.2.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       May 2012

    Purpose
    =======

    SGETRF computes an LU factorization of a general M-by-N matrix A
    using partial pivoting with row interchanges.

    The factorization has the form
       A = P * L * U
    where P is a permutation matrix, L is lower triangular with unit
    diagonal elements (lower trapezoidal if m > n), and U is upper
    triangular (upper trapezoidal if m < n).

    This is the right-looking Level 3 BLAS version of the algorithm.

    Arguments
    =========

    M       (input) INTEGER
            The number of rows of the matrix A.  M >= 0.

    N       (input) INTEGER
            The number of columns of the matrix A.  N >= 0.

    A       (input/output) REAL array on the GPU, dimension (LDA,N).
            On entry, the M-by-N matrix to be factored.
            On exit, the factors L and U from the factorization
            A = P*L*U; the unit diagonal elements of L are not stored.

    LDA     (input) INTEGER
            The leading dimension of the array A.  LDA >= max(1,M).

    IPIV    (output) INTEGER array, dimension (min(M,N))
            The pivot indices; for 1 <= i <= min(M,N), row i of the
            matrix was interchanged with row IPIV(i).

    INFO    (output) INTEGER
            = 0:  successful exit
            < 0:  if INFO = -i, the i-th argument had an illegal value
                  or another error occured, such as memory allocation failed.
            > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
                  has been completed, but the factor U is exactly
                  singular, and division by zero will occur if it is used
                  to solve a system of equations.

    =====================================================================    */

#define AT(i,j) (dAT + (i)*ldda + (j)      )
#define L(i,j)  (dL  + (i)      + (j)*lddl )
#define dL1(j)  (dL1            + (j)*lddl1)

    float c_one     = MAGMA_S_ONE;
    float c_neg_one = MAGMA_S_NEG_ONE;

    int i, s, sb;
    float *dAT;

    /* Check arguments */
    *info = 0;
    if (m < 0)
        *info = -1;
    else if (n < 0)
        *info = -2;
    else if (ldda < max(1,m))
        *info = -4;

    if (*info != 0) {
        magma_xerbla( __func__, -(*info) );
        return *info;
    }

    /* Quick return if possible */
    if (m == 0 || n == 0)
        return *info;

    if ( (storev == 'C') || (storev == 'c') ) {
        magmablas_sgetmo_in( dA, dAT, ldda, m, n );
    } else {
        dAT = dA;
    }

    s = k / ib;
    for(i = 0; i < k; i += ib) {
        sb = min(ib, k-i);

        magmablas_slaswp( n, dAT, ldda, i+1, i+sb, ipiv, 1 );

#ifndef WITHOUTTRTRI
        magma_strmm( MagmaRight, MagmaLower, MagmaTrans, MagmaUnit, 
                     n, sb, 
                     c_one, dL1(i),   lddl1,
                            AT(i, 0), ldda);
#else
        magma_strsm( MagmaRight, MagmaLower, MagmaTrans, MagmaUnit, 
                     n, sb, 
                     c_one, L( i, i), lddl,
                            AT(i, 0), ldda);
#endif

        if ( (i+sb) < m) {
            magma_sgemm( MagmaNoTrans, MagmaTrans, 
                         n, m-(i+sb), sb, 
                         c_neg_one, AT(i,    0), ldda,
                                    L( i+sb, i), lddl, 
                         c_one,     AT(i+sb, 0), ldda );
        }
    }

    if ( (storev == 'C') || (storev == 'c') ) {
        magmablas_sgetmo_in( dA, dAT, ldda, m, n );
    }

    return *info;
    /* End of MAGMA_SGETRF_GPU */
}
