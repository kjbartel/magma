/*
    -- MAGMA (version 1.2.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       May 2012

       @author Hatem Ltaief
       @author Mathieu Faverge

       @generated c Tue May 15 18:17:51 2012

*/
#ifdef MAGMA_WITH_PLASMA

#include <plasma.h>
#include <core_blas.h>
#include "common_magma.h"

#define magma_cgemm magmablas_cgemm
//#define magma_ctrsm magmablas_ctrsm
//#define magma_ctrmm magmablas_ctrmm

extern "C" magma_int_t
magma_cgetrl_gpu( char storev, magma_int_t m, magma_int_t n, magma_int_t ib,
                  cuFloatComplex *hA, magma_int_t ldha, cuFloatComplex *dA, magma_int_t ldda,
                  cuFloatComplex *hL, magma_int_t ldhl, cuFloatComplex *dL, magma_int_t lddl,
                  magma_int_t *ipiv, 
                  cuFloatComplex *dwork, magma_int_t lddwork,
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


    dwork   dimension(LDDWORK, IB) to store a single panel.

    LDDWORK > max(1, ((MB+31)/32)*32)

    INFO    (output) INTEGER
            = 0:  successful exit
            < 0:  if INFO = -i, the i-th argument had an illegal value
                  or another error occured, such as memory allocation failed.
            > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
                  has been completed, but the factor U is exactly
                  singular, and division by zero will occur if it is used
                  to solve a system of equations.

    =====================================================================    */

#define AT(i,j) (dAT + (i)*ib*ldda + (j)*ib)
#define hA(i,j) (hA  + (i)*ib + (j)*ib*ldha)
#define hL(j)   (hL  + (j)*ib*ldhl         )
#define hL2(j)  (hL2 + (j)*ib*ldhl         )
#define dL(j)   (dL  + (j)*ib*lddl         )
#define dL2(j)  (dL2 + (j)*ib*lddl         )

    cuFloatComplex c_one     = MAGMA_C_ONE;
    cuFloatComplex c_neg_one = MAGMA_C_NEG_ONE;

    magma_int_t iinfo;
    magma_int_t maxm, mindim;
    magma_int_t i, rows, cols, s, ii, sb;
    cuFloatComplex *dAT;
#ifndef WITHOUTTRTRI
    cuFloatComplex *dL2 = dL + ib;
    cuFloatComplex *hL2 = hL + ib;
#endif

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

    /* Function Body */
    mindim = min(m, n);
    s      = mindim / ib;

    if ( ib >= mindim ) {
        /* Use CPU code. */
        lapackf77_cgetrf(&m, &n, hA, &ldha, ipiv, info);

#ifndef WITHOUTTRTRI
        CORE_clacpy(PlasmaUpperLower, mindim, mindim, 
                    (PLASMA_Complex32_t*)hA, ldha, 
                    (PLASMA_Complex32_t*)hL2, ldhl );

        CORE_ctrtri( PlasmaLower, PlasmaUnit, mindim, 
                     (PLASMA_Complex32_t*)hL2, ldhl, info );
        if (*info != 0 ) {
          fprintf(stderr, "ERROR, trtri returned with info = %d\n", *info);
        }          

        magma_csetmatrix( mindim, mindim, hL2, ldhl, dL2, lddl );
#endif
            
        if ( (storev == 'R') || (storev == 'r') ) {
            magma_csetmatrix( m, n, hA, ldha, dwork, lddwork );
            magmablas_ctranspose( dA, ldda, dwork, lddwork, m, n );
        } else {
            magma_csetmatrix( m, n, hA, ldha, dA, ldda );
        }
    }
    else {
        /* Use hybrid blocked code. */
        maxm = ((m + 31)/32)*32;

        if ( (storev == 'C') || (storev == 'c') ) {
            magmablas_cgetmo_in( dA, dAT, ldda, m, n );
        } else {
            dAT = dA;
        }
            
        for( i=0; i<s; i++ )
        {
            ii = i * ib;
            sb = min(ib, mindim-ii);
            cols = maxm - ii;

            if ( i>0 ){
                // download i-th panel
                magmablas_ctranspose( dwork, maxm, AT(0, i), ldda, sb, m );
                magma_cgetmatrix( m, sb, dwork, maxm, hA(0, i), ldha );
                
                // make sure that gpu queue is empty
                //magma_device_sync();
#ifndef WITHOUTTRTRI
                magma_ctrmm( MagmaRight, MagmaLower, MagmaTrans, MagmaUnit, 
                             n - (ii+sb), ib, 
                             c_one, dL2(i-1),    lddl, 
                                    AT(i-1,i+1), ldda );
#else
                magma_ctrsm( MagmaRight, MagmaUpper, MagmaNoTrans, MagmaUnit, 
                             n - (ii+sb), ib, 
                             c_one, AT(i-1,i-1), ldda, 
                                    AT(i-1,i+1), ldda );
#endif
                magma_cgemm( MagmaNoTrans, MagmaNoTrans, 
                             n-(ii+sb), m-ii, ib, 
                             c_neg_one, AT(i-1,i+1), ldda, 
                                        AT(i,  i-1), ldda, 
                             c_one,     AT(i,  i+1), ldda );
            }

            // do the cpu part
            rows = m - ii;
            lapackf77_cgetrf( &rows, &sb, hA(i, i), &ldha, ipiv+ii, &iinfo);
            if ( (*info == 0) && (iinfo > 0) )
                *info = iinfo + ii;

            { 
                int j;
                int fin = ii + sb;
                for(j=ii ; j <fin; j++) {
                    ipiv[j] = ii + ipiv[j];
                }
            }
            magmablas_claswp( n-ii, AT(0, i), ldda, ii+1, ii+sb, ipiv, 1 );

#ifndef WITHOUTTRTRI
            CORE_clacpy(PlasmaLower, sb, sb, 
                        (PLASMA_Complex32_t*)hA(i, i), ldha, 
                        (PLASMA_Complex32_t*)hL2(i), ldhl );
            
            CORE_ctrtri( PlasmaLower, PlasmaUnit, sb, 
                         (PLASMA_Complex32_t*)hL2(i), ldhl, info );
            if (*info != 0 ) {
              fprintf(stderr, "ERROR, trtri returned with info = %d\n", *info);
            }
            magma_csetmatrix( sb, sb, hL2(i), ldhl, dL2(i), lddl );
#endif
            // upload i-th panel
            magma_csetmatrix( rows, sb, hA(i, i), ldha, dwork, cols );
            magmablas_ctranspose( AT(i,i), ldda, dwork, cols, rows, sb);

            // do the small non-parallel computations
            if ( s > (i+1) ) {
#ifndef WITHOUTTRTRI
                magma_ctrmm( MagmaRight, MagmaLower, MagmaTrans, MagmaUnit, 
                             sb, sb, 
                             c_one, dL2(i),     lddl,
                                    AT(i, i+1), ldda);
#else
                magma_ctrsm( MagmaRight, MagmaUpper, MagmaNoTrans, MagmaUnit, 
                             sb, sb, 
                             c_one, AT(i, i  ), ldda,
                                    AT(i, i+1), ldda);
#endif
                magma_cgemm( MagmaNoTrans, MagmaNoTrans, 
                             sb, m-(ii+sb), sb, 
                             c_neg_one, AT(i,   i+1), ldda,
                                        AT(i+1, i  ), ldda, 
                             c_one,     AT(i+1, i+1), ldda );
            }
            else {
                /* Update of the last panel */
#ifndef WITHOUTTRTRI
                magma_ctrmm( MagmaRight, MagmaLower, MagmaTrans, MagmaUnit, 
                             n-mindim, sb, 
                             c_one, dL2(i),     lddl,
                                    AT(i, i+1), ldda);
#else
                magma_ctrsm( MagmaRight, MagmaUpper, MagmaNoTrans, MagmaUnit, 
                             n-mindim, sb, 
                             c_one, AT(i, i  ), ldda,
                                    AT(i, i+1), ldda);
#endif
                /* m-(ii+sb) should be always 0 */
                magma_cgemm( MagmaNoTrans, MagmaNoTrans, 
                             n-mindim, m-(ii+sb), sb,
                             c_neg_one, AT(i,   i+1), ldda,
                                        AT(i+1, i  ), ldda, 
                             c_one,     AT(i+1, i+1), ldda );
            }
        }

        if ( (storev == 'C') || (storev == 'c') ) {
            magmablas_cgetmo_out( dA, dAT, ldda, m, n );
        }
    }
    return *info;
}

#endif
