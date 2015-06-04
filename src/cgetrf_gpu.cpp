/*
    -- MAGMA (version 1.2.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       June 2012

       @generated c Thu Jun 28 12:30:37 2012

*/
#include "common_magma.h"

// === Define what BLAS to use ============================================
#define PRECISION_c
#if (defined(PRECISION_s) || defined(PRECISION_d))
  #define magma_cgemm magmablas_cgemm
  #define magma_ctrsm magmablas_ctrsm
#endif
// === End defining what BLAS to use =======================================

extern "C" magma_int_t
magma_cgetrf_gpu(magma_int_t m, magma_int_t n, 
                 cuFloatComplex *dA, magma_int_t ldda,
                 magma_int_t *ipiv, magma_int_t *info)
{
/*  -- MAGMA (version 1.2.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       June 2012

    Purpose
    =======

    CGETRF computes an LU factorization of a general M-by-N matrix A
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

    A       (input/output) COMPLEX array on the GPU, dimension (LDDA,N).
            On entry, the M-by-N matrix to be factored.
            On exit, the factors L and U from the factorization
            A = P*L*U; the unit diagonal elements of L are not stored.

    LDDA     (input) INTEGER
            The leading dimension of the array A.  LDDA >= max(1,M).

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

#define inAT(i,j) (dAT + (i)*nb*lddat + (j)*nb)

    cuFloatComplex c_one     = MAGMA_C_ONE;
    cuFloatComplex c_neg_one = MAGMA_C_NEG_ONE;

    magma_int_t iinfo, nb;
    magma_int_t maxm, maxn, mindim;
    magma_int_t i, rows, cols, s, lddat, lddwork;
    cuFloatComplex *dAT, *dAP, *work;

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
    nb     = magma_get_cgetrf_nb(m);
    s      = mindim / nb;

    if (nb <= 1 || nb >= min(m,n)) {
        /* Use CPU code. */
        magma_cmalloc_cpu( &work, m * n );
        if ( work == NULL ) {
            *info = MAGMA_ERR_HOST_ALLOC;
            return *info;
        }
        magma_cgetmatrix( m, n, dA, ldda, work, m );
        lapackf77_cgetrf(&m, &n, work, &m, ipiv, info);
        magma_csetmatrix( m, n, work, m, dA, ldda );
        free(work);
    }
    else {
        /* Use hybrid blocked code. */
        maxm = ((m + 31)/32)*32;
        maxn = ((n + 31)/32)*32;

        lddat   = maxn;
        lddwork = maxm;

        dAT = dA;

        if (MAGMA_SUCCESS != magma_cmalloc( &dAP, nb*maxm )) {
            *info = MAGMA_ERR_DEVICE_ALLOC;
            return *info;
        }

        if ((m == n) && (m % 32 == 0) && (ldda%32 == 0)){
            lddat = ldda;
            magmablas_cinplace_transpose( dAT, ldda, m);
        }
        else {
            if (MAGMA_SUCCESS != magma_cmalloc( &dAT, maxm*maxn )) {
                magma_free( dAP );
                *info = MAGMA_ERR_DEVICE_ALLOC;
                return *info;
            }
            magmablas_ctranspose2( dAT, lddat, dA, ldda, m, n );
        }

        if (MAGMA_SUCCESS != magma_cmalloc_pinned( &work, maxm*nb )) {
            magma_free( dAP );
            if (! ((m == n) && (m % 32 == 0) && (ldda%32 == 0)) )
                magma_free( dAT );
            *info = MAGMA_ERR_HOST_ALLOC;
            return *info;
        }

        for( i=0; i<s; i++ )
            {
                // download i-th panel
                cols = maxm - i*nb;
                magmablas_ctranspose( dAP, cols, inAT(i,i), lddat, nb, cols );
                magma_cgetmatrix( m-i*nb, nb, dAP, cols, work, lddwork );

                // make sure that gpu queue is empty
                magma_device_sync();

                if ( i>0 ){
                    magma_ctrsm( MagmaRight, MagmaUpper, MagmaNoTrans, MagmaUnit, 
                                 n - (i+1)*nb, nb, 
                                 c_one, inAT(i-1,i-1), lddat, 
                                        inAT(i-1,i+1), lddat );
                    magma_cgemm( MagmaNoTrans, MagmaNoTrans, 
                                 n-(i+1)*nb, m-i*nb, nb, 
                                 c_neg_one, inAT(i-1,i+1), lddat, 
                                            inAT(i,  i-1), lddat, 
                                 c_one,     inAT(i,  i+1), lddat );
                }

                // do the cpu part
                rows = m - i*nb;
                lapackf77_cgetrf( &rows, &nb, work, &lddwork, ipiv+i*nb, &iinfo);
                if ( (*info == 0) && (iinfo > 0) )
                    *info = iinfo + i*nb;

                magmablas_cpermute_long2( n, dAT, lddat, ipiv, nb, i*nb );

                // upload i-th panel
                magma_csetmatrix( m-i*nb, nb, work, lddwork, dAP, maxm );
                magmablas_ctranspose(inAT(i,i), lddat, dAP, maxm, cols, nb);

                // do the small non-parallel computations
                if ( s > (i+1) ) {
                    magma_ctrsm( MagmaRight, MagmaUpper, MagmaNoTrans, MagmaUnit, 
                                 nb, nb, 
                                 c_one, inAT(i, i  ), lddat,
                                        inAT(i, i+1), lddat);
                    magma_cgemm( MagmaNoTrans, MagmaNoTrans, 
                                 nb, m-(i+1)*nb, nb, 
                                 c_neg_one, inAT(i,   i+1), lddat,
                                            inAT(i+1, i  ), lddat, 
                                 c_one,     inAT(i+1, i+1), lddat );
                }
                else {
                    magma_ctrsm( MagmaRight, MagmaUpper, MagmaNoTrans, MagmaUnit, 
                                 n-s*nb, nb, 
                                 c_one, inAT(i, i  ), lddat,
                                        inAT(i, i+1), lddat);
                    magma_cgemm( MagmaNoTrans, MagmaNoTrans, 
                                 n-(i+1)*nb, m-(i+1)*nb, nb,
                                 c_neg_one, inAT(i,   i+1), lddat,
                                            inAT(i+1, i  ), lddat, 
                                 c_one,     inAT(i+1, i+1), lddat );
                }
            }

        magma_int_t nb0 = min(m - s*nb, n - s*nb);
        rows = m - s*nb;
        cols = maxm - s*nb;

        magmablas_ctranspose2( dAP, maxm, inAT(s,s), lddat, nb0, rows);
        magma_cgetmatrix( rows, nb0, dAP, maxm, work, lddwork );

        // make sure that gpu queue is empty
        magma_device_sync();

        // do the cpu part
        lapackf77_cgetrf( &rows, &nb0, work, &lddwork, ipiv+s*nb, &iinfo);
        if ( (*info == 0) && (iinfo > 0) )
            *info = iinfo + s*nb;
        magmablas_cpermute_long2( n, dAT, lddat, ipiv, nb0, s*nb );

        // upload i-th panel
        magma_csetmatrix( rows, nb0, work, lddwork, dAP, maxm );
        magmablas_ctranspose2( inAT(s,s), lddat, dAP, maxm, rows, nb0);

        magma_ctrsm( MagmaRight, MagmaUpper, MagmaNoTrans, MagmaUnit, 
                     n-s*nb-nb0, nb0,
                     c_one, inAT(s,s),     lddat, 
                            inAT(s,s)+nb0, lddat);

        if ((m == n) && (m % 32 == 0) && (ldda%32 == 0)){
            magmablas_cinplace_transpose( dAT, lddat, m );
        }
        else {
            magmablas_ctranspose2( dA, ldda, dAT, lddat, n, m );
            magma_free( dAT );
        }

        magma_free( dAP );
        magma_free_pinned( work );
    }
    return *info;

    /* End of MAGMA_CGETRF_GPU */
}

#undef inAT
