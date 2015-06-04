/*
    -- MAGMA (version 1.4.0-beta2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       June 2013

       @precisions normal z -> s d c

*/
#include "common_magma.h"


// =========================================================================
// definitions of non-GPU-resident multi-GPU subroutines
/* non-gpu-resident interface to multiple GPUs */
extern "C" magma_int_t
magma_zgetrf_m(magma_int_t num_gpus0, magma_int_t m, magma_int_t n, magmaDoubleComplex *a, magma_int_t lda,
               magma_int_t *ipiv, magma_int_t *info);

/* to apply pivoting from the previous big panel on CPU */
extern "C" magma_int_t
magma_zgetrf_piv(magma_int_t num_gpus, magma_int_t m, magma_int_t n, magmaDoubleComplex *a, magma_int_t lda,
                 magma_int_t *ipiv, magma_int_t *info);
// =========================================================================


extern "C" magma_int_t
magma_zgetrf(magma_int_t m, magma_int_t n, magmaDoubleComplex *a, magma_int_t lda,
             magma_int_t *ipiv, magma_int_t *info)
{
/*  -- MAGMA (version 1.4.0-beta2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       June 2013

    Purpose
    =======
    ZGETRF computes an LU factorization of a general M-by-N matrix A
    using partial pivoting with row interchanges.  This version does not
    require work space on the GPU passed as input. GPU memory is allocated
    in the routine.

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

    A       (input/output) COMPLEX_16 array, dimension (LDA,N)
            On entry, the M-by-N matrix to be factored.
            On exit, the factors L and U from the factorization
            A = P*L*U; the unit diagonal elements of L are not stored.

            Higher performance is achieved if A is in pinned memory, e.g.
            allocated using magma_malloc_pinned.

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

#define inAT(i,j) (dAT + (i)*nb*ldda + (j)*nb)

    magmaDoubleComplex *dAT, *dA, *da, *work;
    magmaDoubleComplex c_one     = MAGMA_Z_ONE;
    magmaDoubleComplex c_neg_one = MAGMA_Z_NEG_ONE;
    magma_int_t     iinfo, nb;

    *info = 0;

    if (m < 0)
        *info = -1;
    else if (n < 0)
        *info = -2;
    else if (lda < max(1,m))
        *info = -4;

    if (*info != 0) {
        magma_xerbla( __func__, -(*info) );
        return *info;
    }

    /* Quick return if possible */
    if (m == 0 || n == 0)
        return *info;

    nb = magma_get_zgetrf_nb(m);

    if ( (nb <= 1) || (nb >= min(m,n)) ) {
        /* Use CPU code. */
        lapackf77_zgetrf(&m, &n, a, &lda, ipiv, info);
    } else {
        /* Use hybrid blocked code. */
        magma_int_t maxm, maxn, ldda, maxdim;
        magma_int_t i, rows, cols, s = min(m, n)/nb;
        
        magma_int_t num_gpus = magma_num_gpus();
        if ( num_gpus > 1 ) {
            /* call multi-GPU non-GPU-resident interface  */
            magma_zgetrf_m(num_gpus, m, n, a, lda, ipiv, info);
            if( *info >= 0 ) magma_zgetrf_piv(num_gpus, m, n, a, lda, ipiv, info);
            return *info;
        }

        maxm = ((m + 31)/32)*32;
        maxn = ((n + 31)/32)*32;
        maxdim = max(maxm, maxn);

        ldda = maxn;
        work = a;

        if (maxdim*maxdim < 2*maxm*maxn) {
            // if close to square, allocate square matrix and transpose in-place
            if (MAGMA_SUCCESS != magma_zmalloc( &dA, nb*maxm + maxdim*maxdim )) {
                        /* alloc failed so call non-GPU-resident version */
                        magma_zgetrf_m(num_gpus, m, n, a, lda, ipiv, info);
                        if( *info >= 0 ) magma_zgetrf_piv(num_gpus, m, n, a, lda, ipiv, info);
                        return *info;
            }
            da = dA + nb*maxm;
            
            ldda = maxdim;
            magma_zsetmatrix( m, n, a, lda, da, ldda );
            
            dAT = da;
            magmablas_ztranspose_inplace( ldda, dAT, ldda );
        }
        else {
            // if very rectangular, allocate dA and dAT and transpose out-of-place
            if (MAGMA_SUCCESS != magma_zmalloc( &dA, (nb + maxn)*maxm )) {
                        /* alloc failed so call non-GPU-resident version */
                        magma_zgetrf_m(num_gpus, m, n, a, lda, ipiv, info);
                        if( *info >= 0 ) magma_zgetrf_piv(num_gpus, m, n, a, lda, ipiv, info);
                        return *info;
            }
            da = dA + nb*maxm;
            
            magma_zsetmatrix( m, n, a, lda, da, maxm );
            
            if (MAGMA_SUCCESS != magma_zmalloc( &dAT, maxm*maxn )) {
                        /* alloc failed so call non-GPU-resident version */
                        magma_free( dA );
                        magma_zgetrf_m(num_gpus, m, n, a, lda, ipiv, info);
                        if( *info >= 0 ) magma_zgetrf_piv(num_gpus, m, n, a, lda, ipiv, info);
                        return *info;
            }

            magmablas_ztranspose2( dAT, ldda, da, maxm, m, n );
        }
        
        lapackf77_zgetrf( &m, &nb, work, &lda, ipiv, &iinfo);

        for( i = 0; i < s; i++ )
        {
            // download i-th panel
            cols = maxm - i*nb;
            
            if (i>0){
                magmablas_ztranspose( dA, cols, inAT(i,i), ldda, nb, cols );
                magma_zgetmatrix( m-i*nb, nb, dA, cols, work, lda );
                
                // make sure that gpu queue is empty
                magma_device_sync();
                
                magma_ztrsm( MagmaRight, MagmaUpper, MagmaNoTrans, MagmaUnit,
                             n - (i+1)*nb, nb,
                             c_one, inAT(i-1,i-1), ldda,
                                    inAT(i-1,i+1), ldda );
                magma_zgemm( MagmaNoTrans, MagmaNoTrans,
                             n-(i+1)*nb, m-i*nb, nb,
                             c_neg_one, inAT(i-1,i+1), ldda,
                                        inAT(i,  i-1), ldda,
                             c_one,     inAT(i,  i+1), ldda );

                // do the cpu part
                rows = m - i*nb;
                lapackf77_zgetrf( &rows, &nb, work, &lda, ipiv+i*nb, &iinfo);
            }
            if (*info == 0 && iinfo > 0)
                *info = iinfo + i*nb;
            magmablas_zpermute_long2( ldda, dAT, ldda, ipiv, nb, i*nb );

            // upload i-th panel
            magma_zsetmatrix( m-i*nb, nb, work, lda, dA, cols );
            magmablas_ztranspose( inAT(i,i), ldda, dA, cols, cols, nb);

            // do the small non-parallel computations
            if (s > (i+1)){
                magma_ztrsm( MagmaRight, MagmaUpper, MagmaNoTrans, MagmaUnit,
                             nb, nb,
                             c_one, inAT(i, i  ), ldda,
                                    inAT(i, i+1), ldda);
                magma_zgemm( MagmaNoTrans, MagmaNoTrans,
                             nb, m-(i+1)*nb, nb,
                             c_neg_one, inAT(i,   i+1), ldda,
                                        inAT(i+1, i  ), ldda,
                             c_one,     inAT(i+1, i+1), ldda );
            }
            else{
                magma_ztrsm( MagmaRight, MagmaUpper, MagmaNoTrans, MagmaUnit,
                             n-s*nb, nb,
                             c_one, inAT(i, i  ), ldda,
                                    inAT(i, i+1), ldda);
                magma_zgemm( MagmaNoTrans, MagmaNoTrans,
                             n-(i+1)*nb, m-(i+1)*nb, nb,
                             c_neg_one, inAT(i,   i+1), ldda,
                                        inAT(i+1, i  ), ldda,
                             c_one,     inAT(i+1, i+1), ldda );
            }
        }
        
        magma_int_t nb0 = min(m - s*nb, n - s*nb);
        if ( nb0 > 0 ) {
            rows = m - s*nb;
            cols = maxm - s*nb;
    
            magmablas_ztranspose2( dA, cols, inAT(s,s), ldda, nb0, rows);
            magma_zgetmatrix( rows, nb0, dA, cols, work, lda );
    
            // make sure that gpu queue is empty
            magma_device_sync();
    
            // do the cpu part
            lapackf77_zgetrf( &rows, &nb0, work, &lda, ipiv+s*nb, &iinfo);
            if (*info == 0 && iinfo > 0)
                *info = iinfo + s*nb;
            magmablas_zpermute_long2( ldda, dAT, ldda, ipiv, nb0, s*nb );
    
            magma_zsetmatrix( rows, nb0, work, lda, dA, cols );
            magmablas_ztranspose2( inAT(s,s), ldda, dA, cols, rows, nb0);
    
            magma_ztrsm( MagmaRight, MagmaUpper, MagmaNoTrans, MagmaUnit,
                         n-s*nb-nb0, nb0,
                         c_one, inAT(s, s),     ldda,
                                inAT(s, s)+nb0, ldda);
        }
        
        if (maxdim*maxdim < 2*maxm*maxn) {
            magmablas_ztranspose_inplace( ldda, dAT, ldda );
            magma_zgetmatrix( m, n, da, ldda, a, lda );
        } else {
            magmablas_ztranspose2( da, maxm, dAT, ldda, n, m );
            magma_zgetmatrix( m, n, da, maxm, a, lda );
            magma_free( dAT );
        }

        magma_free( dA );
    }
    
    return *info;
} /* magma_zgetrf */

#undef inAT
