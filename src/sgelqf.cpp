/*
    -- MAGMA (version 1.4.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       August 2013

       @generated s Tue Aug 13 16:44:25 2013

*/
#include "common_magma.h"

extern "C" magma_int_t
magma_sgelqf( magma_int_t m, magma_int_t n,
              float *a,    magma_int_t lda,   float *tau,
              float *work, magma_int_t lwork, magma_int_t *info)
{
/*  -- MAGMA (version 1.4.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       August 2013

    Purpose
    =======
    SGELQF computes an LQ factorization of a REAL M-by-N matrix A:
    A = L * Q.

    Arguments
    =========
    M       (input) INTEGER
            The number of rows of the matrix A.  M >= 0.

    N       (input) INTEGER
            The number of columns of the matrix A.  N >= 0.

    A       (input/output) REAL array, dimension (LDA,N)
            On entry, the M-by-N matrix A.
            On exit, the elements on and below the diagonal of the array
            contain the m-by-min(m,n) lower trapezoidal matrix L (L is
            lower triangular if m <= n); the elements above the diagonal,
            with the array TAU, represent the orthogonal matrix Q as a
            product of elementary reflectors (see Further Details).

            Higher performance is achieved if A is in pinned memory, e.g.
            allocated using magma_malloc_pinned.

    LDA     (input) INTEGER
            The leading dimension of the array A.  LDA >= max(1,M).

    TAU     (output) REAL array, dimension (min(M,N))
            The scalar factors of the elementary reflectors (see Further
            Details).

    WORK    (workspace/output) REAL array, dimension (MAX(1,LWORK))
            On exit, if INFO = 0, WORK(1) returns the optimal LWORK.

            Higher performance is achieved if WORK is in pinned memory, e.g.
            allocated using magma_malloc_pinned.

    LWORK   (input) INTEGER
            The dimension of the array WORK.  LWORK >= max(1,M).
            For optimum performance LWORK >= M*NB, where NB is the
            optimal blocksize.

            If LWORK = -1, then a workspace query is assumed; the routine
            only calculates the optimal size of the WORK array, returns
            this value as the first entry of the WORK array, and no error
            message related to LWORK is issued.

    INFO    (output) INTEGER
            = 0:  successful exit
            < 0:  if INFO = -i, the i-th argument had an illegal value
                  if INFO = -10 internal GPU memory allocation failed.

    Further Details
    ===============
    The matrix Q is represented as a product of elementary reflectors

       Q = H(k) . . . H(2) H(1), where k = min(m,n).

    Each H(i) has the form

       H(i) = I - tau * v * v'

    where tau is a real scalar, and v is a real vector with
    v(1:i-1) = 0 and v(i) = 1; v(i+1:n) is stored on exit in A(i,i+1:n),
    and tau in TAU(i).
    =====================================================================    */

    #define  a_ref(a_1,a_2) ( a+(a_2)*(lda) + (a_1))

    float *dA, *dAT;
    float c_one = MAGMA_S_ONE;
    magma_int_t maxm, maxn, maxdim, nb;
    magma_int_t iinfo, ldda;
    int lquery;

    /* Function Body */
    *info = 0;
    nb = magma_get_sgelqf_nb(m);

    work[0] = MAGMA_S_MAKE( (float)(m*nb), 0 );
    lquery = (lwork == -1);
    if (m < 0) {
        *info = -1;
    } else if (n < 0) {
        *info = -2;
    } else if (lda < max(1,m)) {
        *info = -4;
    } else if (lwork < max(1,m) && ! lquery) {
        *info = -7;
    }
    if (*info != 0) {
        magma_xerbla( __func__, -(*info) );
        return *info;
    }
    else if (lquery) {
        return *info;
    }

    /*  Quick return if possible */
    if (min(m, n) == 0) {
        work[0] = c_one;
        return *info;
    }

    maxm = ((m + 31)/32)*32;
    maxn = ((n + 31)/32)*32;
    maxdim = max(maxm, maxn);

    if (maxdim*maxdim < 2*maxm*maxn)
        {
            ldda = maxdim;

            if (MAGMA_SUCCESS != magma_smalloc( &dA, maxdim*maxdim )) {
                *info = MAGMA_ERR_DEVICE_ALLOC;
                return *info;
            }

            magma_ssetmatrix( m, n, a, lda, dA, ldda );
            dAT = dA;
            magmablas_stranspose_inplace( ldda, dAT, ldda );
        }
    else
        {
            ldda = maxn;

            if (MAGMA_SUCCESS != magma_smalloc( &dA, 2*maxn*maxm )) {
                *info = MAGMA_ERR_DEVICE_ALLOC;
                return *info;
            }

            magma_ssetmatrix( m, n, a, lda, dA, maxm );

            dAT = dA + maxn * maxm;
            magmablas_stranspose2( dAT, ldda, dA, maxm, m, n );
        }

    magma_sgeqrf2_gpu(n, m, dAT, ldda, tau, &iinfo);

    if (maxdim*maxdim < 2*maxm*maxn) {
        magmablas_stranspose_inplace( ldda, dAT, ldda );
        magma_sgetmatrix( m, n, dA, ldda, a, lda );
    } else {
        magmablas_stranspose2( dA, maxm, dAT, ldda, n, m );
        magma_sgetmatrix( m, n, dA, maxm, a, lda );
    }

    magma_free( dA );

    return *info;
} /* magma_sgelqf */

#undef  a_ref
