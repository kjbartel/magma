/*
    -- MAGMA (version 1.5.0-beta2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2014

       @generated from zungqr2.cpp normal z -> c, Fri May 30 10:41:01 2014

       @author Stan Tomov
       @author Mark Gates
*/
#include "common_magma.h"

/**
    Purpose
    -------
    CUNGQR generates an M-by-N COMPLEX matrix Q with orthonormal columns,
    which is defined as the first N columns of a product of K elementary
    reflectors of order M

          Q  =  H(1) H(2) . . . H(k)

    as returned by CGEQRF.

    This version recomputes the T matrices on the CPU and sends them to the GPU.

    Arguments
    ---------
    @param[in]
    m       INTEGER
            The number of rows of the matrix Q. M >= 0.

    @param[in]
    n       INTEGER
            The number of columns of the matrix Q. M >= N >= 0.

    @param[in]
    k       INTEGER
            The number of elementary reflectors whose product defines the
            matrix Q. N >= K >= 0.

    @param[in,out]
    A       COMPLEX array A, dimension (LDDA,N).
            On entry, the i-th column must contain the vector
            which defines the elementary reflector H(i), for
            i = 1,2,...,k, as returned by CGEQRF_GPU in the
            first k columns of its array argument A.
            On exit, the M-by-N matrix Q.

    @param[in]
    lda     INTEGER
            The first dimension of the array A. LDA >= max(1,M).

    @param[in]
    tau     COMPLEX array, dimension (K)
            TAU(i) must contain the scalar factor of the elementary
            reflector H(i), as returned by CGEQRF_GPU.

    @param[out]
    info    INTEGER
      -     = 0:  successful exit
      -     < 0:  if INFO = -i, the i-th argument has an illegal value

    @ingroup magma_cgeqrf_comp
    ********************************************************************/
extern "C" magma_int_t
magma_cungqr2(magma_int_t m, magma_int_t n, magma_int_t k,
              magmaFloatComplex *A, magma_int_t lda,
              magmaFloatComplex *tau,
              magma_int_t *info)
{
#define  A(i,j) ( A + (i) + (j)*lda )
#define dA(i,j) (dA + (i) + (j)*ldda)

    magmaFloatComplex c_zero = MAGMA_C_ZERO;
    magmaFloatComplex c_one  = MAGMA_C_ONE;

    magma_int_t nb = magma_get_cgeqrf_nb(min(m, n));

    magma_int_t  m_kk, n_kk, k_kk, mi;
    magma_int_t lwork, ldda;
    magma_int_t i, ib, ki, kk;  //, iinfo;
    magma_int_t lddwork;
    magmaFloatComplex *dA, *dV, *dW, *dT, *T;
    magmaFloatComplex *work;

    *info = 0;
    if (m < 0) {
        *info = -1;
    } else if ((n < 0) || (n > m)) {
        *info = -2;
    } else if ((k < 0) || (k > n)) {
        *info = -3;
    } else if (lda < max(1,m)) {
        *info = -5;
    }
    if (*info != 0) {
        magma_xerbla( __func__, -(*info) );
        return *info;
    }

    if (n <= 0) {
        return *info;
    }

    // first kk columns are handled by blocked method.
    // ki is start of 2nd-to-last block
    if ((nb > 1) && (nb < k)) {
        ki = (k - nb - 1) / nb * nb;
        kk = min(k, ki + nb);
    } else {
        ki = 0;
        kk = 0;
    }

    // Allocate GPU work space
    // ldda*n     for matrix dA
    // ldda*nb    for dV
    // lddwork*nb for dW larfb workspace
    ldda    = ((m + 31) / 32) * 32;
    lddwork = ((n + 31) / 32) * 32;
    if (MAGMA_SUCCESS != magma_cmalloc( &dA, ldda*n + ldda*nb + lddwork*nb + nb*nb)) {
        *info = MAGMA_ERR_DEVICE_ALLOC;
        return *info;
    }
    dV = dA + ldda*n;
    dW = dA + ldda*n + ldda*nb;
    dT = dA + ldda*n + ldda*nb + lddwork*nb;

    // Allocate CPU work space
    lwork = (n+m+nb) * nb;
    magma_cmalloc_cpu( &work, lwork );

    T = work;

    if (work == NULL) {
        magma_free( dA );
        magma_free_cpu( work );
        *info = MAGMA_ERR_HOST_ALLOC;
        return *info;
    }
    magmaFloatComplex *V = work + (n+nb)*nb;

    magma_queue_t stream;
    magma_queue_create( &stream );

    // Use unblocked code for the last or only block.
    if (kk < n) {
        m_kk = m - kk;
        n_kk = n - kk;
        k_kk = k - kk;
        /*
            lapackf77_cungqr( &m_kk, &n_kk, &k_kk,
                              A(kk, kk), &lda,
                              &tau[kk], work, &lwork, &iinfo );
        */
        lapackf77_clacpy( MagmaUpperLowerStr, &m_kk, &k_kk, A(kk,kk), &lda, V, &m_kk);
        lapackf77_claset( MagmaUpperLowerStr, &m_kk, &n_kk, &c_zero, &c_one, A(kk, kk), &lda );

        lapackf77_clarft( MagmaForwardStr, MagmaColumnwiseStr,
                          &m_kk, &k_kk,
                          V, &m_kk, &tau[kk], work, &k_kk);
        lapackf77_clarfb( MagmaLeftStr, MagmaNoTransStr, MagmaForwardStr, MagmaColumnwiseStr,
                          &m_kk, &n_kk, &k_kk,
                          V, &m_kk, work, &k_kk, A(kk, kk), &lda, work+k_kk*k_kk, &n_kk );
        
        if (kk > 0) {
            magma_csetmatrix( m_kk, n_kk,
                              A(kk, kk),  lda,
                              dA(kk, kk), ldda );
        
            // Set A(1:kk,kk+1:n) to zero.
            magmablas_claset( MagmaUpperLower, kk, n - kk, dA(0, kk), ldda );
        }
    }

    if (kk > 0) {
        // Use blocked code
        // stream: set Aii (V) --> laset --> laset --> larfb --> [next]
        // CPU has no computation
        magmablasSetKernelStream( stream );
        
        for (i = ki; i >= 0; i -= nb) {
            ib = min(nb, k - i);

            // Send current panel to the GPU
            mi = m - i;
            lapackf77_claset( "Upper", &ib, &ib, &c_zero, &c_one, A(i, i), &lda );
            magma_csetmatrix_async( mi, ib,
                                    A(i, i), lda,
                                    dV,      ldda, stream );
            lapackf77_clarft( MagmaForwardStr, MagmaColumnwiseStr,
                              &mi, &ib,
                              A(i,i), &lda, &tau[i], T, &nb);
            magma_csetmatrix_async( ib, ib,
                                    T, nb,
                                    dT, nb, stream );

            // set panel to identity
            magmablas_claset( MagmaUpperLower, i, ib, dA(0, i), ldda );
            magmablas_claset_identity( mi, ib, dA(i, i), ldda );
            
            magma_queue_sync( stream );
            if (i < n) {
                // Apply H to A(i:m,i:n) from the left
                magma_clarfb_gpu( MagmaLeft, MagmaNoTrans, MagmaForward, MagmaColumnwise,
                                  mi, n-i, ib,
                                  dV,       ldda, dT, nb,
                                  dA(i, i), ldda, dW, lddwork );
            }
        }
    
        // copy result back to CPU
        magma_cgetmatrix( m, n,
                          dA(0, 0), ldda, A(0, 0), lda);
    }

    magmablasSetKernelStream( NULL );
    magma_queue_destroy( stream );
    magma_free( dA );
    magma_free_cpu( work );

    return *info;
} /* magma_cungqr */
