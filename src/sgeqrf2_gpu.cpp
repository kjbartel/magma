/*
    -- MAGMA (version 1.4.1-beta2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       December 2013

       @author Stan Tomov
       @generated s Mon Dec  9 17:05:22 2013

*/
#include "common_magma.h"

extern "C" magma_int_t
magma_sgeqrf2_gpu( magma_int_t m, magma_int_t n,
                   float *dA, magma_int_t ldda,
                   float *tau,
                   magma_int_t *info )
{
/*  -- MAGMA (version 1.4.1-beta2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       December 2013

    Purpose
    =======
    SGEQRF computes a QR factorization of a real M-by-N matrix A:
    A = Q * R.
    
    This version has LAPACK-complaint arguments. 
    If the current stream is NULL, this version replaces it with user defined
    stream to overlap computation with communication.    

    Other versions (magma_sgeqrf_gpu and magma_sgeqrf3_gpu) store the
    intermediate T matrices.

    Arguments
    =========
    M       (input) INTEGER
            The number of rows of the matrix A.  M >= 0.

    N       (input) INTEGER
            The number of columns of the matrix A.  N >= 0.

    dA      (input/output) REAL array on the GPU, dimension (LDDA,N)
            On entry, the M-by-N matrix A.
            On exit, the elements on and above the diagonal of the array
            contain the min(M,N)-by-N upper trapezoidal matrix R (R is
            upper triangular if m >= n); the elements below the diagonal,
            with the array TAU, represent the orthogonal matrix Q as a
            product of min(m,n) elementary reflectors (see Further
            Details).

    LDDA    (input) INTEGER
            The leading dimension of the array dA.  LDDA >= max(1,M).
            To benefit from coalescent memory accesses LDDA must be
            dividable by 16.

    TAU     (output) REAL array, dimension (min(M,N))
            The scalar factors of the elementary reflectors (see Further
            Details).

    INFO    (output) INTEGER
            = 0:  successful exit
            < 0:  if INFO = -i, the i-th argument had an illegal value
                  or another error occured, such as memory allocation failed.

    Further Details
    ===============
    The matrix Q is represented as a product of elementary reflectors

       Q = H(1) H(2) . . . H(k), where k = min(m,n).

    Each H(i) has the form

       H(i) = I - tau * v * v'

    where tau is a real scalar, and v is a real vector with
    v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i),
    and tau in TAU(i).
    =====================================================================    */

    #define dA(a_1,a_2)    ( dA+(a_2)*(ldda) + (a_1))
    #define work_ref(a_1)  ( work + (a_1))
    #define hwork          ( work + (nb)*(m))

    float *dwork;
    float *work;
    magma_int_t i, k, ldwork, lddwork, old_i, old_ib, rows;
    magma_int_t nbmin, nx, ib, nb;
    magma_int_t lhwork, lwork;

    /* Function Body */
    *info = 0;
    if (m < 0) {
        *info = -1;
    } else if (n < 0) {
        *info = -2;
    } else if (ldda < max(1,m)) {
        *info = -4;
    }
    if (*info != 0) {
        magma_xerbla( __func__, -(*info) );
        return *info;
    }

    k = min(m,n);
    if (k == 0)
        return *info;

    nb = magma_get_sgeqrf_nb(m);

    lwork  = (m+n) * nb;
    lhwork = lwork - (m)*nb;

    if (MAGMA_SUCCESS != magma_smalloc( &dwork, (n)*nb )) {
        *info = MAGMA_ERR_DEVICE_ALLOC;
        return *info;
    }

    if (MAGMA_SUCCESS != magma_smalloc_pinned( &work, lwork )) {
        magma_free( dwork );
        *info = MAGMA_ERR_HOST_ALLOC;
        return *info;
    }

    /* Define user stream if current stream is NULL */
    magma_queue_t stream[2], current_stream;
    magmablasGetKernelStream(&current_stream);

    magma_queue_create( &stream[0] );
    if (current_stream == NULL) {
        magma_queue_create( &stream[1] );
        magmablasSetKernelStream(stream[1]);
    }
    else
        stream[1] = current_stream;

    nbmin = 2;
    nx    = nb;
    ldwork = m;
    lddwork= n;

    if (nb >= nbmin && nb < k && nx < k) {
        /* Use blocked code initially */
        old_i = 0; old_ib = nb;
        for (i = 0; i < k-nx; i += nb) {
            ib = min(k-i, nb);
            rows = m -i;

            /* download i-th panel */
            magma_queue_sync( stream[1] );
            magma_sgetmatrix_async( rows, ib,
                                    dA(i,i),       ldda,
                                    work_ref(i), ldwork, stream[0] );
            if (i>0){
                /* Apply H' to A(i:m,i+2*ib:n) from the left */
                magma_slarfb_gpu( MagmaLeft, MagmaTrans, MagmaForward, MagmaColumnwise,
                                  m-old_i, n-old_i-2*old_ib, old_ib,
                                  dA(old_i, old_i         ), ldda, dwork,        lddwork,
                                  dA(old_i, old_i+2*old_ib), ldda, dwork+old_ib, lddwork);

                magma_ssetmatrix_async( old_ib, old_ib,
                                        work_ref(old_i),  ldwork,
                                        dA(old_i, old_i), ldda, stream[1] );
            }

            magma_queue_sync( stream[0] );
            lapackf77_sgeqrf(&rows, &ib, work_ref(i), &ldwork, tau+i, hwork, &lhwork, info);
            /* Form the triangular factor of the block reflector
               H = H(i) H(i+1) . . . H(i+ib-1) */
            lapackf77_slarft( MagmaForwardStr, MagmaColumnwiseStr,
                              &rows, &ib,
                              work_ref(i), &ldwork, tau+i, hwork, &ib);

            spanel_to_q( MagmaUpper, ib, work_ref(i), ldwork, hwork+ib*ib );

            /* download the i-th V matrix */
            magma_ssetmatrix_async( rows, ib, work_ref(i), ldwork, dA(i,i), ldda, stream[0] );

            /* download the T matrix */
            magma_queue_sync( stream[1] );
            magma_ssetmatrix_async( ib, ib, hwork, ib, dwork, lddwork, stream[0] );
            magma_queue_sync( stream[0] );

            if (i + ib < n) {

                if (i+nb < k-nx) {
                    /* Apply H' to A(i:m,i+ib:i+2*ib) from the left */
                    magma_slarfb_gpu( MagmaLeft, MagmaTrans, MagmaForward, MagmaColumnwise,
                                      rows, ib, ib,
                                      dA(i, i   ), ldda, dwork,    lddwork,
                                      dA(i, i+ib), ldda, dwork+ib, lddwork);
                    sq_to_panel( MagmaUpper, ib, work_ref(i), ldwork, hwork+ib*ib );
                }
                else {
                    magma_slarfb_gpu( MagmaLeft, MagmaTrans, MagmaForward, MagmaColumnwise,
                                      rows, n-i-ib, ib,
                                      dA(i, i   ), ldda, dwork,    lddwork,
                                      dA(i, i+ib), ldda, dwork+ib, lddwork);
                    sq_to_panel( MagmaUpper, ib, work_ref(i), ldwork, hwork+ib*ib );
                    magma_ssetmatrix_async( ib, ib,
                                            work_ref(i), ldwork,
                                            dA(i,i),     ldda, stream[1] );
                }
                old_i  = i;
                old_ib = ib;
            }
        }
    } else {
        i = 0;
    }
    magma_free( dwork );

    /* Use unblocked code to factor the last or only block. */
    if (i < k) {
        ib   = n-i;
        rows = m-i;
        magma_sgetmatrix_async( rows, ib, dA(i, i), ldda, work, rows, stream[1] );
        magma_queue_sync( stream[1] );
        lhwork = lwork - rows*ib;
        lapackf77_sgeqrf(&rows, &ib, work, &rows, tau+i, work+ib*rows, &lhwork, info);
        
        magma_ssetmatrix_async( rows, ib, work, rows, dA(i, i), ldda, stream[1] );
    }

    magma_free_pinned( work );

    magma_queue_destroy( stream[0] );
    if (current_stream == NULL) {
      magma_queue_destroy( stream[1] );
      magmablasSetKernelStream(NULL);
    }

    return *info;
}   /* magma_sgeqrf2_gpu */
