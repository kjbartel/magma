/*
    -- MAGMA (version 1.2.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       May 2012

       @generated s Tue May 15 18:17:33 2012

*/
#include "common_magma.h"

extern "C" magma_int_t
magma_sorgqr(magma_int_t m, magma_int_t n, magma_int_t k,
             float *a, magma_int_t lda,
             float *tau, float *dT,
             magma_int_t nb, magma_int_t *info)
{
/*  -- MAGMA (version 1.2.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       May 2012

    Purpose
    =======
    SORGQR generates an M-by-N REAL matrix Q with orthonormal columns,
    which is defined as the first N columns of a product of K elementary
    reflectors of order M

          Q  =  H(1) H(2) . . . H(k)

    as returned by SGEQRF.

    Arguments
    =========
    M       (input) INTEGER
            The number of rows of the matrix Q. M >= 0.

    N       (input) INTEGER
            The number of columns of the matrix Q. M >= N >= 0.

    K       (input) INTEGER
            The number of elementary reflectors whose product defines the
            matrix Q. N >= K >= 0.

    A       (input/output) REAL array A, dimension (LDDA,N). 
            On entry, the i-th column must contain the vector
            which defines the elementary reflector H(i), for
            i = 1,2,...,k, as returned by SGEQRF_GPU in the 
            first k columns of its array argument A.
            On exit, the M-by-N matrix Q.

    LDA     (input) INTEGER
            The first dimension of the array A. LDA >= max(1,M).

    TAU     (input) REAL array, dimension (K)
            TAU(i) must contain the scalar factor of the elementary
            reflector H(i), as returned by SGEQRF_GPU.

    DT      (input) REAL array on the GPU device.
            DT contains the T matrices used in blocking the elementary
            reflectors H(i), e.g., this can be the 6th argument of 
            magma_sgeqrf_gpu.

    NB      (input) INTEGER
            This is the block size used in SGEQRF_GPU, and correspondingly
            the size of the T matrices, used in the factorization, and
            stored in DT.

    INFO    (output) INTEGER
            = 0:  successful exit
            < 0:  if INFO = -i, the i-th argument has an illegal value
    =====================================================================    */

    #define  a_ref(i,j)     ( a + (j)*lda  + (i))
    #define da_ref(i,j)     (da + (j)*ldda + (i))
    #define t_ref(a_1)      (dT+(a_1)*nb)

    magma_int_t  i__1, i__2, i__3;
    magma_int_t lwork, ldda;
    static magma_int_t i, ib, ki, kk, iinfo;
    magma_int_t lddwork = min(m, n);
    float *da, *work, *dwork;
    static cudaStream_t stream;

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

    if (n <= 0)
      return *info;

    /* Allocate GPU work space */
    ldda = ((m+31)/32)*32;
    lddwork = ((lddwork+31)/32)*32;
    if (MAGMA_SUCCESS != magma_smalloc( &da, (n)*ldda + nb*lddwork )) {
        *info = MAGMA_ERR_DEVICE_ALLOC;
        return *info;
    }
    dwork = da + (n)*ldda;

    /* Allocate CPU work space */
    lwork = n * nb;
    work = (float *)malloc(lwork*sizeof(float));
    if( work == NULL ) {
        magma_free( da );
        *info = MAGMA_ERR_HOST_ALLOC;
        return *info;
    }

    magma_queue_create( &stream );

    if ( (nb > 1) && (nb < k) )
      {
        /*  Use blocked code after the last block.
            The first kk columns are handled by the block method. */
        ki = (k - nb - 1) / nb * nb;
        kk = min(k, ki + nb);

        /* Set A(1:kk,kk+1:n) to zero. */
        magmablas_slaset(MagmaUpperLower, kk, n-kk, da_ref(0,kk), ldda);
      }
    else
      kk = 0;

    /* Use unblocked code for the last or only block. */
    if (kk < n)
      {
        i__1 = m - kk;
        i__2 = n - kk;
        i__3 = k - kk;
        lapackf77_sorgqr(&i__1, &i__2, &i__3, 
                         a_ref(kk, kk), &lda,
                         &tau[kk], work, &lwork, &iinfo);

        magma_ssetmatrix( i__1, i__2,
                          a_ref(kk, kk),  lda,
                          da_ref(kk, kk), ldda );
      }

    if (kk > 0)
      {
        /* Use blocked code */
        for (i = ki; i >= 0; i-=nb)
          {
            ib = min(nb, k - i);

            /* Send the current panel to the GPU */
            i__2 = m - i;
            spanel_to_q(MagmaUpper, ib, a_ref(i,i), lda, work);
            magma_ssetmatrix( i__2, ib, a_ref(i, i), lda, da_ref(i, i), ldda );
                             
            if (i + ib < n)
              {
                /* Apply H to A(i:m,i+ib:n) from the left */
                i__3 = n - i - ib;
                magma_slarfb_gpu( MagmaLeft, MagmaNoTrans, MagmaForward, MagmaColumnwise,
                                  i__2, i__3, ib,
                                  da_ref(i, i   ), ldda, t_ref(i),      nb,
                                  da_ref(i, i+ib), ldda,    dwork, lddwork);
              }

            /* Apply H to rows i:m of current block on the CPU */
            lapackf77_sorgqr(&i__2, &ib, &ib, 
                             a_ref(i, i), &lda, 
                             &tau[i], work, &lwork, &iinfo);
            magma_ssetmatrix_async( i__2, ib,
                                    a_ref(i,i),  lda,
                                    da_ref(i,i), ldda, stream );

            /* Set rows 1:i-1 of current block to zero */
            i__2 = i + ib;
            magmablas_slaset(MagmaUpperLower, i, i__2 - i, da_ref(0,i), ldda);
          }
      }

    magma_sgetmatrix( m, n, da_ref(0, 0), ldda, a_ref(0, 0), lda );

    
    magma_queue_destroy( stream );
    magma_free( da );
    free(work);

    return *info;
} /* magma_sorgqr */

#undef da_ref
#undef a_ref
#undef t_ref
