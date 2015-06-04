/*
    -- MAGMA (version 1.6.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date November 2014

       @generated from zgeqr2x_gpu-v4.cu normal z -> s, Sat Nov 15 19:53:59 2014

*/
#include "common_magma.h"
#include "commonblas_s.h"

// 512 is maximum number of threads for CUDA capability 1.x
#define BLOCK_SIZE 512


/**
    Purpose
    -------
    SGEQR2 computes a QR factorization of a real m by n matrix A:
    A = Q * R.

    This expert routine requires two more arguments than the standard
    sgeqr2, namely, dT and ddA, explained below. The storage for A is
    also not as in the LAPACK's sgeqr2 routine (see below).

    The first is used to output the triangular
    n x n factor T of the block reflector used in the factorization.
    The second holds the diagonal nxn blocks of A, i.e., the diagonal
    submatrices of R. This routine implements the left looking QR.

    This version adds internal blocking.

    Arguments
    ---------
    @param[in]
    m       INTEGER
            The number of rows of the matrix A.  M >= 0.

    @param[in]
    n       INTEGER
            The number of columns of the matrix A.  N >= 0.

    @param[in,out]
    dA      REAL array, dimension (LDA,N)
            On entry, the m by n matrix A.
            On exit, the unitary matrix Q as a
            product of elementary reflectors (see Further Details).
    \n
            the elements on and above the diagonal of the array
            contain the min(m,n) by n upper trapezoidal matrix R (R is
            upper triangular if m >= n); the elements below the diagonal,
            with the array TAU, represent the unitary matrix Q as a
            product of elementary reflectors (see Further Details).

    @param[in]
    ldda    INTEGER
            The leading dimension of the array A.  LDA >= max(1,M).

    @param[out]
    dtau    REAL array, dimension (min(M,N))
            The scalar factors of the elementary reflectors (see Further
            Details).

    @param[out]
    dT      REAL array, dimension N x N.
            Stores the triangular N x N factor T of the block reflector
            used in the factorization. The lower triangular part is 0.

    @param[out]
    ddA     REAL array, dimension N x N.
            Stores the elements of the upper N x N diagonal block of A.
            LAPACK stores this array in A. There are 0s below the diagonal.

    @param
    dwork   (workspace) DOUBLE_PRECISION array, dimension (3 N)

    @param[out]
    info    INTEGER
      -     = 0: successful exit
      -     < 0: if INFO = -i, the i-th argument had an illegal value

    Further Details
    ---------------
    The matrix Q is represented as a product of elementary reflectors

        Q = H(1) H(2) . . . H(k), where k = min(m,n).

    Each H(i) has the form

        H(i) = I - tau * v * v**H

    where tau is a real scalar, and v is a real vector with
    v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i),
    and tau in TAU(i).

    @ingroup magma_sgeqrf_comp
    ********************************************************************/
extern "C" magma_int_t
magma_sgeqr2x4_gpu(
    magma_int_t m, magma_int_t n,
    magmaFloat_ptr dA, magma_int_t ldda,
    magmaFloat_ptr dtau,
    magmaFloat_ptr dT,
    magmaFloat_ptr ddA,
    magmaFloat_ptr dwork,
    magma_queue_t queue,
    magma_int_t *info)
{
    #define dA(i_,j_) (dA + (j_)*(ldda) + (i_))
    #define dT(i_,j_) (dT + (j_)*(k)    + (i_))
    #define BS 32

    magma_int_t i, k;

    float *dnorm = (float *)dwork;
    float *work = (float *)(dwork+2*n);

    magma_queue_t cstream;
    magmablasGetKernelStream(&cstream);
    magmablasSetKernelStream(queue);

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

    /* Compute the norms of the trailing columns */
    k = min(m,n);
    magmablas_snrm2_cols(m, k, dA(0,0), ldda, dnorm);

    for (magma_int_t b=0; b < k; b += BS) {
        for (i = b; i < min(k, b+BS); ++i) {

            /*   Apply H**H to A(:,i) from the left */
            if (i-b > 0) {
                /* Compute the (i-1)th column of T */
                if ( i-1 > 0 ) {
                    magma_sgemv_kernel3<<< i-1, BLOCK_SIZE, 0, magma_stream >>>
                        ( m-i+1, dA(i-1,0), ldda, dA(i-1, i-1), work, dtau+i-1);
                    magma_strmv_kernel2<<< i-1, i-1, 0, magma_stream >>>
                        ( dT(0,0), k, work, dT(0,i-1), dtau+i-1);
                }

                /* dwork = V**H c */
                magma_sgemv_kernel1<<< i-b, BLOCK_SIZE, 0, magma_stream >>>
                    (m-b, dA(b, b),  ldda, dA(b,i), work);

                /* dwork = T**H work */
                magma_strmv_tkernel<<< i-b, i-b, 0, magma_stream >>>
                    (dT(b,b), k, work, work+i-b);

                /* c = c - V work */
                if ( m-b > 0 ) {
                    dim3  blocks3( (m-b + BLOCK_SIZE-1) / BLOCK_SIZE );
                    dim3 threads3( BLOCK_SIZE );
                    magma_sgemv_kernel2<<< blocks3, threads3, 0, magma_stream >>>
                        (m-b, i-b, dA(b,b), ldda,  work+i-b, dA(b, i));
                }
            }

            /*   Adjust the dnorm[i] to hold the norm of A(i:m,i) */
            if ( i > 0 ) {
                magma_snrm2_adjust_kernel<<< 1, i, 0, magma_stream >>>(dnorm+i, dA(0, i));
            }

            /*  Generate elementary reflector H(i) to annihilate A(i+1:m,i)
                1. 1 is not yet put on the diagonal of A
                2. Elements above the diagonal are copied in ddA and
                   the ones in A are set to zero
                3. update T */
            magma_slarfgx_gpu(m-i, dA(i, i), dA(min(i+1,m),i), dtau+i,
                              dnorm+i, ddA + i + i*n, i);

            if (i == 0) {
                float tt = MAGMA_S_ONE;
                magmablas_slacpy(MagmaUpperLower, 1, 1, dtau, 1, dT(0,0), 1);
                magma_ssetmatrix_async(1, 1, &tt, 1, dA(i, i), 1, magma_stream);
            }
        }
        if ( i-1 > 0 ) {
            magma_sgemv_kernel3<<< i-1, BLOCK_SIZE, 0, magma_stream >>>
                ( m-i+1, dA(i-1,0), ldda, dA(i-1, i-1), work, dtau+i-1);
            magma_strmv_kernel2<<< i-1, i-1, 0, magma_stream >>>
                ( dT(0,0), k, work, dT(0,i-1), dtau+i-1);
        }

        /* Apply the transformations to the trailing matrix. */
        //magma_slarfb2_gpu( MagmaLeft, MagmaConjTrans, MagmaForward, MagmaColumnwise,
        magma_slarfb2_gpu(
                           m-b, k-i, BS,
                           dA(b, b), ldda, dT+b+b*k, k,
                           dA(b, i), ldda, work, k-i);
    }

    magmablasSetKernelStream(cstream);

    return *info;
} /* magma_sgeqr2 */
