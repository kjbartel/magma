/*
    -- MAGMA (version 1.4.1-beta2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       December 2013

       @generated d Mon Dec  9 17:05:22 2013

*/
#include "common_magma.h"

extern "C" magma_int_t
magma_dlarfb2_gpu( magma_int_t m, magma_int_t n, magma_int_t k,
                   const double *dV,    magma_int_t ldv,
                   const double *dT,    magma_int_t ldt,
                   double *dC,          magma_int_t ldc,
                   double *dwork,       magma_int_t ldwork )
{
    double c_zero    = MAGMA_D_ZERO;
    double c_one     = MAGMA_D_ONE;
    double c_neg_one = MAGMA_D_NEG_ONE;

    if (m <= 0 || n <= 0)
        return MAGMA_SUCCESS;

    // W = C^H V
    // magma_dgemm( MagmaTrans, MagmaNoTrans,
    magmablas_dgemm_reduce(
                           n, k, m,
                           c_one,  dC,    ldc,
                           dV,    ldv,
                           c_zero, dwork, ldwork);

    // W = W T^H = C^H V T^H
    magma_dtrmm( MagmaRight, MagmaUpper, MagmaNoTrans, MagmaNonUnit,
                 n, k,
                 c_one, dT,    ldt,
                 dwork, ldwork);

    // C = C - V W^H = C - V T V^H C = (I - V T V^H) C = H C
    magma_dgemm( MagmaNoTrans, MagmaTrans,
                 m, n, k,
                 c_neg_one, dV,    ldv,
                 dwork, ldwork,
                 c_one,     dC,    ldc);
    
    return MAGMA_SUCCESS;
}

//////////////////////////////////////////////////////////////////////////////

extern "C" magma_int_t
magma_dgeqr2x3_gpu(magma_int_t *m, magma_int_t *n, double *dA,
                   magma_int_t *ldda, double *dtau,
                   double *dT, double *ddA,
                   double *dwork, magma_int_t *info)
{
/*  -- MAGMA (version 1.4.1-beta2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       December 2013

    Purpose
    =======
    DGEQR2 computes a QR factorization of a real m by n matrix A:
    A = Q * R.

    This expert routine requires two more arguments than the standard
    dgeqr2, namely, dT and ddA, explained below. The storage for A is
    also not as in the LAPACK's dgeqr2 routine (see below).

    The first is used to output the triangular
    n x n factor T of the block reflector used in the factorization.
    The second holds the diagonal nxn blocks of A, i.e., the diagonal
    submatrices of R. This routine implements the left looking QR.

    This version adds internal blocking.

    Arguments
    =========
    M       (input) INTEGER
            The number of rows of the matrix A.  M >= 0.

    N       (input) INTEGER
            The number of columns of the matrix A.  N >= 0.

    A       (input/output) DOUBLE_PRECISION array, dimension (LDA,N)
            On entry, the m by n matrix A.
            On exit, the unitary matrix Q as a
            product of elementary reflectors (see Further Details).

            the elements on and above the diagonal of the array
            contain the min(m,n) by n upper trapezoidal matrix R (R is
            upper triangular if m >= n); the elements below the diagonal,
            with the array TAU, represent the unitary matrix Q as a
            product of elementary reflectors (see Further Details).

    LDA     (input) INTEGER
            The leading dimension of the array A.  LDA >= max(1,M).

    TAU     (output) DOUBLE_PRECISION array, dimension (min(M,N))
            The scalar factors of the elementary reflectors (see Further
            Details).

    dT      (output) DOUBLE_PRECISION array, dimension N x N.
            Stores the triangular N x N factor T of the block reflector
            used in the factorization. The lower triangular part is 0.

    ddA     (output) DOUBLE_PRECISION array, dimension N x N.
            Stores the elements of the upper N x N diagonal block of A.
            LAPACK stores this array in A. There are 0s below the diagonal.

    RWORK   (workspace) DOUBLE_PRECISION array, dimension (3 N)

    INFO    (output) INTEGER
            = 0: successful exit
            < 0: if INFO = -i, the i-th argument had an illegal value

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

    #define da_ref(a_1,a_2) ( dA+(a_2)*(*ldda) + (a_1))
    #define BLOCK_SIZE 32

    magma_int_t i, k;

    double *dnorm = dwork;
    double *work = (double *)(dwork+2*(*n));

    *info = 0;
    if (*m < 0) {
        *info = -1;
    } else if (*n < 0) {
        *info = -2;
    } else if (*ldda < max(1,*m)) {
        *info = -4;
    }
    if (*info != 0) {
        magma_xerbla( __func__, -(*info) );
        return *info;
    }

    /* Compute the norms of the trailing columns */
    k = min(*m,*n);
    magmablas_dnrm2_cols(*m, k, da_ref(0,0), *ldda, dnorm);

    for (int b=0; b < k; b += BLOCK_SIZE) {
        for (i = b; i < min(k, b+BLOCK_SIZE); ++i) {

            /*   Apply H' to A(:,i) from the left                           */
            if ( i-b > 0)
                magma_dlarfbx_gpu(*m-b, i-b, da_ref(b, b), *ldda,
                                  dT+b+b*k, k, da_ref(b, i), work);

            /*   Adjust the dnorm[i] to hold the norm of A(i:m,i)           */
            if ( i > 0 )
                magmablas_dnrm2_adjust(i, dnorm+i, da_ref(0, i));
            
            /*  Generate elementary reflector H(i) to annihilate A(i+1:m,i)
                1. 1 is not yet put on the diagonal of A
                2. Elements above the diagonal are copied in ddA and
                   the ones in A are set to zero
                3. update T                                                 */
            magma_dlarfgtx_gpu(*m-i, da_ref(i, i), da_ref(min(i+1,*m), i), dtau+i,
                               dnorm+i, ddA + i + i*(*n), i,
                               da_ref(i,0), *ldda,  dT, k, work);
        }
        
        /* Apply the transformations to the trailing matrix. */
        //magma_dlarfb2_gpu( MagmaLeft, MagmaTrans, MagmaForward, MagmaColumnwise,
        magma_dlarfb2_gpu(
                           *m-b, k-i, BLOCK_SIZE,
                           da_ref(b, b), *ldda, dT+b+b*k, k,
                           da_ref(b, i), *ldda, work, k-i);
    }

    return *info;
} /* magma_dgeqr2 */
