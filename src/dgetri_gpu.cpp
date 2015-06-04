/*
    -- MAGMA (version 1.6.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date November 2014

       @generated from zgetri_gpu.cpp normal z -> d, Sat Nov 15 19:54:09 2014

*/
#include "common_magma.h"

#define PRECISION_d

// === Define what BLAS to use ============================================
#define magma_dtrsm magmablas_dtrsm
// === End defining what BLAS to use ======================================

/**
    Purpose
    -------
    DGETRI computes the inverse of a matrix using the LU factorization
    computed by DGETRF. This method inverts U and then computes inv(A) by
    solving the system inv(A)*L = inv(U) for inv(A).
    
    Note that it is generally both faster and more accurate to use DGESV,
    or DGETRF and DGETRS, to solve the system AX = B, rather than inverting
    the matrix and multiplying to form X = inv(A)*B. Only in special
    instances should an explicit inverse be computed with this routine.

    Arguments
    ---------
    @param[in]
    n       INTEGER
            The order of the matrix A.  N >= 0.

    @param[in,out]
    dA      DOUBLE_PRECISION array on the GPU, dimension (LDDA,N)
            On entry, the factors L and U from the factorization
            A = P*L*U as computed by DGETRF_GPU.
            On exit, if INFO = 0, the inverse of the original matrix A.

    @param[in]
    ldda    INTEGER
            The leading dimension of the array A.  LDDA >= max(1,N).

    @param[in]
    ipiv    INTEGER array, dimension (N)
            The pivot indices from DGETRF; for 1 <= i <= N, row i of the
            matrix was interchanged with row IPIV(i).

    @param[out]
    dwork   (workspace) DOUBLE_PRECISION array on the GPU, dimension (MAX(1,LWORK))
  
    @param[in]
    lwork   INTEGER
            The dimension of the array DWORK.  LWORK >= N*NB, where NB is
            the optimal blocksize returned by magma_get_dgetri_nb(n).
    \n
            Unlike LAPACK, this version does not currently support a
            workspace query, because the workspace is on the GPU.

    @param[out]
    info    INTEGER
      -     = 0:  successful exit
      -     < 0:  if INFO = -i, the i-th argument had an illegal value
      -     > 0:  if INFO = i, U(i,i) is exactly zero; the matrix is
                  singular and its cannot be computed.

    @ingroup magma_dgesv_comp
    ********************************************************************/
extern "C" magma_int_t
magma_dgetri_gpu(
    magma_int_t n,
    magmaDouble_ptr dA, magma_int_t ldda, magma_int_t *ipiv,
    magmaDouble_ptr dwork, magma_int_t lwork,
    magma_int_t *info )
{
    #define dA(i, j)  (dA + (i) + (j)*ldda)
    #define dL(i, j)  (dL + (i) + (j)*lddl)
    
    /* Local variables */
    double c_zero    = MAGMA_D_ZERO;
    double c_one     = MAGMA_D_ONE;
    double c_neg_one = MAGMA_D_NEG_ONE;
    magmaDouble_ptr dL = dwork;
    magma_int_t lddl = n;
    magma_int_t nb   = magma_get_dgetri_nb(n);
    magma_int_t j, jmax, jb, jp;
    
    *info = 0;
    if (n < 0)
        *info = -1;
    else if (ldda < max(1,n))
        *info = -3;
    else if ( lwork < n*nb )
        *info = -6;

    if (*info != 0) {
        magma_xerbla( __func__, -(*info) );
        return *info;
    }

    /* Quick return if possible */
    if ( n == 0 )
        return *info;
    
    /* Invert the triangular factor U */
    magma_dtrtri_gpu( MagmaUpper, MagmaNonUnit, n, dA, ldda, info );
    if ( *info != 0 )
        return *info;
    
    jmax = ((n-1) / nb)*nb;
    for( j = jmax; j >= 0; j -= nb ) {
        jb = min( nb, n-j );
        
        // copy current block column of A to work space dL
        // (only needs lower trapezoid, but we also copy upper triangle),
        // then zero the strictly lower trapezoid block column of A.
        magmablas_dlacpy( MagmaFull, n-j, jb,
                          dA(j,j), ldda,
                          dL(j,0), lddl );
        magmablas_dlaset( MagmaLower, n-j-1, jb, c_zero, c_zero, dA(j+1,j), ldda );
        
        // compute current block column of Ainv
        // Ainv(:, j:j+jb-1)
        //   = ( U(:, j:j+jb-1) - Ainv(:, j+jb:n) L(j+jb:n, j:j+jb-1) )
        //   * L(j:j+jb-1, j:j+jb-1)^{-1}
        // where L(:, j:j+jb-1) is stored in dL.
        if ( j+jb < n ) {
            magma_dgemm( MagmaNoTrans, MagmaNoTrans, n, jb, n-j-jb,
                         c_neg_one, dA(0,j+jb), ldda,
                                    dL(j+jb,0), lddl,
                         c_one,     dA(0,j),    ldda );
        }
        // TODO use magmablas work interface
        magma_dtrsm( MagmaRight, MagmaLower, MagmaNoTrans, MagmaUnit,
                     n, jb, c_one,
                     dL(j,0), lddl,
                     dA(0,j), ldda );
    }

    // Apply column interchanges
    for( j = n-2; j >= 0; --j ) {
        jp = ipiv[j] - 1;
        if ( jp != j ) {
            magmablas_dswap( n, dA(0,j), 1, dA(0,jp), 1 );
        }
    }
    
    return *info;
}
