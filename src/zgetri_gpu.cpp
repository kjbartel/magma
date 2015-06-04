/*
    -- MAGMA (version 1.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2011

       @precisions normal z -> s d c

*/
#include "common_magma.h"

// === Define what BLAS to use ============================================
#define PRECISION_z
#if (defined(PRECISION_s) || defined(PRECISION_d))
  #define cublasZgemm magmablas_zgemm
  #define cublasZtrsm magmablas_ztrsm
#endif

#if (GPUSHMEM >= 200)
#if (defined(PRECISION_s))
    #undef  cublasSgemm
    #define cublasSgemm magmablas_sgemm_fermi80
#endif
#endif
// === End defining what BLAS to use ======================================

extern "C" magma_int_t
magma_zgetri_gpu( magma_int_t n, cuDoubleComplex *dA, magma_int_t lda,
                  magma_int_t *ipiv, cuDoubleComplex *dwork, magma_int_t lwork,
                  magma_int_t *info )
{
/*  -- MAGMA (version 1.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2011

    Purpose
    =======

        ZGETRI computes the inverse of a matrix using the LU factorization
        computed by ZGETRF. This method inverts U and then computes inv(A) by
        solving the system inv(A)*L = inv(U) for inv(A).
        
        Note that it is generally both faster and more accurate to use ZGESV,
        or ZGETRF and ZGETRS, to solve the system AX = B, rather than inverting
        the matrix and multiplying to form X = inv(A)*B. Only in special
        instances should an explicit inverse be computed with this routine.

    Arguments
    =========

        N       (input) INTEGER
                The order of the matrix A.  N >= 0.

        dA      (input/output) COMPLEX_16 array on the GPU, dimension (LDA,N)
                On entry, the factors L and U from the factorization
                A = P*L*U as computed by ZGETRF_GPU.
                On exit, if INFO = 0, the inverse of the original matrix A.

        LDA     (input) INTEGER
                The leading dimension of the array A.  LDA >= max(1,N).

        IPIV    (input) INTEGER array, dimension (N)
                The pivot indices from ZGETRF; for 1<=i<=N, row i of the
                matrix was interchanged with row IPIV(i).

        DWORK    (workspace/output) COMPLEX*16 array on the GPU, dimension (MAX(1,LWORK))
      
        LWORK   (input) INTEGER
                The dimension of the array DWORK.  LWORK >= N*NB, where NB is
                the optimal blocksize returned by magma_get_zgetri_nb(n).
                
                Unlike LAPACK, this version does not currently support a
                workspace query, because the workspace is on the GPU.

        INFO    (output) INTEGER
                = 0:  successful exit
                < 0:  if INFO = -i, the i-th argument had an illegal value
                > 0:  if INFO = i, U(i,i) is exactly zero; the matrix is
                      singular and its cannot be computed.

  ===================================================================== */

    /* Local variables */
    magma_int_t ret;
    cuDoubleComplex c_one = MAGMA_Z_ONE;
    cuDoubleComplex c_neg_one = MAGMA_Z_NEG_ONE;
    cuDoubleComplex *dL = dwork;
    magma_int_t     ldl = n;
    magma_int_t      nb = magma_get_zgetri_nb(n);
    magma_int_t j, jmax, jb, jp;
    
    *info = 0;
    if (n < 0)
        *info = -1;
    else if (lda < max(1,n))
        *info = -3;
    else if ( lwork < n*nb )
        *info = -6;

    if (*info != 0) {
        magma_xerbla( __func__, -(*info) );
        return MAGMA_ERR_ILLEGAL_VALUE;
    }

    /* Quick return if possible */
    if ( n == 0 )
        return MAGMA_SUCCESS;
    
    /* Invert the triangular factor U */
    ret = magma_ztrtri_gpu( MagmaUpper, MagmaNonUnit, n, dA, lda, info );
    if ( *info != 0 )
        return ret;
    
    jmax = ((n-1) / nb)*nb;
    for( j = jmax; j >= 0; j -= nb ) {
        jb = min( nb, n-j );
        
        // copy current block column of L to work space,
        // then replace with zeros in A.
        magmablas_zlacpy( MagmaUpperLower, n-j, jb,
                          &dA[j + j*lda], lda,
                          &dL[j        ], ldl );
        magmablas_zlaset( MagmaLower, n-j, jb, &dA[j + j*lda], lda );
        
        // compute current block column of Ainv
        // Ainv(:, j:j+jb-1)
        //   = ( U(:, j:j+jb-1) - Ainv(:, j+jb:n) L(j+jb:n, j:j+jb-1) )
        //   * L(j:j+jb-1, j:j+jb-1)^{-1}
        // where L(:, j:j+jb-1) is stored in dL.
        if ( j+jb < n ) {
            cublasZgemm( MagmaNoTrans, MagmaNoTrans, n, jb, n-j-jb,
                         c_neg_one, &dA[(j+jb)*lda], lda,
                                    &dL[ j+jb     ], ldl,
                         c_one,     &dA[     j*lda], lda );
        }
        cublasZtrsm( MagmaRight, MagmaLower, MagmaNoTrans, MagmaUnit,
                     n, jb, c_one,
                     &dL[j    ], ldl,
                     &dA[j*lda], lda );
    }

    // Apply column interchanges
    for( j = n-2; j >= 0; --j ) {
        jp = ipiv[j] - 1;
        if ( jp != j ) {
            magmablas_zswap( n, &dA[ j*lda ], 1, &dA[ jp*lda ], 1 );
        }
    }
    
    return MAGMA_SUCCESS;
}
