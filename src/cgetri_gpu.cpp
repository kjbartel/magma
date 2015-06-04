/*
    -- MAGMA (version 1.2.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       May 2012

       @generated c Tue May 15 18:17:31 2012

*/
#include "common_magma.h"

// === Define what BLAS to use ============================================
#define PRECISION_c
#if (defined(PRECISION_s) || defined(PRECISION_d))
  #define magma_cgemm magmablas_cgemm
  #define magma_ctrsm magmablas_ctrsm
#endif

#if (GPUSHMEM >= 200)
#if (defined(PRECISION_s))
    #undef  magma_sgemm
    #define magma_sgemm magmablas_sgemm_fermi80
#endif
#endif
// === End defining what BLAS to use ======================================

extern "C" magma_int_t
magma_cgetri_gpu( magma_int_t n, cuFloatComplex *dA, magma_int_t lda,
                  magma_int_t *ipiv, cuFloatComplex *dwork, magma_int_t lwork,
                  magma_int_t *info )
{
/*  -- MAGMA (version 1.2.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       May 2012

    Purpose
    =======

        CGETRI computes the inverse of a matrix using the LU factorization
        computed by CGETRF. This method inverts U and then computes inv(A) by
        solving the system inv(A)*L = inv(U) for inv(A).
        
        Note that it is generally both faster and more accurate to use CGESV,
        or CGETRF and CGETRS, to solve the system AX = B, rather than inverting
        the matrix and multiplying to form X = inv(A)*B. Only in special
        instances should an explicit inverse be computed with this routine.

    Arguments
    =========

        N       (input) INTEGER
                The order of the matrix A.  N >= 0.

        dA      (input/output) COMPLEX array on the GPU, dimension (LDA,N)
                On entry, the factors L and U from the factorization
                A = P*L*U as computed by CGETRF_GPU.
                On exit, if INFO = 0, the inverse of the original matrix A.

        LDA     (input) INTEGER
                The leading dimension of the array A.  LDA >= max(1,N).

        IPIV    (input) INTEGER array, dimension (N)
                The pivot indices from CGETRF; for 1<=i<=N, row i of the
                matrix was interchanged with row IPIV(i).

        DWORK    (workspace/output) COMPLEX*16 array on the GPU, dimension (MAX(1,LWORK))
      
        LWORK   (input) INTEGER
                The dimension of the array DWORK.  LWORK >= N*NB, where NB is
                the optimal blocksize returned by magma_get_cgetri_nb(n).
                
                Unlike LAPACK, this version does not currently support a
                workspace query, because the workspace is on the GPU.

        INFO    (output) INTEGER
                = 0:  successful exit
                < 0:  if INFO = -i, the i-th argument had an illegal value
                > 0:  if INFO = i, U(i,i) is exactly zero; the matrix is
                      singular and its cannot be computed.

  ===================================================================== */

    /* Local variables */
    cuFloatComplex c_one     = MAGMA_C_ONE;
    cuFloatComplex c_neg_one = MAGMA_C_NEG_ONE;
    cuFloatComplex *dL = dwork;
    magma_int_t     ldl = n;
    magma_int_t      nb = magma_get_cgetri_nb(n);
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
        return *info;
    }

    /* Quick return if possible */
    if ( n == 0 )
        return *info;
    
    /* Invert the triangular factor U */
    magma_ctrtri_gpu( MagmaUpper, MagmaNonUnit, n, dA, lda, info );
    if ( *info != 0 )
        return *info;
    
    jmax = ((n-1) / nb)*nb;
    for( j = jmax; j >= 0; j -= nb ) {
        jb = min( nb, n-j );
        
        // copy current block column of L to work space,
        // then replace with zeros in A.
        magmablas_clacpy( MagmaUpperLower, n-j, jb,
                          &dA[j + j*lda], lda,
                          &dL[j        ], ldl );
        magmablas_claset( MagmaLower, n-j, jb, &dA[j + j*lda], lda );
        
        // compute current block column of Ainv
        // Ainv(:, j:j+jb-1)
        //   = ( U(:, j:j+jb-1) - Ainv(:, j+jb:n) L(j+jb:n, j:j+jb-1) )
        //   * L(j:j+jb-1, j:j+jb-1)^{-1}
        // where L(:, j:j+jb-1) is stored in dL.
        if ( j+jb < n ) {
            magma_cgemm( MagmaNoTrans, MagmaNoTrans, n, jb, n-j-jb,
                         c_neg_one, &dA[(j+jb)*lda], lda,
                                    &dL[ j+jb     ], ldl,
                         c_one,     &dA[     j*lda], lda );
        }
        magma_ctrsm( MagmaRight, MagmaLower, MagmaNoTrans, MagmaUnit,
                     n, jb, c_one,
                     &dL[j    ], ldl,
                     &dA[j*lda], lda );
    }

    // Apply column interchanges
    for( j = n-2; j >= 0; --j ) {
        jp = ipiv[j] - 1;
        if ( jp != j ) {
            magmablas_cswap( n, &dA[ j*lda ], 1, &dA[ jp*lda ], 1 );
        }
    }
    
    return *info;
}
