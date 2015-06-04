/*
    -- MAGMA (version 1.3.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2012

       @generated d Wed Nov 14 22:53:02 2012

*/
#include "common_magma.h"

// === Define what BLAS to use ============================================
#define PRECISION_d

#if (defined(PRECISION_s) || defined(PRECISION_d))
  #define magma_dgemm magmablas_dgemm
  #define magma_dtrsm magmablas_dtrsm
#endif

#if (GPUSHMEM >= 200) && defined(PRECISION_s)
  #undef  magma_sgemm
  #define magma_sgemm magmablas_sgemm_fermi80
#endif
// === End defining what BLAS to use ======================================

#define A(i, j)  (a   +(j)*lda  + (i))
#define dA(i, j) (work+(j)*ldda + (i))

extern "C" magma_int_t
magma_dtrtri(char uplo, char diag, magma_int_t n,
              double *a, magma_int_t lda, magma_int_t *info)
{
/*  -- MAGMA (version 1.3.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2012

    Purpose
    =======

    DTRTRI computes the inverse of a real upper or lower triangular
    matrix A.

    This is the Level 3 BLAS version of the algorithm.

    Arguments
    =========

    UPLO    (input) CHARACTER*1
            = 'U':  A is upper triangular;
            = 'L':  A is lower triangular.

    DIAG    (input) CHARACTER*1
            = 'N':  A is non-unit triangular;
            = 'U':  A is unit triangular.

    N       (input) INTEGER
            The order of the matrix A.  N >= 0.

    A       (input/output) DOUBLE_PRECISION array, dimension (LDA,N)
            On entry, the triangular matrix A.  If UPLO = 'U', the
            leading N-by-N upper triangular part of the array A contains
            the upper triangular matrix, and the strictly lower
            triangular part of A is not referenced.  If UPLO = 'L', the
            leading N-by-N lower triangular part of the array A contains
            the lower triangular matrix, and the strictly upper
            triangular part of A is not referenced.  If DIAG = 'U', the
            diagonal elements of A are also not referenced and are
            assumed to be 1.
            On exit, the (triangular) inverse of the original matrix, in
            the same storage format.

    LDA     (input) INTEGER
            The leading dimension of the array A.  LDA >= max(1,N).

    INFO    (output) INTEGER
            = 0: successful exit
            < 0: if INFO = -i, the i-th argument had an illegal value
            > 0: if INFO = i, A(i,i) is exactly zero.  The triangular
                    matrix is singular and its inverse cannot be computed.

    ===================================================================== */

    /* Local variables */
    char uplo_[2] = {uplo, 0};
    char diag_[2] = {diag, 0};
    magma_int_t     ldda, nb, nn, j, jb;
    double c_zero     = MAGMA_D_ZERO;
    double c_one      = MAGMA_D_ONE;
    double c_neg_one  = MAGMA_D_NEG_ONE;
    double *work;

    int upper  = lapackf77_lsame(uplo_, "U");
    int nounit = lapackf77_lsame(diag_, "N");

    *info = 0;

    if ((! upper) && (! lapackf77_lsame(uplo_, "L")))
        *info = -1;
    else if ((! nounit) && (! lapackf77_lsame(diag_, "U")))
        *info = -2;
    else if (n < 0)
        *info = -3;
    else if (lda < max(1,n))
        *info = -5;

    if (*info != 0) {
        magma_xerbla( __func__, -(*info) );
        return *info;
    }

    /* Quick return */
    if ( n == 0 )
        return *info;

    /* Check for singularity if non-unit */
    if (nounit) {
        for ( j=0; j<n; ++j ) {
            if ( MAGMA_D_EQUAL( *A(j,j), c_zero )) {
                *info = j+1;  // Fortran index
                return *info;
            }
        }
    }

    /* Determine the block size for this environment */
    nb = magma_get_dpotrf_nb(n);

    ldda = ((n+31)/32)*32;
    if (MAGMA_SUCCESS != magma_dmalloc( &work, (n)*ldda )) {
        *info = MAGMA_ERR_DEVICE_ALLOC;
        return *info;
    }

    cudaStream_t stream[2];
    magma_queue_create( &stream[0] );
    magma_queue_create( &stream[1] );

    if (nb <= 1 || nb >= n)
        lapackf77_dtrtri(uplo_, diag_, &n, a, &lda, info);
    else {
        if (upper) {
            /* Compute inverse of upper triangular matrix */
            for (j=0; j<n; j=j+nb) {
                jb = min(nb, (n-j));
                magma_dsetmatrix( jb, (n-j),
                                  A(j, j),  lda,
                                  dA(j, j), ldda );

                /* Compute rows 1:j-1 of current block column */
                magma_dtrmm( MagmaLeft, MagmaUpper,
                             MagmaNoTrans, MagmaNonUnit, j, jb,
                             c_one, dA(0,0), ldda, dA(0, j),ldda);

                magma_dtrsm( MagmaRight, MagmaUpper,
                             MagmaNoTrans, MagmaNonUnit, j, jb,
                             c_neg_one, dA(j,j), ldda, dA(0, j),ldda);

                //cublasGetMatrix(j ,jb, sizeof( double),
                //dA(0, j), ldda, A(0, j), lda);

                magma_dgetmatrix_async( jb, jb,
                                        dA(j, j), ldda,
                                        A(j, j),  lda, stream[1] );

                magma_dgetmatrix_async( j, jb,
                                        dA(0, j), ldda,
                                        A(0, j),  lda, stream[0] );

                magma_queue_sync( stream[1] );

                /* Compute inverse of current diagonal block */
                lapackf77_dtrtri(MagmaUpperStr, diag_, &jb, A(j,j), &lda, info);

                magma_dsetmatrix( jb, jb,
                                  A(j, j),  lda,
                                  dA(j, j), ldda );
            }
        }
        else {
            /* Compute inverse of lower triangular matrix */
            nn=((n-1)/nb)*nb+1;

            for(j=nn-1; j>=0; j=j-nb) {
                jb=min(nb,(n-j));

                if((j+jb) < n) {
                    magma_dsetmatrix( (n-j), jb,
                                      A(j, j),  lda,
                                      dA(j, j), ldda );

                    /* Compute rows j+jb:n of current block column */
                    magma_dtrmm( MagmaLeft, MagmaLower,
                                 MagmaNoTrans, MagmaNonUnit, (n-j-jb), jb,
                                 c_one, dA(j+jb,j+jb), ldda, dA(j+jb, j), ldda );

                    magma_dtrsm( MagmaRight, MagmaLower,
                                 MagmaNoTrans, MagmaNonUnit, (n-j-jb), jb,
                                 c_neg_one, dA(j,j), ldda, dA(j+jb, j), ldda );

                    //cublasGetMatrix((n-j), jb, sizeof( double),dA(j, j), ldda, A(j, j), lda);

                    magma_dgetmatrix_async( n-j-jb, jb,
                                            dA(j+jb, j), ldda,
                                            A(j+jb, j),  lda, stream[1] );

                    magma_dgetmatrix_async( jb, jb,
                                            dA(j,j), ldda,
                                            A(j,j),  lda, stream[0] );

                    magma_queue_sync( stream[0] );
                }

                /* Compute inverse of current diagonal block */
                lapackf77_dtrtri(MagmaLowerStr, diag_, &jb, A(j,j), &lda, info);

                magma_dsetmatrix( jb, jb,
                                  A(j, j),  lda,
                                  dA(j, j), ldda );
            }
        }
    }

    magma_queue_destroy( stream[0] );
    magma_queue_destroy( stream[1] );
    magma_free( work );

    return *info;
}
