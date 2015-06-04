/*
    -- MAGMA (version 1.5.0-beta3) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date July 2014

       @author Stan Tomov
       @generated from zpotrf.cpp normal z -> d, Fri Jul 18 17:34:15 2014
*/
#include "common_magma.h"

#define PRECISION_d

// === Define what BLAS to use ============================================
//#if defined(PRECISION_s) || defined(PRECISION_d)
    #define magma_dtrsm magmablas_dtrsm
//#endif
// === End defining what BLAS to use ======================================

/**
    Purpose
    -------
    DPOTRF computes the Cholesky factorization of a real symmetric
    positive definite matrix A. This version does not require work
    space on the GPU passed as input. GPU memory is allocated in the
    routine.

    The factorization has the form
        A = U**T * U,  if uplo = MagmaUpper, or
        A = L  * L**T, if uplo = MagmaLower,
    where U is an upper triangular matrix and L is lower triangular.

    This is the block version of the algorithm, calling Level 3 BLAS.
    If the current stream is NULL, this version replaces it with user defined
    stream to overlap computation with communication.

    Arguments
    ---------
    @param[in]
    uplo    magma_uplo_t
      -     = MagmaUpper:  Upper triangle of A is stored;
      -     = MagmaLower:  Lower triangle of A is stored.

    @param[in]
    n       INTEGER
            The order of the matrix A.  N >= 0.

    @param[in,out]
    A       DOUBLE_PRECISION array, dimension (LDA,N)
            On entry, the symmetric matrix A.  If uplo = MagmaUpper, the leading
            N-by-N upper triangular part of A contains the upper
            triangular part of the matrix A, and the strictly lower
            triangular part of A is not referenced.  If uplo = MagmaLower, the
            leading N-by-N lower triangular part of A contains the lower
            triangular part of the matrix A, and the strictly upper
            triangular part of A is not referenced.
    \n
            On exit, if INFO = 0, the factor U or L from the Cholesky
            factorization A = U**T * U or A = L * L**T.
    \n
            Higher performance is achieved if A is in pinned memory, e.g.
            allocated using magma_malloc_pinned.

    @param[in]
    lda     INTEGER
            The leading dimension of the array A.  LDA >= max(1,N).

    @param[out]
    info    INTEGER
      -     = 0:  successful exit
      -     < 0:  if INFO = -i, the i-th argument had an illegal value
                  or another error occured, such as memory allocation failed.
      -     > 0:  if INFO = i, the leading minor of order i is not
                  positive definite, and the factorization could not be
                  completed.

    @ingroup magma_dposv_comp
    ********************************************************************/
extern "C" magma_int_t
magma_dpotrf(magma_uplo_t uplo, magma_int_t n,
             double *A, magma_int_t lda, magma_int_t *info)
{
#define A(i, j)  (A    + (j)*lda  + (i))
#define dA(i, j) (work + (j)*ldda + (i))

    /* Local variables */
    const char* uplo_ = lapack_uplo_const( uplo );
    magma_int_t        ldda, nb;
    magma_int_t j, jb;
    double    c_one     = MAGMA_D_ONE;
    double    c_neg_one = MAGMA_D_NEG_ONE;
    double   *work;
    double             d_one     =  1.0;
    double             d_neg_one = -1.0;
    int upper = (uplo == MagmaUpper);

    *info = 0;
    if (! upper && uplo != MagmaLower) {
        *info = -1;
    } else if (n < 0) {
        *info = -2;
    } else if (lda < max(1,n)) {
        *info = -4;
    }
    if (*info != 0) {
        magma_xerbla( __func__, -(*info) );
        return *info;
    }

    /* Quick return */
    if ( n == 0 )
        return *info;

    magma_int_t num_gpus = magma_num_gpus();
    if ( num_gpus > 1 ) {
        /* call multiple-GPU interface  */
        return magma_dpotrf_m(num_gpus, uplo, n, A, lda, info);
    }

    ldda = ((n+31)/32)*32;
    
    if (MAGMA_SUCCESS != magma_dmalloc( &work, (n)*ldda )) {
        /* alloc failed so call the non-GPU-resident version */
        return magma_dpotrf_m(num_gpus, uplo, n, A, lda, info);
    }

    /* Define user stream if current stream is NULL */
    cudaStream_t stream[3], current_stream;
    magmablasGetKernelStream(&current_stream);

    magma_queue_create( &stream[0] );
    magma_queue_create( &stream[2] );

    if (current_stream == NULL) {
        magma_queue_create( &stream[1] );
        magmablasSetKernelStream(stream[1]);
    }
    else
        stream[1] = current_stream;

    nb = magma_get_dpotrf_nb(n);

    if (nb <= 1 || nb >= n) {
        lapackf77_dpotrf(uplo_, &n, A, &lda, info);
    } else {
        /* Use hybrid blocked code. */
        if (upper) {
            /* Compute the Cholesky factorization A = U'*U. */
            for (j=0; j < n; j += nb) {
                /* Update and factorize the current diagonal block and test
                   for non-positive-definiteness. Computing MIN */
                jb = min(nb, (n-j));
                magma_dsetmatrix_async( jb, (n-j), A(j, j), lda, dA(j, j), ldda, stream[1]);
                
                magma_dsyrk(MagmaUpper, MagmaTrans, jb, j,
                            d_neg_one, dA(0, j), ldda,
                            d_one,     dA(j, j), ldda);
                magma_queue_sync( stream[1] );

                magma_dgetmatrix_async( jb, jb,
                                        dA(j, j), ldda,
                                        A(j, j),  lda, stream[0] );
                
                if ( (j+jb) < n) {
                    magma_dgemm(MagmaTrans, MagmaNoTrans,
                                jb, (n-j-jb), j,
                                c_neg_one, dA(0, j   ), ldda,
                                           dA(0, j+jb), ldda,
                                c_one,     dA(j, j+jb), ldda);
                }
                
                magma_dgetmatrix_async( j, jb,
                                        dA(0, j), ldda,
                                        A (0, j),  lda, stream[2] );

                magma_queue_sync( stream[0] );
                lapackf77_dpotrf(MagmaUpperStr, &jb, A(j, j), &lda, info);
                if (*info != 0) {
                    *info = *info + j;
                    break;
                }
                magma_dsetmatrix_async( jb, jb,
                                        A(j, j),  lda,
                                        dA(j, j), ldda, stream[0] );
                magma_queue_sync( stream[0] );

                if ( (j+jb) < n ) {
                    magma_dtrsm(MagmaLeft, MagmaUpper, MagmaTrans, MagmaNonUnit,
                                jb, (n-j-jb),
                                c_one, dA(j, j   ), ldda,
                                       dA(j, j+jb), ldda);
                }
            }
        }
        else {
            //=========================================================
            // Compute the Cholesky factorization A = L*L'.
            for (j=0; j < n; j += nb) {
                //  Update and factorize the current diagonal block and test
                //  for non-positive-definiteness. Computing MIN
                jb = min(nb, (n-j));
                magma_dsetmatrix_async( (n-j), jb, A(j, j), lda, dA(j, j), ldda, stream[1]);

                magma_dsyrk(MagmaLower, MagmaNoTrans, jb, j,
                            d_neg_one, dA(j, 0), ldda,
                            d_one,     dA(j, j), ldda);
                magma_queue_sync( stream[1] );

                magma_dgetmatrix_async( jb, jb,
                                        dA(j,j), ldda,
                                        A(j,j),  lda, stream[0] );

                if ( (j+jb) < n) {
                    magma_dgemm( MagmaNoTrans, MagmaTrans,
                                 (n-j-jb), jb, j,
                                 c_neg_one, dA(j+jb, 0), ldda,
                                            dA(j,    0), ldda,
                                 c_one,     dA(j+jb, j), ldda);
                }
                
                magma_dgetmatrix_async( jb, j,
                                        dA(j, 0), ldda,
                                        A(j, 0),  lda, stream[2] );

                magma_queue_sync( stream[0] );
                lapackf77_dpotrf(MagmaLowerStr, &jb, A(j, j), &lda, info);
                if (*info != 0) {
                    *info = *info + j;
                    break;
                }
                magma_dsetmatrix_async( jb, jb,
                                        A(j, j),  lda,
                                        dA(j, j), ldda, stream[0] );
                magma_queue_sync( stream[0] );

                if ( (j+jb) < n) {
                    magma_dtrsm(MagmaRight, MagmaLower, MagmaTrans, MagmaNonUnit,
                                (n-j-jb), jb,
                                c_one, dA(j,    j), ldda,
                                       dA(j+jb, j), ldda);
                }
            }
        }
    }
    
    magma_queue_destroy( stream[0] );
    magma_queue_destroy( stream[2] );
    if (current_stream == NULL) {
        magma_queue_destroy( stream[1] );
        magmablasSetKernelStream(NULL);
    }

    magma_free( work );
    
    return *info;
} /* magma_dpotrf */
