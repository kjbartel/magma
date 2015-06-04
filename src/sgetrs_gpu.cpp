/*
    -- MAGMA (version 1.4.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       December 2013

       @generated s Tue Dec 17 13:18:36 2013

*/
#include "common_magma.h"

extern "C" magma_int_t
magma_sgetrs_gpu(char trans, magma_int_t n, magma_int_t nrhs,
                 float *dA, magma_int_t ldda,
                 magma_int_t *ipiv,
                 float *dB, magma_int_t lddb,
                 magma_int_t *info)
{
/*  -- MAGMA (version 1.4.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       December 2013

    Purpose
    =======
    Solves a system of linear equations
      A * X = B  or  A' * X = B
    with a general N-by-N matrix A using the LU factorization computed by SGETRF_GPU.

    Arguments
    =========
    TRANS   (input) CHARACTER*1
            Specifies the form of the system of equations:
            = 'N':  A * X = B  (No transpose)
            = 'T':  A'* X = B  (Transpose)
            = 'C':  A'* X = B  (Conjugate transpose = Transpose)

    N       (input) INTEGER
            The order of the matrix A.  N >= 0.

    NRHS    (input) INTEGER
            The number of right hand sides, i.e., the number of columns
            of the matrix B.  NRHS >= 0.

    A       (input) REAL array on the GPU, dimension (LDA,N)
            The factors L and U from the factorization A = P*L*U as computed
            by SGETRF_GPU.

    LDA     (input) INTEGER
            The leading dimension of the array A.  LDA >= max(1,N).

    IPIV    (input) INTEGER array, dimension (N)
            The pivot indices from SGETRF; for 1<=i<=N, row i of the
            matrix was interchanged with row IPIV(i).

    B       (input/output) REAL array on the GPU, dimension (LDB,NRHS)
            On entry, the right hand side matrix B.
            On exit, the solution matrix X.

    LDB     (input) INTEGER
            The leading dimension of the array B.  LDB >= max(1,N).

    INFO    (output) INTEGER
            = 0:  successful exit
            < 0:  if INFO = -i, the i-th argument had an illegal value

    HWORK   (workspace) REAL array, dimension N*NRHS
    =====================================================================    */


    float c_one = MAGMA_S_ONE;
    float *work = NULL;
    char            trans_[2] = {trans, 0};
    int notran = lapackf77_lsame(trans_, "N");
    magma_int_t i1, i2, inc;

    *info = 0;
    if ( (! notran) &&
         (! lapackf77_lsame(trans_, "T")) &&
         (! lapackf77_lsame(trans_, "C")) ) {
        *info = -1;
    } else if (n < 0) {
        *info = -2;
    } else if (nrhs < 0) {
        *info = -3;
    } else if (ldda < max(1,n)) {
        *info = -5;
    } else if (lddb < max(1,n)) {
        *info = -8;
    }
    if (*info != 0) {
        magma_xerbla( __func__, -(*info) );
        return *info;
    }

    /* Quick return if possible */
    if (n == 0 || nrhs == 0) {
        return *info;
    }

    magma_smalloc_cpu( &work, n * nrhs );
    if ( work == NULL ) {
        *info = MAGMA_ERR_HOST_ALLOC;
        return *info;
    }
      
    i1 = 1;
    i2 = n;
    if (notran) {
        inc = 1;

        /* Solve A * X = B. */
        magma_sgetmatrix( n, nrhs, dB, lddb, work, n );
        lapackf77_slaswp(&nrhs, work, &n, &i1, &i2, ipiv, &inc);
        magma_ssetmatrix( n, nrhs, work, n, dB, lddb );

        if ( nrhs == 1) {
            magma_strsv(MagmaLower, MagmaNoTrans, MagmaUnit,    n, dA, ldda, dB, 1 );
            magma_strsv(MagmaUpper, MagmaNoTrans, MagmaNonUnit, n, dA, ldda, dB, 1 );
        } else {
            magma_strsm(MagmaLeft, MagmaLower, MagmaNoTrans, MagmaUnit,    n, nrhs, c_one, dA, ldda, dB, lddb );
            magma_strsm(MagmaLeft, MagmaUpper, MagmaNoTrans, MagmaNonUnit, n, nrhs, c_one, dA, ldda, dB, lddb );
        }
    } else {
        inc = -1;

        /* Solve A' * X = B. */
        if ( nrhs == 1) {
            magma_strsv(MagmaUpper, trans, MagmaNonUnit, n, dA, ldda, dB, 1 );
            magma_strsv(MagmaLower, trans, MagmaUnit,    n, dA, ldda, dB, 1 );
        } else {
            magma_strsm(MagmaLeft, MagmaUpper, trans, MagmaNonUnit, n, nrhs, c_one, dA, ldda, dB, lddb );
            magma_strsm(MagmaLeft, MagmaLower, trans, MagmaUnit,    n, nrhs, c_one, dA, ldda, dB, lddb );
        }

        magma_sgetmatrix( n, nrhs, dB, lddb, work, n );
        lapackf77_slaswp(&nrhs, work, &n, &i1, &i2, ipiv, &inc);
        magma_ssetmatrix( n, nrhs, work, n, dB, lddb );
    }
    magma_free_cpu(work);

    return *info;
}

