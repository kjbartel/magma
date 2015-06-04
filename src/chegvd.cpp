/*
    -- MAGMA (version 1.2.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       June 2012

       @author Raffaele Solca
       @author Stan Tomov

       @generated c Thu Jun 28 12:31:01 2012

*/
#include "common_magma.h"

/* This ctrmm interface is used for TAU profiling */
void Mymagma_ctrmm(char side, char uplo, char trans, char unit,
                   magma_int_t n, magma_int_t m,
                   cuFloatComplex alpha, cuFloatComplex *db, magma_int_t lddb,
                   cuFloatComplex *dz, magma_int_t lddz)
{
  magma_ctrmm(side, uplo, trans, unit, n, m, alpha, db, lddb, dz, lddz);
  magma_device_sync();
}

/* This ctrsm interface is used for TAU profiling */
void Mymagma_ctrsm(char side, char uplo, char trans, char unit,
                   magma_int_t n, magma_int_t m,
                   cuFloatComplex alpha, cuFloatComplex *db, magma_int_t lddb,
                   cuFloatComplex *dz, magma_int_t lddz)
{
  magma_ctrsm(side, uplo, trans, unit, n, m, alpha, db, lddb, dz, lddz);
  magma_device_sync();
}

extern "C" magma_int_t
magma_chegvd(magma_int_t itype, char jobz, char uplo, magma_int_t n,
             cuFloatComplex *a, magma_int_t lda, cuFloatComplex *b, magma_int_t ldb,
             float *w, cuFloatComplex *work, magma_int_t lwork,
             float *rwork, magma_int_t lrwork,
             magma_int_t *iwork, magma_int_t liwork, magma_int_t *info)
{
/*  -- MAGMA (version 1.2.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       June 2012

    Purpose
    =======
    CHEGVD computes all the eigenvalues, and optionally, the eigenvectors
    of a complex generalized Hermitian-definite eigenproblem, of the form
    A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.  Here A and
    B are assumed to be Hermitian and B is also positive definite.
    If eigenvectors are desired, it uses a divide and conquer algorithm.

    The divide and conquer algorithm makes very mild assumptions about
    floating point arithmetic. It will work on machines with a guard
    digit in add/subtract, or on those binary machines without guard
    digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or
    Cray-2. It could conceivably fail on hexadecimal or decimal machines
    without guard digits, but we know of none.

    Arguments
    =========
    ITYPE   (input) INTEGER
            Specifies the problem type to be solved:
            = 1:  A*x = (lambda)*B*x
            = 2:  A*B*x = (lambda)*x
            = 3:  B*A*x = (lambda)*x

    JOBZ    (input) CHARACTER*1
            = 'N':  Compute eigenvalues only;
            = 'V':  Compute eigenvalues and eigenvectors.

    UPLO    (input) CHARACTER*1
            = 'U':  Upper triangles of A and B are stored;
            = 'L':  Lower triangles of A and B are stored.

    N       (input) INTEGER
            The order of the matrices A and B.  N >= 0.

    A       (input/output) COMPLEX*16 array, dimension (LDA, N)
            On entry, the Hermitian matrix A.  If UPLO = 'U', the
            leading N-by-N upper triangular part of A contains the
            upper triangular part of the matrix A.  If UPLO = 'L',
            the leading N-by-N lower triangular part of A contains
            the lower triangular part of the matrix A.

            On exit, if JOBZ = 'V', then if INFO = 0, A contains the
            matrix Z of eigenvectors.  The eigenvectors are normalized
            as follows:
            if ITYPE = 1 or 2, Z**H*B*Z = I;
            if ITYPE = 3, Z**H*inv(B)*Z = I.
            If JOBZ = 'N', then on exit the upper triangle (if UPLO='U')
            or the lower triangle (if UPLO='L') of A, including the
            diagonal, is destroyed.

    LDA     (input) INTEGER
            The leading dimension of the array A.  LDA >= max(1,N).

    B       (input/output) COMPLEX*16 array, dimension (LDB, N)
            On entry, the Hermitian matrix B.  If UPLO = 'U', the
            leading N-by-N upper triangular part of B contains the
            upper triangular part of the matrix B.  If UPLO = 'L',
            the leading N-by-N lower triangular part of B contains
            the lower triangular part of the matrix B.

            On exit, if INFO <= N, the part of B containing the matrix is
            overwritten by the triangular factor U or L from the Cholesky
            factorization B = U**H*U or B = L*L**H.

    LDB     (input) INTEGER
            The leading dimension of the array B.  LDB >= max(1,N).

    W       (output) DOUBLE PRECISION array, dimension (N)
            If INFO = 0, the eigenvalues in ascending order.

    WORK    (workspace/output) COMPLEX*16 array, dimension (MAX(1,LWORK))
            On exit, if INFO = 0, WORK(1) returns the optimal LWORK.

    LWORK   (input) INTEGER
            The length of the array WORK.
            If N <= 1,                LWORK >= 1.
            If JOBZ  = 'N' and N > 1, LWORK >= N + 1.
            If JOBZ  = 'V' and N > 1, LWORK >= 2*N*nb + N**2.

            If LWORK = -1, then a workspace query is assumed; the routine
            only calculates the optimal sizes of the WORK, RWORK and
            IWORK arrays, returns these values as the first entries of
            the WORK, RWORK and IWORK arrays, and no error message
            related to LWORK or LRWORK or LIWORK is issued by XERBLA.

    RWORK   (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LRWORK))
            On exit, if INFO = 0, RWORK(1) returns the optimal LRWORK.

    LRWORK  (input) INTEGER
            The dimension of the array RWORK.
            If N <= 1,                LRWORK >= 1.
            If JOBZ  = 'N' and N > 1, LRWORK >= N.
            If JOBZ  = 'V' and N > 1, LRWORK >= 1 + 5*N + 2*N**2.

            If LRWORK = -1, then a workspace query is assumed; the
            routine only calculates the optimal sizes of the WORK, RWORK
            and IWORK arrays, returns these values as the first entries
            of the WORK, RWORK and IWORK arrays, and no error message
            related to LWORK or LRWORK or LIWORK is issued by XERBLA.

    IWORK   (workspace/output) INTEGER array, dimension (MAX(1,LIWORK))
            On exit, if INFO = 0, IWORK(1) returns the optimal LIWORK.

    LIWORK  (input) INTEGER
            The dimension of the array IWORK.
            If N <= 1,                LIWORK >= 1.
            If JOBZ  = 'N' and N > 1, LIWORK >= 1.
            If JOBZ  = 'V' and N > 1, LIWORK >= 3 + 5*N.

            If LIWORK = -1, then a workspace query is assumed; the
            routine only calculates the optimal sizes of the WORK, RWORK
            and IWORK arrays, returns these values as the first entries
            of the WORK, RWORK and IWORK arrays, and no error message
            related to LWORK or LRWORK or LIWORK is issued by XERBLA.

    INFO    (output) INTEGER
            = 0:  successful exit
            < 0:  if INFO = -i, the i-th argument had an illegal value
            > 0:  CPOTRF or CHEEVD returned an error code:
               <= N:  if INFO = i and JOBZ = 'N', then the algorithm
                      failed to converge; i off-diagonal elements of an
                      intermediate tridiagonal form did not converge to
                      zero;
                      if INFO = i and JOBZ = 'V', then the algorithm
                      failed to compute an eigenvalue while working on
                      the submatrix lying in rows and columns INFO/(N+1)
                      through mod(INFO,N+1);
               > N:   if INFO = N + i, for 1 <= i <= N, then the leading
                      minor of order i of B is not positive definite.
                      The factorization of B could not be completed and
                      no eigenvalues or eigenvectors were computed.

    Further Details
    ===============

    Based on contributions by
       Mark Fahey, Department of Mathematics, Univ. of Kentucky, USA

    Modified so that no backsubstitution is performed if CHEEVD fails to
    converge (NEIG in old code could be greater than N causing out of
    bounds reference to A - reported by Ralf Meyer).  Also corrected the
    description of INFO and the test on ITYPE. Sven, 16 Feb 05.
    =====================================================================  */

    char uplo_[2] = {uplo, 0};
    char jobz_[2] = {jobz, 0};

    cuFloatComplex c_one = MAGMA_C_ONE;

    cuFloatComplex *da;
    cuFloatComplex *db;
    magma_int_t ldda = n;
    magma_int_t lddb = n;

    magma_int_t lower;
    char trans[1];
    magma_int_t wantz;
    magma_int_t lquery;

    //    magma_int_t lopt;
    magma_int_t lwmin;
    //    magma_int_t liopt;
    magma_int_t liwmin;
    //    magma_int_t lropt;
    magma_int_t lrwmin;

    cudaStream_t stream;
    magma_queue_create( &stream );

    wantz = lapackf77_lsame(jobz_, MagmaVectorsStr);
    lower = lapackf77_lsame(uplo_, MagmaLowerStr);
    lquery = lwork == -1 || lrwork == -1 || liwork == -1;

    *info = 0;
    if (itype < 1 || itype > 3) {
        *info = -1;
    } else if (! (wantz || lapackf77_lsame(jobz_, MagmaNoVectorsStr))) {
        *info = -2;
    } else if (! (lower || lapackf77_lsame(uplo_, MagmaUpperStr))) {
        *info = -3;
    } else if (n < 0) {
        *info = -4;
    } else if (lda < max(1,n)) {
        *info = -6;
    } else if (ldb < max(1,n)) {
        *info = -8;
    }

    magma_int_t nb = magma_get_chetrd_nb(n);

    if (wantz) {
        lwmin = 2 * n + n * n;
        lrwmin = 1 + 5 * n + 2 * n * n;
        liwmin = 5 * n + 3;
    } else {
        lwmin = n * (nb + 1);
        lrwmin = n;
        liwmin = 1;
    }

    MAGMA_C_SET2REAL(work[0],(float)lwmin);
    rwork[0] = lrwmin;
    iwork[0] = liwmin;

    if (lwork < lwmin && ! lquery) {
        *info = -11;
    } else if (lrwork < lrwmin && ! lquery) {
        *info = -13;
    } else if (liwork < liwmin && ! lquery) {
         *info = -15;
    }

    if (*info != 0) {
        magma_xerbla( __func__, -(*info) );
        return *info;
    }
    else if (lquery) {
        return *info;
    }

    /*     Quick return if possible */
    if (n == 0) {
        return *info;
    }

    if (MAGMA_SUCCESS != magma_cmalloc( &da, n*ldda ) ||
        MAGMA_SUCCESS != magma_cmalloc( &db, n*lddb )) {
      *info = MAGMA_ERR_DEVICE_ALLOC;
      return *info;
    }

    /*     Form a Cholesky factorization of B. */
    magma_csetmatrix( n, n, b, ldb, db, lddb );

    magma_csetmatrix_async( n, n,
                            a,  lda,
                            da, ldda, stream );

    magma_cpotrf_gpu(uplo_[0], n, db, lddb, info);
    if (*info != 0) {
        *info = n + *info;
        return *info;
    }

    magma_queue_sync( stream );

    magma_cgetmatrix_async( n, n,
                            db, lddb,
                            b,  ldb, stream );

    /*  Transform problem to standard eigenvalue problem and solve. */
    magma_chegst_gpu(itype, uplo_[0], n, da, ldda, db, lddb, info);

    magma_cheevd_gpu(jobz_[0], uplo_[0], n, da, ldda, w, a, lda,
                     work, lwork, rwork, lrwork, iwork, liwork, info);
    /* Computing MAX */
    //    d__1 = (float) lopt, d__2 = work[1].r;
    //    lopt = (magma_int_t) max(d__1,d__2);
    /* Computing MAX */
    //    d__1 = (float) lropt;
    //    lropt = (magma_int_t) max(d__1,rwork[1]);
    /* Computing MAX */
    //    d__1 = (float) liopt, d__2 = (floatreal) iwork[1];
    //    liopt = (magma_int_t) max(d__1,d__2);

    if (wantz && *info == 0)
      {
        /* Backtransform eigenvectors to the original problem. */
        if (itype == 1 || itype == 2)
          {
            /* For A*x=(lambda)*B*x and A*B*x=(lambda)*x;
               backtransform eigenvectors: x = inv(L)'*y or inv(U)*y */
            if (lower) {
                *(unsigned char *)trans = MagmaConjTrans;
            } else {
                *(unsigned char *)trans = MagmaNoTrans;
            }

            Mymagma_ctrsm(MagmaLeft, uplo_[0], *trans, MagmaNonUnit,
                          n, n, c_one, db, lddb, da, ldda);

        } else if (itype == 3)
          {
            /*  For B*A*x=(lambda)*x;
                backtransform eigenvectors: x = L*y or U'*y */
            if (lower) {
                *(unsigned char *)trans = MagmaNoTrans;
            } else {
                *(unsigned char *)trans = MagmaConjTrans;
            }

            Mymagma_ctrmm(MagmaLeft, uplo_[0], *trans, MagmaNonUnit,
                          n, n, c_one, db, lddb, da, ldda);
        }

        magma_cgetmatrix( n, n, da, ldda, a, lda );

    }

    magma_queue_sync( stream );
    magma_queue_destroy( stream );

    /*work[0].r = (floatreal) lopt, work[0].i = 0.;
    rwork[0] = (floatreal) lropt;
    iwork[0] = liopt;*/

    magma_free( da );
    magma_free( db );

    return *info;
} /* magma_chegvd */
