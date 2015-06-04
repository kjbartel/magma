/*
    -- MAGMA (version 1.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2011

       @author Stan Tomov
       @author Raffaele Solca

       @precisions normal d -> s

*/ 
#include "common_magma.h"

extern "C" magma_int_t
magma_dsygvd(magma_int_t itype, char jobz, char uplo, magma_int_t n,
             double *a, magma_int_t lda, double *b, magma_int_t ldb, 
             double *w, double *work, magma_int_t lwork, 
             magma_int_t *iwork, magma_int_t liwork, magma_int_t *info)
{
/*  -- MAGMA (version 1.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2011

    Purpose   
    =======   
    DSYGVD computes all the eigenvalues, and optionally, the eigenvectors   
    of a real generalized symmetric-definite eigenproblem, of the form   
    A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.  Here A and   
    B are assumed to be symmetric and B is also positive definite.   
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
            On entry, the symmetric matrix A.  If UPLO = 'U', the   
            leading N-by-N upper triangular part of A contains the   
            upper triangular part of the matrix A.  If UPLO = 'L',   
            the leading N-by-N lower triangular part of A contains   
            the lower triangular part of the matrix A.   

            On exit, if JOBZ = 'V', then if INFO = 0, A contains the   
            matrix Z of eigenvectors.  The eigenvectors are normalized   
            as follows:   
            if ITYPE = 1 or 2, Z\*\*T*B*Z = I;   
            if ITYPE = 3, Z\*\*T*inv(B)*Z = I.   
            If JOBZ = 'N', then on exit the upper triangle (if UPLO='U')   
            or the lower triangle (if UPLO='L') of A, including the   
            diagonal, is destroyed.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= max(1,N).   

    B       (input/output) COMPLEX*16 array, dimension (LDB, N)   
            On entry, the symmetric matrix B.  If UPLO = 'U', the   
            leading N-by-N upper triangular part of B contains the   
            upper triangular part of the matrix B.  If UPLO = 'L',   
            the leading N-by-N lower triangular part of B contains   
            the lower triangular part of the matrix B.   

            On exit, if INFO <= N, the part of B containing the matrix is   
            overwritten by the triangular factor U or L from the Cholesky   
            factorization B = U\*\*T*U or B = L*L\*\*T.   

    LDB     (input) INTEGER   
            The leading dimension of the array B.  LDB >= max(1,N).   

    W       (output) DOUBLE PRECISION array, dimension (N)   
            If INFO = 0, the eigenvalues in ascending order.   

    WORK    (workspace/output) COMPLEX*16 array, dimension (MAX(1,LWORK))   
            On exit, if INFO = 0, WORK(1) returns the optimal LWORK.   

    LWORK   (input) INTEGER   
            The length of the array WORK.   
            If N <= 1,                LWORK >= 1.   
            If JOBZ  = 'N' and N > 1, LWORK >= 2*N*nb + 1.   
            If JOBZ  = 'V' and N > 1, LWORK >= 1 + 6*N*nb + 2*N**2.   

            If LWORK = -1, then a workspace query is assumed; the routine   
            only calculates the optimal sizes of the WORK and   
            IWORK arrays, returns these values as the first entries of   
            the WORK and IWORK arrays, and no error message   
            related to LWORK or LIWORK is issued by XERBLA.   

    IWORK   (workspace/output) INTEGER array, dimension (MAX(1,LIWORK))   
            On exit, if INFO = 0, IWORK(1) returns the optimal LIWORK.   

    LIWORK  (input) INTEGER   
            The dimension of the array IWORK.   
            If N <= 1,                LIWORK >= 1.   
            If JOBZ  = 'N' and N > 1, LIWORK >= 1.   
            If JOBZ  = 'V' and N > 1, LIWORK >= 3 + 5*N.   

            If LIWORK = -1, then a workspace query is assumed; the   
            routine only calculates the optimal sizes of the WORK   
            and IWORK arrays, returns these values as the first entries   
            of the WORK and IWORK arrays, and no error message   
            related to LWORK or LIWORK is issued by XERBLA.   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   
            > 0:  DPOTRF or DSYEVD returned an error code:   
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

    Modified so that no backsubstitution is performed if DSYEVD fails to   
    converge (NEIG in old code could be greater than N causing out of   
    bounds reference to A - reported by Ralf Meyer).  Also corrected the   
    description of INFO and the test on ITYPE. Sven, 16 Feb 05.   
    =====================================================================  */

    char uplo_[2] = {uplo, 0};
    char jobz_[2] = {jobz, 0};

    static double zone = MAGMA_D_ONE;
    
    double *da;
    double *db;
    magma_int_t ldda = n;
    magma_int_t lddb = n;

    static magma_int_t lower;
    static char trans[1];
    static magma_int_t wantz, lquery;

    static magma_int_t lopt, lwmin, liopt, liwmin;
  
    static cudaStream_t stream;
    cudaStreamCreate(&stream);

    wantz = lapackf77_lsame(jobz_, MagmaVectorsStr);
    lower = lapackf77_lsame(uplo_, MagmaLowerStr);
    lquery = lwork == -1 || liwork == -1;

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

    magma_int_t nb = magma_get_dsytrd_nb(n); 
  
    if (n < 1) {
      liwmin = 1;
      lwmin = 1;
    } else if (wantz) {
      lwmin = 1 + 6 * n * nb + 2* n * n;
      liwmin = 5 * n + 3;
    } else {
        lwmin = 2 * n * nb + 1;
        liwmin = 1;
    }

    lopt = lwmin;
    liopt = liwmin;

    work[ 0] =  lopt;
    iwork[0] = liopt;

    if (lwork < lwmin && ! lquery) {
        *info = -11;
    } else if (liwork < liwmin && ! lquery) {
         *info = -13;
    }

    if (*info != 0) {
        magma_xerbla( __func__, -(*info) );
        return MAGMA_ERR_ILLEGAL_VALUE;
    }
    else if (lquery) {
        return MAGMA_SUCCESS;
    }

    /*  Quick return if possible */
    if (n == 0) {
        return 0;
    }

    if (cudaSuccess != cudaMalloc( (void**)&da, n*ldda*sizeof(double) ) ||
        cudaSuccess != cudaMalloc( (void**)&db, n*lddb*sizeof(double) ) ) {
      *info = -17;
      return MAGMA_ERR_CUBLASALLOC;
    }
  
    /* Form a Cholesky factorization of B. */
    cublasSetMatrix(n, n, sizeof(double), b, ldb, db, lddb);

    cudaMemcpy2DAsync(da, ldda*sizeof(double),
                       a,  lda*sizeof(double),
                      sizeof(double)*n, n,
                      cudaMemcpyHostToDevice, stream);  
  
    magma_dpotrf_gpu(uplo_[0], n, db, lddb, info);
    if (*info != 0) {
        *info = n + *info;
        return 0;
    }

    cudaStreamSynchronize(stream);
  
    cudaMemcpy2DAsync( b,  ldb*sizeof(double),
                      db, lddb*sizeof(double),
                      sizeof(double)*n, n,
                      cudaMemcpyDeviceToHost, stream);

    /*  Transform problem to standard eigenvalue problem and solve. */
    magma_dsygst_gpu(itype, uplo_[0], n, da, ldda, db, lddb, info);
  
    magma_dsyevd_gpu(jobz_[0], uplo_[0], n, da, ldda, w, a, lda, 
                     work, lwork, iwork, liwork, info);

    lopt  = max( lopt, (magma_int_t) work[0]);
    liopt = max(liopt, iwork[0]);

    if (wantz && *info == 0) 
      {
        /* Backtransform eigenvectors to the original problem. */
        if (itype == 1 || itype == 2) 
          {
            /* For A*x=(lambda)*B*x and A*B*x=(lambda)*x;   
               backtransform eigenvectors: x = inv(L)'*y or inv(U)*y */
            if (lower) {
                *(unsigned char *)trans = MagmaTrans;
            } else {
                *(unsigned char *)trans = MagmaNoTrans;
            }

            cublasDtrsm(MagmaLeft, uplo_[0], *trans, MagmaNonUnit,
                        n, n, zone, db, lddb, da, ldda);

        } else if (itype == 3) 
          {
            /*  For B*A*x=(lambda)*x;   
                backtransform eigenvectors: x = L*y or U'*y */
            if (lower) {
                *(unsigned char *)trans = MagmaNoTrans;
            } else {
                *(unsigned char *)trans = MagmaTrans;
            }

            cublasDtrmm(MagmaLeft, uplo_[0], *trans, MagmaNonUnit, 
                        n, n, zone, db, lddb, da, ldda);
        }

        cublasGetMatrix(n, n, sizeof(double), da, ldda, a, lda);

    }

    cudaStreamSynchronize(stream);
    cudaStreamDestroy(stream);
  
    work[0] = (double) lopt;
    iwork[0] = liopt;

    cublasFree(da);
    cublasFree(db);
  
    return MAGMA_SUCCESS;
} /* magma_dsygvd */
