/*    
    -- MAGMA (version 1.2.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       June 2012

       @author Stan Tomov
       @author Raffaele Solca

       @generated s Thu Jun 28 12:30:52 2012

*/
#include "common_magma.h"

extern "C" magma_int_t 
magma_ssyevd_gpu(char jobz, char uplo, 
                 magma_int_t n, 
                 float *da, magma_int_t ldda, 
                 float *w, 
                 float *wa,  magma_int_t ldwa,
                 float *work, magma_int_t lwork,
                 magma_int_t *iwork, magma_int_t liwork,
                 magma_int_t *info)
{
/*  -- MAGMA (version 1.2.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       June 2012

    Purpose   
    =======
    SSYEVD_GPU computes all eigenvalues and, optionally, eigenvectors of
    a real symmetric matrix A.  If eigenvectors are desired, it uses a   
    divide and conquer algorithm.   

    The divide and conquer algorithm makes very mild assumptions about   
    floating point arithmetic. It will work on machines with a guard   
    digit in add/subtract, or on those binary machines without guard   
    digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or   
    Cray-2. It could conceivably fail on hexadecimal or decimal machines   
    without guard digits, but we know of none.   

    Arguments   
    =========   
    JOBZ    (input) CHARACTER*1   
            = 'N':  Compute eigenvalues only;   
            = 'V':  Compute eigenvalues and eigenvectors.   

    UPLO    (input) CHARACTER*1   
            = 'U':  Upper triangle of A is stored;   
            = 'L':  Lower triangle of A is stored.   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    DA      (device input/output) REAL array on the GPU, 
            dimension (LDDA, N).
            On entry, the symmetric matrix A.  If UPLO = 'U', the   
            leading N-by-N upper triangular part of A contains the   
            upper triangular part of the matrix A.  If UPLO = 'L',   
            the leading N-by-N lower triangular part of A contains   
            the lower triangular part of the matrix A.   
            On exit, if JOBZ = 'V', then if INFO = 0, A contains the   
            orthonormal eigenvectors of the matrix A.   
            If JOBZ = 'N', then on exit the lower triangle (if UPLO='L')   
            or the upper triangle (if UPLO='U') of A, including the   
            diagonal, is destroyed.   

    LDDA    (input) INTEGER   
            The leading dimension of the array DA.  LDDA >= max(1,N).   

    W       (output) DOUBLE PRECISION array, dimension (N)   
            If INFO = 0, the eigenvalues in ascending order.

    WA      (workspace) DOUBLE PRECISION array, dimension (LDWA, N)

    LDWA    (input) INTEGER
            The leading dimension of the array WA.  LDWA >= max(1,N).

    WORK    (workspace/output) REAL array, dimension (MAX(1,LWORK))   
            On exit, if INFO = 0, WORK(1) returns the optimal LWORK.   

    LWORK   (input) INTEGER   
            The length of the array WORK.   
            If N <= 1,                LWORK must be at least 1.   
            If JOBZ  = 'N' and N > 1, LWORK must be at least 2*N + 1.   
            If JOBZ  = 'V' and N > 1, LWORK must be at least 1 + 6*N + 2*N**2.   

            If LWORK = -1, then a workspace query is assumed; the routine
            only calculates the optimal sizes of the WORK and IWORK
            arrays, returns these values as the first entries of the WORK
            and IWORK arrays, and no error message related to LWORK or
            LIWORK is issued by XERBLA.

    IWORK   (workspace/output) INTEGER array, dimension (MAX(1,LIWORK))   
            On exit, if INFO = 0, IWORK(1) returns the optimal LIWORK.   

    LIWORK  (input) INTEGER   
            The dimension of the array IWORK.   
            If N <= 1,                LIWORK must be at least 1.   
            If JOBZ  = 'N' and N > 1, LIWORK must be at least 1.   
            If JOBZ  = 'V' and N > 1, LIWORK must be at least 3 + 5*N.   

            If LIWORK = -1, then a workspace query is assumed; the
            routine only calculates the optimal sizes of the WORK and
            IWORK arrays, returns these values as the first entries of
            the WORK and IWORK arrays, and no error message related to
            LWORK or LIWORK is issued by XERBLA.

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   
            > 0:  if INFO = i and JOBZ = 'N', then the algorithm failed   
                  to converge; i off-diagonal elements of an intermediate   
                  tridiagonal form did not converge to zero;   
                  if INFO = i and JOBZ = 'V', then the algorithm failed   
                  to compute an eigenvalue while working on the submatrix   
                  lying in rows and columns INFO/(N+1) through   
                  mod(INFO,N+1).   

    Further Details   
    ===============   
    Based on contributions by   
       Jeff Rutter, Computer Science Division, University of California   
       at Berkeley, USA   

    Modified description of INFO. Sven, 16 Feb 05.   
    =====================================================================   */

    char uplo_[2] = {uplo, 0};
    char jobz_[2] = {jobz, 0};
    magma_int_t c__1 = 1;
    
    float d__1;

    float eps;
    magma_int_t inde;
    float anrm;
    float rmin, rmax;
    magma_int_t lopt;
    float sigma;
    magma_int_t iinfo, lwmin, liopt;
    magma_int_t lower;
    magma_int_t wantz;
    magma_int_t indwk2, llwrk2;
    magma_int_t iscale;
    float safmin;
    float bignum;
    magma_int_t indtau;
    magma_int_t indwrk, liwmin;
    magma_int_t llwork;
    float smlnum;
    magma_int_t lquery;

    bool dc_freed = false;

    float *dwork;
    float *dc;
    magma_int_t lddc = ldda;

    wantz = lapackf77_lsame(jobz_, MagmaVectorsStr);
    lower = lapackf77_lsame(uplo_, MagmaLowerStr);
    lquery = lwork == -1 || liwork == -1;

    *info = 0;
    if (! (wantz || lapackf77_lsame(jobz_, MagmaNoVectorsStr))) {
        *info = -1;
    } else if (! (lower || lapackf77_lsame(uplo_, MagmaUpperStr))) {
        *info = -2;
    } else if (n < 0) {
        *info = -3;
    } else if (ldda < max(1,n)) {
        *info = -5;
    }

    magma_int_t nb = magma_get_ssytrd_nb(n);
    
    if (wantz) {
      lwmin = 2 *n + n*n;
      liwmin = 5 * n + 3;
    } else {
      lwmin = n * (nb + 1);
      liwmin = 1;
    }

    work[0]  = lwmin;
    iwork[0] = liwmin; 

    if ((lwork < lwmin) && !lquery) {
        *info = -10;
    } else if ((liwork < liwmin) && ! lquery) {
        *info = -12;
    } else if ((liwork < liwmin) && ! lquery) {
      *info = -14;
    }

    if (*info != 0) {
        magma_xerbla( __func__, -(*info) );
        return MAGMA_ERR_ILLEGAL_VALUE;
    }
    else if (lquery) {
        return MAGMA_SUCCESS;
    }

    /* Quick return if possible */
    if (n == 0) {
        return MAGMA_SUCCESS;
    }

    if (n == 1) {
        float tmp;
        magma_sgetvector( 1, da, 1, &tmp, 1 );
        w[0] = tmp;
        if (wantz) {
            tmp = 1.;
            magma_ssetvector( 1, &tmp, 1, da, 1 );
        }
        return MAGMA_SUCCESS;
    }

    cudaStream_t stream;
    magma_queue_create( &stream );

    if (MAGMA_SUCCESS != magma_smalloc( &dc, n*lddc )) {
      fprintf(stderr, "!!!! device memory allocation error (magma_ssyevd_gpu)\n");
      return MAGMA_ERR_DEVICE_ALLOC;
    }
    if (MAGMA_SUCCESS != magma_smalloc( &dwork, n )) {
      fprintf(stderr, "!!!! device memory allocation error (magma_ssyevd_gpu)\n");
      return MAGMA_ERR_DEVICE_ALLOC;
    }

    --w;
    --work;
    --iwork;

    /* Get machine constants. */
    safmin = lapackf77_slamch("Safe minimum");
    eps = lapackf77_slamch("Precision");
    smlnum = safmin / eps;
    bignum = 1. / smlnum;
    rmin = magma_ssqrt(smlnum);
    rmax = magma_ssqrt(bignum);

    /* Scale matrix to allowable range, if necessary. */
    anrm = magmablas_slansy('M', uplo, n, da, ldda, dwork);
    magma_free( dwork );
    iscale = 0;
    if (anrm > 0. && anrm < rmin) {
        iscale = 1;
        sigma = rmin / anrm;
    } else if (anrm > rmax) {
        iscale = 1;
        sigma = rmax / anrm;
    }
    if (iscale == 1)
      magmablas_slascl(uplo, 0, 0, 1., sigma, n, n, da, ldda, info);

    /* Call SSYTRD to reduce symmetric matrix to tridiagonal form. */
    inde = 1;
    indtau = inde + n;
    indwrk = indtau + n;
    llwork = lwork - indwrk + 1;
    indwk2 = indwrk + n * n;
    llwrk2 = lwork - indwk2 + 1;
  
//#define ENABLE_TIMER
#ifdef ENABLE_TIMER 
    magma_timestr_t start, end;
    
    start = get_current_time();
#endif

#ifdef FAST_SYMV
    magma_ssytrd2_gpu(uplo, n, da, ldda, &w[1], &work[inde],
                      &work[indtau], wa, ldwa, &work[indwrk], llwork, 
                      dc, lddc*n, &iinfo);
#else
    magma_ssytrd_gpu(uplo, n, da, ldda, &w[1], &work[inde],
                     &work[indtau], wa, ldwa, &work[indwrk], llwork, 
                     &iinfo);
#endif

#ifdef ENABLE_TIMER    
    end = get_current_time();
    
    printf("time ssytrd = %6.2f\n", GetTimerValue(start,end)/1000.);
#endif        

    /* For eigenvalues only, call SSTERF.  For eigenvectors, first call   
       SSTEDC to generate the eigenvector matrix, WORK(INDWRK), of the   
       tridiagonal matrix, then call SORMTR to multiply it to the Householder 
       transformations represented as Householder vectors in A. */
    if (! wantz) {
        lapackf77_ssterf(&n, &w[1], &work[inde], info);
    } else {

#ifdef ENABLE_TIMER
        start = get_current_time();
#endif
        
        if (MAGMA_SUCCESS != magma_smalloc( &dwork, 3*n*(n/2 + 1) )) {
            magma_free( dc );  // if not enough memory is available free dc to be able do allocate dwork
            dc_freed=true;
#ifdef ENABLE_TIMER
            printf("dc deallocated\n");
#endif
            if (MAGMA_SUCCESS != magma_smalloc( &dwork, 3*n*(n/2 + 1) )) {
                fprintf (stderr, "!!!! device memory allocation error (magma_ssyevd_gpu)\n");
                return MAGMA_ERR_DEVICE_ALLOC;
            }
        }
        
        magma_sstedx('A', n, 0., 0., 0, 0, &w[1], &work[inde],
                     &work[indwrk], n, &work[indwk2],
                     llwrk2, &iwork[1], liwork, dwork, info);
        
        magma_free( dwork );

#ifdef ENABLE_TIMER  
        end = get_current_time();
        
        printf("time sstedx = %6.2f\n", GetTimerValue(start,end)/1000.);
#endif

        if(dc_freed){
            dc_freed = false;
            if (MAGMA_SUCCESS != magma_smalloc( &dc, n*lddc )) {
                fprintf (stderr, "!!!! device memory allocation error (magma_ssyevd_gpu)\n");
                return MAGMA_ERR_DEVICE_ALLOC;
            }
        }

        magma_ssetmatrix( n, n, &work[indwrk], n, dc, lddc );
        
#ifdef ENABLE_TIMER  
        start = get_current_time();
#endif

        magma_sormtr_gpu(MagmaLeft, uplo, MagmaNoTrans, n, n, da, ldda, &work[indtau],
                         dc, lddc, wa, ldwa, &iinfo);
        
        magma_scopymatrix( n, n,
                           dc, lddc,
                           da, ldda );

#ifdef ENABLE_TIMER    
        end = get_current_time();
        
        printf("time sormtr + copy = %6.2f\n", GetTimerValue(start,end)/1000.);
#endif        

    }

    /* If matrix was scaled, then rescale eigenvalues appropriately. */
    if (iscale == 1) {
        d__1 = 1. / sigma;
        blasf77_sscal(&n, &d__1, &w[1], &c__1);
    }

    MAGMA_S_SET2REAL(work[1], (float) lopt);
    iwork[1] = liopt;

    magma_queue_destroy( stream );
    if (!dc_freed)
        magma_free( dc );

    return MAGMA_SUCCESS;
} /* magma_ssyevd_gpu */