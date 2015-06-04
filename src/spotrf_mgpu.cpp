/*
    -- MAGMA (version 1.3.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2012

       @generated s Wed Nov 14 22:52:59 2012

*/
#include "common_magma.h"

/* use two streams; one for comm, and one for comp */
extern "C" magma_int_t
magma_spotrf2_mgpu(int num_gpus, char uplo, magma_int_t m, magma_int_t n, 
                   magma_int_t off_i, magma_int_t off_j, magma_int_t nb,
                   float **d_lA, magma_int_t ldda, 
                   float **d_lP, magma_int_t lddp,
                   float  *a,    magma_int_t lda, magma_int_t h,
                   cudaStream_t stream[][3], cudaEvent_t event[][5],
                   magma_int_t *info );

/* use three streams; seems to be faster on Keeneland, but has problem on Pluto */
extern "C" magma_int_t
magma_spotrf3_mgpu(int num_gpus, char uplo, magma_int_t m, magma_int_t n, 
                   magma_int_t off_i, magma_int_t off_j, magma_int_t nb,
                   float **d_lA, magma_int_t ldda, 
                   float **d_lP, magma_int_t lddp,
                   float  *a,    magma_int_t lda, magma_int_t h,
                   cudaStream_t stream[][3], cudaEvent_t event[][5],
                   magma_int_t *info );

extern "C" magma_int_t
magma_spotrf_mgpu(int num_gpus, char uplo, magma_int_t n, 
                  float **d_lA, magma_int_t ldda, magma_int_t *info)
{
/*  -- MAGMA (version 1.3.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2012

    Purpose   
    =======   
    SPOTRF computes the Cholesky factorization of a real symmetric   
    positive definite matrix dA.   

    The factorization has the form   
       dA = U**T * U,  if UPLO = 'U', or   
       dA = L  * L**T,  if UPLO = 'L',   
    where U is an upper triangular matrix and L is lower triangular.   

    This is the block version of the algorithm, calling Level 3 BLAS.   

    Arguments   
    =========   
    UPLO    (input) CHARACTER*1   
            = 'U':  Upper triangle of dA is stored;   
            = 'L':  Lower triangle of dA is stored.   

    N       (input) INTEGER   
            The order of the matrix dA.  N >= 0.   

    dA      (input/output) REAL array on the GPU, dimension (LDDA,N)   
            On entry, the symmetric matrix dA.  If UPLO = 'U', the leading   
            N-by-N upper triangular part of dA contains the upper   
            triangular part of the matrix dA, and the strictly lower   
            triangular part of dA is not referenced.  If UPLO = 'L', the   
            leading N-by-N lower triangular part of dA contains the lower   
            triangular part of the matrix dA, and the strictly upper   
            triangular part of dA is not referenced.   

            On exit, if INFO = 0, the factor U or L from the Cholesky   
            factorization dA = U**T * U or dA = L * L**T.   

    LDDA     (input) INTEGER   
            The leading dimension of the array dA.  LDDA >= max(1,N).
            To benefit from coalescent memory accesses LDDA must be
            dividable by 16.

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   
            > 0:  if INFO = i, the leading minor of order i is not   
                  positive definite, and the factorization could not be   
                  completed.   
    =====================================================================   */


    magma_int_t     j, jb, nb, nb0, nb2, d, id, j_local, j_local2, lddp, h;
    char            uplo_[2] = {uplo, 0};
    float c_one     = MAGMA_S_ONE;
    float c_neg_one = MAGMA_S_NEG_ONE;
    float *work;
    float          d_one     =  1.0;
    float          d_neg_one = -1.0;
    int upper = lapackf77_lsame(uplo_, "U");
    float *dwork[MagmaMaxGPUs];
    cudaStream_t    stream[MagmaMaxGPUs][3];
    cudaEvent_t     event[MagmaMaxGPUs][5];

    *info = 0;
    nb = magma_get_spotrf_nb(n);
    if ( (! upper) && (! lapackf77_lsame(uplo_, "L")) ) {
        *info = -1;
    } else if (n < 0) {
        *info = -2;
    } else if (!upper) {
        lddp = nb*(n/(nb*num_gpus));
        if( n%(nb*num_gpus) != 0 ) lddp+=min(nb,n-num_gpus*lddp);
        if( ldda < lddp ) *info = -4;
    } else if( ldda < n ) {
        *info = -4;
    }
    if (*info != 0) {
        magma_xerbla( __func__, -(*info) );
        return *info;
    }

    if (num_gpus == 1 && ((nb <= 1) || (nb >= n)) ) {
      /*  Use unblocked code. */
      magma_setdevice(0);
      if (MAGMA_SUCCESS != magma_smalloc_pinned( &work, n*nb )) {
        *info = MAGMA_ERR_HOST_ALLOC;
        return *info;
      }
      magma_sgetmatrix( n, n, d_lA[0], ldda, work, n );
      lapackf77_spotrf(uplo_, &n, work, &n, info);
      magma_ssetmatrix( n, n, work, n, d_lA[0], ldda );
      magma_free_pinned( work );
    } else {
      lddp = 32*((n+31)/32);
      for( d=0; d<num_gpus; d++ ) {
        magma_setdevice(d);
        if (MAGMA_SUCCESS != magma_smalloc( &dwork[d], num_gpus*nb*lddp )) {
          for( j=0; j<d; j++ ) {
            magma_setdevice(j);
            magma_free( dwork[j] );
          }
          *info = MAGMA_ERR_DEVICE_ALLOC;
          return *info;
        }
        for( j=0; j<3; j++ ) magma_queue_create( &stream[d][j] );
        for( j=0; j<5; j++ ) magma_event_create( &event[d][j]  );
      }
      magma_setdevice(0);
      h = 1; //num_gpus; //(n+nb-1)/nb;
      if (MAGMA_SUCCESS != magma_smalloc_pinned( &work, n*nb*h )) {
        *info = MAGMA_ERR_HOST_ALLOC;
        return *info;
      }
      if (upper) {
          /* with two streams */
          //magma_spotrf2_mgpu(num_gpus, uplo, n, n, 0, 0, nb, d_lA, ldda, dwork, lddp, work, n,  
          //                   h, stream, event, info);
          /* with three streams */
          magma_spotrf3_mgpu(num_gpus, uplo, n, n, 0, 0, nb, d_lA, ldda, dwork, lddp, work, n,  
                             h, stream, event, info);
      } else {
          /* with two streams */
          //magma_spotrf2_mgpu(num_gpus, uplo, n, n, 0, 0, nb, d_lA, ldda, dwork, lddp, work, nb*h, 
          //                   h, stream, event, info);
          /* with three streams */
          magma_spotrf3_mgpu(num_gpus, uplo, n, n, 0, 0, nb, d_lA, ldda, dwork, lddp, work, nb*h, 
                             h, stream, event, info);
      }

      /* clean up */
      for( d=0; d<num_gpus; d++ ) {
        magma_setdevice(d);
        for( j=0; j<3; j++ ) {
            magma_queue_sync( stream[d][j] );
            magma_queue_destroy( stream[d][j] ); 
        }
        magmablasSetKernelStream(NULL);

        for( j=0; j<5; j++ ) 
            magma_event_destroy( event[d][j] );

        magma_free( dwork[d] );
      }
      magma_setdevice(0);
      magma_free_pinned( work );
    } /* end of not lapack */


    return *info;
} /* magma_spotrf_mgpu */
