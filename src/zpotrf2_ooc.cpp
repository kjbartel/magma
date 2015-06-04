/*
    -- MAGMA (version 1.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2011

       @precisions normal z -> s d c

*/
#include "common_magma.h"

/* === Define what BLAS to use ============================================ */
#define PRECISION_z
#if (defined(PRECISION_s) || defined(PRECISION_d))
  #define cublasDgemm magmablas_dgemm
  #define cublasDtrsm magmablas_dtrsm
#endif

#if (GPUSHMEM >= 200)
#if (defined(PRECISION_s))
     #undef  cublasSgemm
     #define cublasSgemm magmablas_sgemm_fermi80
  #endif
#endif
/* === End defining what BLAS to use ====================================== */
// Flops formula
#include "../testing/flops.h"
#define PRECISION_z
#if defined(PRECISION_z) || defined(PRECISION_c)
#define FLOPS(n) ( 6. * FMULS_POTRF(n) + 2. * FADDS_POTRF(n) )
#else
#define FLOPS(n) (      FMULS_POTRF(n) +      FADDS_POTRF(n) )
#endif

extern "C" magma_int_t 
magma_zhtodpo(int num_gpus, char *uplo, magma_int_t m, magma_int_t n, magma_int_t off_i, magma_int_t off_j, magma_int_t nb,
              cuDoubleComplex *h_A, magma_int_t lda, cuDoubleComplex **d_lA, magma_int_t ldda, cudaStream_t **stream,
              magma_int_t *info);

extern "C" magma_int_t 
magma_zdtohpo(int num_gpus, char *uplo, magma_int_t m, magma_int_t n, magma_int_t off_i, magma_int_t off_j, magma_int_t nb, magma_int_t NB,
              cuDoubleComplex *a, magma_int_t lda, cuDoubleComplex **work, magma_int_t ldda, cudaStream_t **stream,
              magma_int_t *info);

extern "C" magma_int_t
magma_zpotrf2_mgpu(int num_gpus, char uplo, magma_int_t m, magma_int_t n, magma_int_t off_i, magma_int_t off_j, magma_int_t nb,
                   cuDoubleComplex **d_lA, magma_int_t ldda, cuDoubleComplex **d_lP, magma_int_t lddlp, 
                   cuDoubleComplex *work, magma_int_t ldwrk, cudaStream_t **streaml, magma_int_t *info);

#define A(i, j)  (a   +(j)*lda  + (i))
#define dA(d, i, j) (dwork[(d)]+(j)*lddla + (i))
#define dT(d, i, j) (dt[(d)]   +(j)*ldda  + (i))
#define dAup(d, i, j) (dwork[(d)]+(j)*NB + (i))
#define dTup(d, i, j) (dt[(d)]   +(j)*nb + (i))

extern "C" magma_int_t 
magma_zpotrf2_ooc(magma_int_t num_gpus0, char uplo, magma_int_t n, 
                  cuDoubleComplex *a, magma_int_t lda, magma_int_t *info)
{
/*  -- MAGMA (version 1.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2011

    Purpose   
    =======   

    ZPOTRF_OOC computes the Cholesky factorization of a complex Hermitian   
    positive definite matrix A. This version does not require work
    space on the GPU passed as input. GPU memory is allocated in the
    routine. The matrix A may not fit entirely in the GPU memory.

    The factorization has the form   
       A = U\*\*H * U,  if UPLO = 'U', or   
       A = L  * L\*\*H, if UPLO = 'L',   
    where U is an upper triangular matrix and L is lower triangular.   

    This is the block version of the algorithm, calling Level 3 BLAS.   

    Arguments   
    =========   

    UPLO    (input) CHARACTER*1   
            = 'U':  Upper triangle of A is stored;   
            = 'L':  Lower triangle of A is stored.   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    A       (input/output) COMPLEX_16 array, dimension (LDA,N)   
            On entry, the symmetric matrix A.  If UPLO = 'U', the leading   
            N-by-N upper triangular part of A contains the upper   
            triangular part of the matrix A, and the strictly lower   
            triangular part of A is not referenced.  If UPLO = 'L', the   
            leading N-by-N lower triangular part of A contains the lower   
            triangular part of the matrix A, and the strictly upper   
            triangular part of A is not referenced.   

            On exit, if INFO = 0, the factor U or L from the Cholesky   
            factorization A = U\*\*H*U or A = L*L\*\*H.   

            Higher performance is achieved if A is in pinned memory, e.g.
            allocated using cudaMallocHost.

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= max(1,N).   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value 
                  if INFO = -6, the GPU memory allocation failed 
            > 0:  if INFO = i, the leading minor of order i is not   
                  positive definite, and the factorization could not be   
                  completed.   

    =====================================================================    */


    /* Local variables */
    cuDoubleComplex        zone  = MAGMA_Z_ONE;
    cuDoubleComplex        mzone = MAGMA_Z_NEG_ONE;
    cuDoubleComplex        *dwork[4], *dt[4], *work;

    char                uplo_[2] = {uplo, 0};
    magma_int_t            ldda, lddla, ldwrk, nb, iinfo, n_local[4], J2, d, num_gpus;
    static magma_int_t    j, jj, jb, jb1, jb2, jb3, J, JB, NB, MB;
    double                done  = (double) 1.0;
    double                mdone = (double)-1.0;
    long int            upper = lapackf77_lsame(uplo_, "U");
#if CUDA_VERSION > 3010
    size_t totalMem;
#else
    unsigned int totalMem;
#endif
    CUdevice dev;
    static cudaStream_t stream[4][4];
    magma_timestr_t start, end;
    double chol_time = 0.0;

    *info = 0;
    if ((! upper) && (! lapackf77_lsame(uplo_, "L"))) {
      *info = -1;
    } else if (n < 0) {
      *info = -2;
    } else if (lda < max(1,n)) {
      *info = -4;
    }
    if (*info != 0) {
        magma_xerbla( __func__, -(*info) );
        return MAGMA_ERR_ILLEGAL_VALUE;
    }

    /* Quick return */
    if ( n == 0 )
      return MAGMA_SUCCESS;

    nb = magma_get_dpotrf_nb(n);
    if( num_gpus0 > n/nb ) {
      num_gpus = n/nb;
      if( n%nb != 0 ) num_gpus ++;
    } else {
      num_gpus = num_gpus0;
    }
    ldda = n/(nb*num_gpus);
    if( n%(nb*num_gpus) != 0 ) ldda++;
    ldda = num_gpus*((nb*ldda+31)/32)*32;

    /* figure out NB */
    cuDeviceGet( &dev, 0);
    cuDeviceTotalMem( &totalMem, dev );
    totalMem /= sizeof(cuDoubleComplex);
    MB = n;  /* number of rows in the big panel    */
    NB = (magma_int_t)(num_gpus*(0.8*totalMem/ldda-nb)); /* number of columns in the big panel */
    if( NB >= n ) {
#ifdef CHECK_ZPOTRF_OOC
      printf( "      * still fit in GPU memory.\n" );
#endif
      NB = n;
    } else {
#ifdef CHECK_ZPOTRF_OOC
      printf( "      * don't fit in GPU memory.\n" );
#endif
      NB = (NB / nb) * nb;   /* making sure it's devisable by nb   */
    }
#ifdef CHECK_ZPOTRF_OOC
    if( NB != n ) printf( "      * running in out-core mode (n=%d, NB=%d, nb=%d).\n",n,NB,nb );
    else          printf( "      * running in in-core mode  (n=%d, NB=%d, nb=%d).\n",n,NB,nb );
    fflush(stdout);
#endif
    ldda  = ((n+31)/32)*32;
    lddla = ((nb*(1+n/(nb*num_gpus))+31)/32)*32;
    for (d=0; d<num_gpus; d++ ) {
      cudaSetDevice(d);
      if (CUBLAS_STATUS_SUCCESS != cublasAlloc(NB*lddla+nb*ldda, sizeof(cuDoubleComplex), (void**)&dt[d])) {
        *info = -6;
        return MAGMA_ERR_CUBLASALLOC;
      }
      dwork[d] = &dt[d][nb*ldda];
      if ( (cudaSuccess != cudaStreamCreate(&stream[d][0])) ||
           (cudaSuccess != cudaStreamCreate(&stream[d][1])) ||
           (cudaSuccess != cudaStreamCreate(&stream[d][2])) ||
           (cudaSuccess != cudaStreamCreate(&stream[d][3])) ) {
        *info = -6;
        return cudaErrorInvalidValue;
      }
    }
    cudaSetDevice(0);
    ldwrk = n;
    if (cudaSuccess != cudaMallocHost( (void**)&work, ldwrk*nb*sizeof(cuDoubleComplex) ) ) {
      *info = -6;
      return MAGMA_ERR_HOSTALLOC;
    }

    if (nb <= 1 || nb >= n) {
      lapackf77_zpotrf(uplo_, &n, a, &lda, info);
    } else {

    /* Use hybrid blocked code. */
    if (upper) {
      /* =========================================================== *
       * Compute the Cholesky factorization A = U'*U.                *
       * big panel is divided by block-row and distributed in block  *
       * column cyclic format                                        */

      /* for each big-panel */
      for( J=0; J<n; J+=NB ) {
        JB = min(NB,n-J);
        jb = min(JB,nb);
        if( num_gpus0 > (n-J)/nb ) {
          num_gpus = (n-J)/nb;
          if( (n-J)%nb != 0 ) num_gpus ++;
        } else {
          num_gpus = num_gpus0;
        }

        /* load the new big-panel by block-rows */
        magma_zhtodpo( num_gpus, &uplo, JB, n, J, J, nb, a, lda, dwork, NB, (cudaStream_t **)stream, &iinfo);

#ifdef  ROW_MAJOR_PROFILE
        start = get_current_time();
#endif
        /* update with the previous big-panels */
        for( j=0; j<J; j+=nb ) {
          /* upload the diagonal of big panel */
          for( d=0; d<num_gpus; d++ ) {
            cudaSetDevice(d);
            cudaMemcpy2DAsync( dTup(d, 0, J), nb *sizeof(cuDoubleComplex), 
                               A(j, J),       lda*sizeof(cuDoubleComplex), 
                               sizeof(cuDoubleComplex)*nb, JB, 
                               cudaMemcpyHostToDevice,stream[d][0]);
            n_local[d] = 0;
          }

          /* upload off-diagonals */
          for( jj=J+JB; jj<n; jj+=nb ) {
            d  = ((jj-J)/nb)%num_gpus;
            cudaSetDevice(d);

            jb2 = min(nb, n-jj);
            cudaMemcpy2DAsync( dTup(d, 0, J+JB+n_local[d]), nb *sizeof(cuDoubleComplex), 
                               A(j, jj),                    lda*sizeof(cuDoubleComplex), 
                               sizeof(cuDoubleComplex)*nb, jb2, 
                               cudaMemcpyHostToDevice,stream[d][0]);
            n_local[d] += jb2;
          }
          for( d=0; d<num_gpus; d++ ) {
              cudaSetDevice(d);
            cudaStreamSynchronize(stream[d][1]);
          }

          /* update the current big-panel using the previous block-row */
          for( jj=0; jj<JB; jj+=nb ) { /* diagonal */
            d = (jj/nb)%num_gpus;
            cudaSetDevice(d);

            J2 = (jj/(nb*num_gpus))*nb;
            jb1 = min(JB,jj+nb); 
            jb2 = min(nb,JB-jj); 
            jb3 = min(nb,J-j);
            jb  = jb1-jb2;
            cublasZgemm( MagmaConjTrans, MagmaNoTrans, 
                         jb, jb2, nb, 
                         mzone, dTup(d, 0, J   ),  nb, 
                                dTup(d, 0, J+jb),  nb,
                         zone,  dAup(d, 0, J2), NB);

            cublasZherk(MagmaUpper, MagmaConjTrans, jb2, jb3,
                        mdone, dTup(d, 0,  J+jb),  nb,
                        done,  dAup(d, jb, J2), NB);
          }

          if( n > J+JB ) { /* off-diagonal */
            for( d=0; d<num_gpus; d++ ) {
              cudaSetDevice(d);
              /* local number of columns in the big panel */
              n_local[d] = (((n-J)/nb)/num_gpus)*nb;
              if (d < ((n-J)/nb)%num_gpus)
                n_local[d] += nb;
              else if (d == ((n-J)/nb)%num_gpus)
                n_local[d] += (n-J)%nb;

              /* local number of columns in diagonal */
              n_local[d] -= ((JB/nb)/num_gpus)*nb;
              if (d < (JB/nb)%num_gpus)
                n_local[d] -= nb;

              J2 = nb*(JB/(nb*num_gpus));
              if( d < (JB/nb)%num_gpus ) J2+=nb;

              cublasZgemm( MagmaConjTrans, MagmaNoTrans, 
                           JB, n_local[d], nb, 
                           mzone, dTup(d, 0, J   ),  nb, 
                                  dTup(d, 0, J+JB),  nb,
                           zone,  dAup(d, 0, J2), NB);
            }
          } 
        } /* end of updates with previous rows */

        /* factor the big panel */
        magma_zpotrf2_mgpu(num_gpus, uplo, JB, n-J, J, J, nb, dwork, NB, dt, ldda, a, lda, (cudaStream_t **)stream, &iinfo);
        if( iinfo != 0 ) {
            *info = J+iinfo;
            break;
        }
#ifdef  ROW_MAJOR_PROFILE
        end = get_current_time();
        chol_time += GetTimerValue(start, end);
#endif

        /* upload the off-diagonal (and diagonal!!!) big panel */
        //if( J+JB < n )
        //cudaMemcpy2DAsync(  A  (   J, J+JB), lda*sizeof(cuDoubleComplex),
        //                   dAup(d, 0, J+JB), NB *sizeof(cuDoubleComplex),
        //                   sizeof(cuDoubleComplex)*JB, n-J-JB, 
        //                   cudaMemcpyDeviceToHost,stream[2]);
        //magma_zdtohpo(num_gpus, &uplo, JB, n, J, J, nb, a, lda, dwork, NB, (cudaStream_t **)stream, &iinfo);
        magma_zdtohpo(num_gpus, &uplo, JB, n, J, J, nb, NB, a, lda, dwork, NB, (cudaStream_t **)stream, &iinfo);

      }
    } else {
      /* ========================================================= *
       * Compute the Cholesky factorization A = L*L'.              */

      /* for each big-panel */
      for( J=0; J<n; J+=NB ) {
        JB = min(NB,n-J);
        if( num_gpus0 > (n-J)/nb ) {
          num_gpus = (n-J)/nb;
          if( (n-J)%nb != 0 ) num_gpus ++;
        } else {
          num_gpus = num_gpus0;
        }

        /* load the new big-panel by block-columns */
        magma_zhtodpo( num_gpus, &uplo, n, JB, J, J, nb, a, lda, dwork, lddla, (cudaStream_t **)stream, &iinfo);

        /* update with the previous big-panels */
#ifdef  ROW_MAJOR_PROFILE
        start = get_current_time();
#endif
        for( j=0; j<J; j+=nb ) {
          /* upload the diagonal of big panel */
          for( d=0; d<num_gpus; d++ ) {
            cudaSetDevice(d);
            cudaMemcpy2DAsync( dT(d, J, 0), ldda *sizeof(cuDoubleComplex), 
                               A(J, j),     lda  *sizeof(cuDoubleComplex), 
                               sizeof(cuDoubleComplex)*JB, nb, 
                               cudaMemcpyHostToDevice,stream[d][0]);
            n_local[d] = 0;
          }

          /* upload off-diagonals */
          for( jj=J+JB; jj<n; jj+=nb ) {
            d  = ((jj-J)/nb)%num_gpus;
            cudaSetDevice(d);

            jb2 = min(nb, n-jj);
            cudaMemcpy2DAsync( dT(d, J+JB+n_local[d], 0), ldda *sizeof(cuDoubleComplex), 
                               A(jj, j),                  lda*sizeof(cuDoubleComplex), 
                               sizeof(cuDoubleComplex)*jb2, nb, 
                               cudaMemcpyHostToDevice,stream[d][0]);
            n_local[d] += jb2;
          }
          for( d=0; d<num_gpus; d++ ) {
              cudaSetDevice(d);
            cudaStreamSynchronize(stream[d][1]);
          }

          /* update the current big-panel using the previous block-row */
          for( jj=0; jj<JB; jj+=nb ) { /* diagonal */
            d = (jj/nb)%num_gpus;
            cudaSetDevice(d);

            J2 = (jj/(nb*num_gpus))*nb;
            jb1 = min(JB,jj+nb); 
            jb2 = min(nb,JB-jj); 
            jb3 = min(nb,J-j);
            jb  = jb1-jb2;
            cublasZgemm( MagmaNoTrans, MagmaConjTrans, 
                         jb2, jb, nb, 
                         mzone, dT(d, J+jb, 0), ldda, 
                                dT(d, J,    0), ldda,
                         zone,  dA(d, J2,   0), lddla);

            cublasZherk(MagmaLower, MagmaNoTrans, jb2, jb3,
                        mdone, dT(d, J+jb, 0),  ldda,
                        done,  dA(d, J2,  jb ), lddla);
          }

          if( n > J+JB ) { /* off-diagonal */
            for( d=0; d<num_gpus; d++ ) {
              cudaSetDevice(d);
              /* local number of columns in the big panel */
              n_local[d] = (((n-J)/nb)/num_gpus)*nb;
              if (d < ((n-J)/nb)%num_gpus)
                n_local[d] += nb;
              else if (d == ((n-J)/nb)%num_gpus)
                n_local[d] += (n-J)%nb;

              /* local number of columns in diagonal */
              n_local[d] -= ((JB/nb)/num_gpus)*nb;
              if (d < (JB/nb)%num_gpus)
                n_local[d] -= nb;

              J2 = nb*(JB/(nb*num_gpus));
              if( d < (JB/nb)%num_gpus ) J2+=nb;

              cublasZgemm( MagmaNoTrans, MagmaConjTrans, 
                           n_local[d], JB, nb, 
                           mzone, dT(d, J+JB, 0), ldda, 
                                  dT(d, J,    0), ldda,
                           zone,  dA(d, J2,   0), lddla);
            }
          }
        }

        /* factor the big panel */
        magma_zpotrf2_mgpu(num_gpus, uplo, n-J, JB, J, J, nb, dwork, lddla, dt, ldda, a, lda, (cudaStream_t **)stream, &iinfo);
        if( iinfo != 0 ) {
            *info = J+iinfo;
            break;
        }
#ifdef  ROW_MAJOR_PROFILE
        end = get_current_time();
        chol_time += GetTimerValue(start, end);
#endif

        /* upload the off-diagonal big panel */
        //if( J+JB < n )
        //cudaMemcpy2DAsync(  A(J+JB, J),    lda *sizeof(cuDoubleComplex),
        //                   dA(d, J+JB, 0), ldda*sizeof(cuDoubleComplex),
        //                   sizeof(cuDoubleComplex)*(n-J-JB), JB,
        //                   cudaMemcpyDeviceToHost,stream[0][2]);
        magma_zdtohpo( num_gpus, &uplo, n, JB, J, J, nb, NB, a, lda, dwork, lddla, (cudaStream_t **)stream, &iinfo);

      } /* end of for J */
    } /* if upper */
    } /* if nb */

    if( num_gpus0 > n/nb ) {
      num_gpus = n/nb;
      if( n%nb != 0 ) num_gpus ++;
    } else {
      num_gpus = num_gpus0;
    }
    for (d=0; d<num_gpus; d++ ) {
        cudaSetDevice(d);

        cublasFree(dt[d]);
        cudaStreamDestroy(stream[d][0]);
        cudaStreamDestroy(stream[d][1]);
        cudaStreamDestroy(stream[d][2]);
        cudaStreamDestroy(stream[d][3]);
    }
    cudaSetDevice(0);
    cudaFreeHost(work);

#ifdef  ROW_MAJOR_PROFILE
    printf("\n Performance %f / %f = %f GFlop/s\n", FLOPS((double)n)/1000000, chol_time, FLOPS( (double)n ) / (1000000*chol_time));
#endif
    return MAGMA_SUCCESS;
} /* magma_zpotrf_ooc */

#undef A
#undef dA
#undef dT
#undef dAup
#undef dTup
