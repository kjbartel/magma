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

#define A(i, j)  (a   +(j)*lda  + (i))
#define dA(i, j) (work+(j)*ldda + (i))
#define dT(i, j) (dt  +(j)*ldda + (i))
#define dAup(i, j) (work+(j)*NB + (i))
#define dTup(i, j) (dt  +(j)*nb + (i))

extern "C" magma_int_t 
magma_zpotrf_ooc(char uplo, magma_int_t n, 
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
    cuDoubleComplex            zone  = MAGMA_Z_ONE;
    cuDoubleComplex            mzone = MAGMA_Z_NEG_ONE;
    cuDoubleComplex            *work, *dt;

    char                    uplo_[2] = {uplo, 0};
    magma_int_t                ldda, nb;
    static magma_int_t        j, jj, jb, J, JB, NB, MB;
    double                    done  = (double) 1.0;
    double                    mdone = (double)-1.0;
    long int                upper = lapackf77_lsame(uplo_, "U");
#if CUDA_VERSION > 3010
        size_t totalMem;
#else
        unsigned int totalMem;
#endif
    CUdevice dev;
    static cudaStream_t stream[3];

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

    ldda = ((n+31)/32)*32;
    
    cudaStreamCreate(&stream[0]);
    cudaStreamCreate(&stream[1]);
    cudaStreamCreate(&stream[2]);

    nb = magma_get_dpotrf_nb(n);
    /* figure out NB */
    cuDeviceGet( &dev, 0);
        cuDeviceTotalMem( &totalMem, dev );
        totalMem /= sizeof(cuDoubleComplex);
        MB = n;                                /* number of rows in the big panel    */
        NB = (magma_int_t)(0.8*totalMem/n-nb); /* number of columns in the big panel */
        if( NB >= n ) {
#ifdef CHECK_ZPOTRF_OOC
          printf( "      * still fit in GPU memory.\n" );
#endif
          NB = n;
        }
#ifdef CHECK_ZPOTRF_OOC
          else {
          printf( "      * don't fit in GPU memory.\n" );
        }
#endif
        NB = (NB / nb) * nb;   /* making sure it's devisable by nb   */
    if (CUBLAS_STATUS_SUCCESS != cublasAlloc((NB+nb)*ldda, sizeof(cuDoubleComplex), (void**)&dt)) {
          *info = -6;
          return MAGMA_ERR_CUBLASALLOC;
    }
        work = &dt[nb*ldda];
#ifdef CHECK_ZPOTRF_OOC
        if( NB != n ) printf( "      * running in out-core mode (n=%d, NB=%d, nb=%d).\n",n,NB,nb );
        else          printf( "      * running in in-core mode  (n=%d, NB=%d, nb=%d).\n",n,NB,nb );
        fflush(stdout);
#endif


    if (nb <= 1 || nb >= n) {
          lapackf77_zpotrf(uplo_, &n, a, &lda, info);
    } else {

        /* Use hybrid blocked code. */
        if (upper) {
          /* ========================================================= *
           * Compute the Cholesky factorization A = U'*U.              */

          /* for each big-panel */
          for( J=0; J<n; J+=NB ) {
                JB = min(NB,n-J);

                /* load the new big-panel by block-rows */
            for (jj=0; jj<JB; jj+=nb) {
                  j  = J+jj;
                  jb = min(nb, (n-j));
              cublasSetMatrix(jb, (n-j), sizeof(cuDoubleComplex), 
                              A(j, j), lda, dAup(jj,j), NB);
                }
                /* load the panel in one-shot */
        //jb = min(nb, (n-J));
        //cublasSetMatrix(JB, n-J, sizeof(cuDoubleComplex),
        //                A(J, J), lda, dAup(0,J), NB);

                /* update with the previous big-panels */
                for( j=0; j<J; j+=nb ) {
                  /* upload the block-rows */
              cublasSetMatrix(nb, (n-J), sizeof(cuDoubleComplex), 
                              A(j, J), lda, dTup(0, J), nb);

                  /* update the current big-panel *
                   * using the previous block-row */
              cublasZherk(MagmaUpper, MagmaConjTrans, JB, nb,
                          mdone, dTup(0, J), nb, 
                          done,  dAup(0, J), NB);
                  if( (J+JB) < n ) 
              cublasZgemm( MagmaConjTrans, MagmaNoTrans, 
                           JB, (n-J-JB), nb, 
                           mzone, dTup(0, J   ), nb, 
                                  dTup(0, J+JB), nb,
                           zone,  dAup(0, J+JB), NB);
                }

                /* for each block-column in the big panel */
            for (jj=0; jj<JB; jj+=nb) {
                  j  = J+jj;
              jb = min(nb, (n-j));

              /* Update the current diagonal block */
              cublasZherk(MagmaUpper, MagmaConjTrans, jb, jj, 
                          mdone, dAup(0,  j), NB, 
                          done,  dAup(jj, j), NB);

                  /* send the diagonal-block to CPU */
              cudaMemcpy2DAsync(  A  (J, j), lda*sizeof(cuDoubleComplex), 
                                 dAup(0, j), NB *sizeof(cuDoubleComplex), 
                                 sizeof(cuDoubleComplex)*(jj+jb), jb,
                                 cudaMemcpyDeviceToHost, stream[1]);
              //cudaMemcpy2DAsync(  A  ( j, j), lda*sizeof(cuDoubleComplex), 
              //                   dAup(jj, j), NB *sizeof(cuDoubleComplex), 
              //                   sizeof(cuDoubleComplex)*jb, jb,
              //                   cudaMemcpyDeviceToHost, stream[1]);
                
              if ( (j+jb) < n) {
                        /* update the current off-diagonal blocks with the previous rows */
                cublasZgemm(MagmaConjTrans, MagmaNoTrans, 
                            jb, (n-j-jb), jj,
                            mzone, dAup(0,  j   ), NB, 
                                   dAup(0,  j+jb), NB,
                            zone,  dAup(jj, j+jb), NB);
              }
             
                  /* factor the diagonal block */
                  cudaStreamSynchronize(stream[1]);
                  lapackf77_zpotrf(MagmaUpperStr, &jb, A(j, j), &lda, info);
                  if (*info != 0) {
                    *info = *info + j;
                    break;
                  }

              if ( (j+jb) < n ) {
                    /* send the diagonal block to GPU */
                    cudaMemcpy2DAsync(dAup(jj, j), NB  * sizeof(cuDoubleComplex), 
                                   A  (j,  j), lda * sizeof(cuDoubleComplex), 
                                   sizeof(cuDoubleComplex)*jb, jb, 
                                   cudaMemcpyHostToDevice,stream[0]);

                        /* do the solves on GPU */
                cublasZtrsm(MagmaLeft, MagmaUpper, MagmaConjTrans, MagmaNonUnit, 
                            jb, (n-j-jb),
                            zone, dAup(jj, j   ), NB, 
                                  dAup(jj, j+jb), NB);

                    /* send off-diagonal block to CPU */
            //cudaMemcpy2DAsync(  A  (j,  j+jb), lda*sizeof(cuDoubleComplex),
                    //                   dAup(jj, j+jb), NB *sizeof(cuDoubleComplex),
                    //                   sizeof(cuDoubleComplex)*jb, n-j-jb,
                    //                   cudaMemcpyDeviceToHost, stream[2]);
                  }

                } /* end for jj */

                /* upload the off-diagonal big panel */
                if( J+JB < n )
            cudaMemcpy2DAsync(  A  (J, J+JB), lda*sizeof(cuDoubleComplex),
                               dAup(0, J+JB), NB *sizeof(cuDoubleComplex),
                               sizeof(cuDoubleComplex)*JB, n-J-JB, 
                               cudaMemcpyDeviceToHost,stream[2]);
          }
        } else {
          /* ========================================================= *
           * Compute the Cholesky factorization A = L*L'.              */

          /* for each big-panel */
          for( J=0; J<n; J+=NB ) {
                JB = min(NB,n-J);

                /* load the new big-panel by block-columns*/
            for (jj=0; jj<JB; jj+=nb) {
                  j  = J+jj;
                  jb = min(nb, (n-j));
              cublasSetMatrix((n-j), jb, sizeof(cuDoubleComplex), 
                              A(j, j), lda, dA(j, jj), ldda);
                }

                /* update with the previous big-panels */
                for( j=0; j<J; j+=nb ) {

                  /* upload the block-column */
              cublasSetMatrix((n-J), nb, sizeof(cuDoubleComplex), 
                              A(J, j), lda, dT(J, 0), ldda);

                  /* update the current big-panel    *
                   * using the previous block-column */
              cublasZherk(MagmaLower, MagmaNoTrans, JB, nb,
                          mdone, dT(J, 0), ldda, 
                          done,  dA(J, 0), ldda);
                  if( J+JB < n )
              cublasZgemm( MagmaNoTrans, MagmaConjTrans, 
                           (n-J-JB), JB, nb,
                           mzone, dT(J+JB, 0), ldda, 
                                  dT(J,    0), ldda,
                           zone,  dA(J+JB, 0), ldda);
                }

                /* for each block-column in the big panel */
            for (jj=0; jj<JB; jj+=nb) {
                  j  = J+jj;
                  jb = min(nb, (n-j));

              /* Update the current diagonal block */
              cublasZherk(MagmaLower, MagmaNoTrans, jb, jj,
                          mdone, dA(j, 0), ldda, 
                          done,  dA(j, jj), ldda);

                  /* upload the current diagonal block to CPU for factorization *
                   * this requires the synchronization before factorization     */
              cudaMemcpy2DAsync(  A(j,j),  lda *sizeof(cuDoubleComplex),
                                 dA(j,jj), ldda*sizeof(cuDoubleComplex),
                                 sizeof(cuDoubleComplex)*jb, jb,
                                 cudaMemcpyDeviceToHost,stream[1]);
                  /* upload the corresponding off-diagonal block-row from previous itrs   *
                   * to CPU. this can wait till end.                                      */
              cudaMemcpy2DAsync(  A(j, J), lda *sizeof(cuDoubleComplex),
                                 dA(j, 0), ldda*sizeof(cuDoubleComplex),
                                 sizeof(cuDoubleComplex)*jb, jj,
                                 cudaMemcpyDeviceToHost,stream[0]);

              if ( (j+jb) < n) {
                        /* update the off-diagonal blocks of the current block-column *
                         * using the previous columns                                 */
                cublasZgemm( MagmaNoTrans, MagmaConjTrans, 
                             (n-j-jb), jb, jj,
                             mzone, dA(j+jb, 0),  ldda, 
                                    dA(j,    0),  ldda,
                             zone,  dA(j+jb, jj), ldda);
              }
                
                  /* CPU wait for the diagonal-block and factor */
              cudaStreamSynchronize(stream[1]);
              lapackf77_zpotrf(MagmaLowerStr, &jb, A(j, j), &lda, info);
              if (*info != 0){
                *info = *info + j;
                break;
              }

              if ( (j+jb) < n) {
                    /* send the diagonal-block to GPU */
                cudaMemcpy2DAsync( dA(j, jj), ldda*sizeof(cuDoubleComplex), 
                                   A(j,   j), lda *sizeof(cuDoubleComplex), 
                                   sizeof(cuDoubleComplex)*jb, jb, 
                                   cudaMemcpyHostToDevice,stream[0]);
                
                        /* GPU do the solves with the current diagonal-block */
                cublasZtrsm( MagmaRight, MagmaLower, MagmaConjTrans, MagmaNonUnit, 
                             (n-j-jb), jb, zone, 
                                                 dA(j,    jj), ldda, 
                             dA(j+jb, jj), ldda);
                  }
                } /* end of for jj */

                /* upload the off-diagonal big panel */
                if( J+JB < n )
            cudaMemcpy2DAsync(  A(J+JB, J), lda *sizeof(cuDoubleComplex),
                               dA(J+JB, 0), ldda*sizeof(cuDoubleComplex),
                               sizeof(cuDoubleComplex)*(n-J-JB), JB,
                               cudaMemcpyDeviceToHost,stream[2]);

          } /* end of for J */
    } /* if upper */
        } /* if nb */
    cudaStreamDestroy(stream[0]);
    cudaStreamDestroy(stream[1]);
    cudaStreamDestroy(stream[2]);

    cublasFree(dt);
    
    return MAGMA_SUCCESS;
} /* magma_zpotrf_ooc */

