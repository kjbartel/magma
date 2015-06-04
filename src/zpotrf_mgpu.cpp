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
  #define cublasZgemm magmablas_zgemm
  #define cublasZtrsm magmablas_ztrsm
#endif

#if (GPUSHMEM >= 200)
  #if (defined(PRECISION_s))
     #undef  cublasSgemm
     #define cublasSgemm magmablas_sgemm_fermi80
  #endif
#endif
/* === End defining what BLAS to use ======================================= */

#define dlA(id, i, j)  (d_lA[id] + (j)*ldda + (i))
#define dlP(id, i, j)  (d_lP[id] + (j)*ldda + (i))

#define dlAT(id, i, j)  (d_lA[id] + (j)*lddat + (i))
#define dlPT(id, i, j)  (d_lP[id] + (j)*nb    + (i))

extern "C" magma_int_t
magma_zpotrf_mgpu(int num_gpus, char uplo, magma_int_t n, 
                 cuDoubleComplex **d_lA, magma_int_t ldda, magma_int_t *info)
{
/*  -- MAGMA (version 1.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2011

    Purpose   
    =======   
    ZPOTRF computes the Cholesky factorization of a complex Hermitian   
    positive definite matrix dA.   

    The factorization has the form   
       dA = U\*\*H * U,  if UPLO = 'U', or   
       dA = L  * L\*\*H,  if UPLO = 'L',   
    where U is an upper triangular matrix and L is lower triangular.   

    This is the block version of the algorithm, calling Level 3 BLAS.   

    Arguments   
    =========   
    UPLO    (input) CHARACTER*1   
            = 'U':  Upper triangle of dA is stored;   
            = 'L':  Lower triangle of dA is stored.   

    N       (input) INTEGER   
            The order of the matrix dA.  N >= 0.   

    dA      (input/output) COMPLEX_16 array on the GPU, dimension (LDDA,N)   
            On entry, the Hermitian matrix dA.  If UPLO = 'U', the leading   
            N-by-N upper triangular part of dA contains the upper   
            triangular part of the matrix dA, and the strictly lower   
            triangular part of dA is not referenced.  If UPLO = 'L', the   
            leading N-by-N lower triangular part of dA contains the lower   
            triangular part of the matrix dA, and the strictly upper   
            triangular part of dA is not referenced.   

            On exit, if INFO = 0, the factor U or L from the Cholesky   
            factorization dA = U\*\*H*U or dA = L*L\*\*H.   

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


    magma_int_t     j, jb, nb, nb0, nb2, d, id, j_local, j_local2;
    char            uplo_[2] = {uplo, 0};
    cuDoubleComplex zone  = MAGMA_Z_ONE;
    cuDoubleComplex mzone = MAGMA_Z_NEG_ONE;
    cuDoubleComplex *work;
    double          done  = (double) 1.0;
    double          mdone = (double)-1.0;
    long int        upper = lapackf77_lsame(uplo_, "U");
    cuDoubleComplex *d_lP[4], *dlpanel;
        magma_int_t n_local[4], lddat_local[4], lddat, ldpanel;
    static cudaStream_t stream[4][4];

    *info = 0;
    if ( (! upper) && (! lapackf77_lsame(uplo_, "L")) ) {
        *info = -1;
    } else if (n < 0) {
        *info = -2;
    } else if (ldda < max(1,n)) {
        *info = -4;
    }
    if (*info != 0) {
                magma_xerbla( __func__, -(*info) );
        return MAGMA_ERR_ILLEGAL_VALUE;
        }
    nb = magma_get_zpotrf_nb(n);

    if (cudaSuccess != cudaMallocHost( (void**)&work, n*nb*sizeof(cuDoubleComplex) ) ) {
          *info = -6;
          return MAGMA_ERR_HOSTALLOC;
    }

    if ((nb <= 1) || (nb >= n)) {
      /*  Use unblocked code. */
          cudaSetDevice(0);
      cublasGetMatrix(n, n, sizeof(cuDoubleComplex), dlA(0,0,0), ldda, work, n);
      lapackf77_zpotrf(uplo_, &n, work, &n, info);
      cublasSetMatrix(n, n, sizeof(cuDoubleComplex), work, n, dlA(0,0,0), ldda);
    } else {

      for( d=0; d<num_gpus; d++ ) {
        cudaSetDevice(d);

        /* local-n and local-ld */
        n_local[d] = ((n/nb)/num_gpus)*nb;
        if (d < (n/nb)%num_gpus)
          n_local[d] += nb;
        else if (d == (n/nb)%num_gpus)
          n_local[d] += n%nb;
        lddat_local[d] = ((n_local[d]+31)/32)*32;
        
        if ( CUBLAS_STATUS_SUCCESS != cublasAlloc(nb*ldda, sizeof(cuDoubleComplex), (void**)&d_lP[d]) ) {
          for( j=0; j<d; j++ ) {
            cudaSetDevice(j);
            cublasFree(d_lP[d]);
          }
          return MAGMA_ERR_CUBLASALLOC;
        }
        cudaStreamCreate(&stream[d][0]);
        cudaStreamCreate(&stream[d][1]);
        cudaStreamCreate(&stream[d][2]);
        cudaStreamCreate(&stream[d][3]);
      }

      /* Use blocked code. */
      if (upper) 
        {     
          /* Compute the Cholesky factorization A = U'*U. */
          for (j=0; j<n; j+=nb) {

            /* Set the GPU number that holds the current panel */
            id = (j/nb)%num_gpus;
            cudaSetDevice(id);

            /* Set the local index where the current panel is */
            j_local = j/(nb*num_gpus);
            jb = min(nb, (n-j));
              
            if( j>0 && (j+jb)<n) {
              /* wait for the off-diagonal column off the current diagonal *
               * and send it to gpus                                       */
              cudaStreamSynchronize(stream[id][2]);
              for( d=0; d<num_gpus; d++ ) {
                cudaSetDevice(d);
                if( d != id )
                  cudaMemcpy2DAsync(dlP(d,jb,0), ldda*sizeof(cuDoubleComplex), 
                                    &work[jb],   n   *sizeof(cuDoubleComplex), 
                                    j*sizeof(cuDoubleComplex), jb, 
                                    cudaMemcpyHostToDevice,stream[d][3]);
              }
            }
            
            /* Update the current diagonal block */
            cudaSetDevice(id);
            cublasZherk(MagmaUpper, MagmaConjTrans, jb, j, 
                        mdone, dlA(id, 0, nb*j_local), ldda, 
                        done,  dlA(id, j, nb*j_local), ldda);

            /* send the diagonal to cpu */
            cudaMemcpy2DAsync(work, n *sizeof(cuDoubleComplex), 
                              dlA(id, j, nb*j_local), ldda*sizeof(cuDoubleComplex), 
                              jb*sizeof(cuDoubleComplex), jb, 
                              cudaMemcpyDeviceToHost,stream[id][0]);

            if ( j>0 && (j+jb)<n) {
              /* Compute the local block column of the panel. */
              for( d=0; d<num_gpus; d++ ) {
                cudaSetDevice(d);
                j_local2 = j_local+1;
                if( d > id ) j_local2 --;

                /* wait for the off-diagonal */
                if( d != id ) {
                  //cudaStreamSynchronize(stream[id][3]);
                  dlpanel = dlP(d, jb, 0);
                } else {
                  dlpanel = dlA(d, 0, nb*j_local);
                }
                
                /* update the panel */
                cublasZgemm(MagmaConjTrans, MagmaNoTrans, 
                            jb, (n_local[d]-nb*(j_local2-1)-jb), j, 
                            mzone, dlpanel,                ldda, 
                            dlA(d, 0, nb*j_local2), ldda,
                            zone,  dlA(d, j, nb*j_local2), ldda);
              }
            }
          
            /* factor the diagonal */
            cudaSetDevice(id);
            cudaStreamSynchronize(stream[id][0]);
            lapackf77_zpotrf(MagmaUpperStr, &jb, work, &n, info);
            if (*info != 0) {
              *info = *info + j;
              break;
            }

            /* send the diagonal to gpus */
            if ( (j+jb) < n) {
              for( d=0; d<num_gpus; d++ ) {
                cudaSetDevice(d);
                if( d == id ) dlpanel = dlA(d, j, nb*j_local);
                else          dlpanel = dlP(d, 0, 0);
                
                cudaMemcpy2DAsync( dlpanel, ldda*sizeof(cuDoubleComplex), 
                                   work,    n   *sizeof(cuDoubleComplex), 
                                   sizeof(cuDoubleComplex)*jb, jb, 
                                   cudaMemcpyHostToDevice,stream[d][1]);
              }
            } else {
              cudaSetDevice(id);
              cudaMemcpy2DAsync( dlA(id, j, nb*j_local), ldda*sizeof(cuDoubleComplex), 
                                 work,                   n   *sizeof(cuDoubleComplex), 
                                 sizeof(cuDoubleComplex)*jb, jb, 
                                 cudaMemcpyHostToDevice,stream[id][1]);
            }

            /* panel-factorize the off-diagonal */
            if ( (j+jb) < n) {
              for( d=0; d<num_gpus; d++ ) {
                cudaSetDevice(d);
                
                /* next column */
                j_local2 = j_local+1;
                if( d > id ) j_local2--;
                if( d == id ) dlpanel = dlA(d, j, nb*j_local);
                else          dlpanel = dlP(d, 0, 0);
                nb0 = min(nb, n_local[d]-nb*j_local2 );
                
                cudaStreamSynchronize(stream[d][1]);
                if( d == (j/nb+1)%num_gpus ) {
                  /* owns the next column, look-ahead the column */
                  /*
                  cublasZtrsm( MagmaLeft, MagmaUpper, MagmaConjTrans, MagmaNonUnit, 
                               jb, nb0, zone,
                               dlpanel,                ldda, 
                               dlA(d, j, nb*j_local2), ldda);
                  */
                  nb2 = n_local[d] - j_local2*nb;
                  cublasZtrsm( MagmaLeft, MagmaUpper, MagmaConjTrans, MagmaNonUnit,
                               jb, nb2, zone,
                               dlpanel,                    ldda,
                               dlA(d, j, nb*j_local2), ldda);

                  
                  /* send the column to cpu */
                  if( j+jb+nb0 < n ) {
                    cudaMemcpy2DAsync(&work[nb0],             n   *sizeof(cuDoubleComplex), 
                                      dlA(d, 0, nb*j_local2), ldda*sizeof(cuDoubleComplex), 
                                      (j+jb)*sizeof(cuDoubleComplex), nb0, 
                                      cudaMemcpyDeviceToHost,stream[d][2]);
                  }

                  /* update the remaining blocks */
                  /*
                  nb2 = n_local[d] - j_local2*nb - nb0;
                  cublasZtrsm( MagmaLeft, MagmaUpper, MagmaConjTrans, MagmaNonUnit, 
                               jb, nb2, zone,
                               dlpanel,                    ldda, 
                               dlA(d, j, nb*j_local2+nb0), ldda);
                  */
                } else {
                  /* update the entire trailing matrix */
                  nb2 = n_local[d] - j_local2*nb;
                  cublasZtrsm( MagmaLeft, MagmaUpper, MagmaConjTrans, MagmaNonUnit, 
                               jb, nb2, zone,
                               dlpanel,                ldda, 
                               dlA(d, j, nb*j_local2), ldda);
                }
              }
            } /* end of ztrsm */
          } /* end of for j=1, .., n */
        } else { 
            /* Compute the Cholesky factorization A = L*L'. */

            for (j=0; j<n; j+=nb) {

              /* Set the GPU number that holds the current panel */
              id = (j/nb)%num_gpus;
              cudaSetDevice(id);

              /* Set the local index where the current panel is */
              j_local = j/(nb*num_gpus);
          jb = min(nb, (n-j));

                  if( j>0 && (j+jb)<n) {
                    /* wait for the off-diagonal row off the current diagonal *
                         * and send it to gpus                                    */
                cudaStreamSynchronize(stream[id][2]);
            for( d=0; d<num_gpus; d++ ) {
                  cudaSetDevice(d);
                      lddat = lddat_local[d];
                          if( d != id ) 
                  cudaMemcpy2DAsync(dlPT(d,0,jb), nb*sizeof(cuDoubleComplex), 
                                                    &work[jb*jb], jb*sizeof(cuDoubleComplex), 
                                    jb*sizeof(cuDoubleComplex), j, 
                                    cudaMemcpyHostToDevice,stream[d][3]);
                        }
                  }

              /* Update the current diagonal block */
              cudaSetDevice(id);
                  lddat = lddat_local[id];
                  cublasZherk(MagmaLower, MagmaNoTrans, jb, j,
                              mdone, dlAT(id, nb*j_local, 0), lddat,
                              done,  dlAT(id, nb*j_local, j), lddat);

                  /* send the diagonal to cpu */
              cudaMemcpy2DAsync(work,                    jb   *sizeof(cuDoubleComplex), 
                                dlAT(id, nb*j_local, j), lddat*sizeof(cuDoubleComplex), 
                                jb*sizeof(cuDoubleComplex), jb, 
                                cudaMemcpyDeviceToHost,stream[id][0]);

                  if ( j > 0 && (j+jb) < n) {
                        /* compute the block-rows of the panel */
                        for( d=0; d<num_gpus; d++ ) {
                  cudaSetDevice(d);
                          lddat = lddat_local[d];
                          j_local2 = j_local+1;
                          if( d > id ) j_local2 --;

                          /* wait for the off-diagonal */
                  if( d != id ) {
                                  cudaStreamSynchronize(stream[id][3]);
                                  dlpanel = dlPT(d, 0, jb);
                                  ldpanel = nb;
                          } else {
                                  dlpanel = dlAT(d, nb*j_local, 0);
                                  ldpanel = lddat;
                          }

                          /* update the panel */
                      cublasZgemm( MagmaNoTrans, MagmaConjTrans,
                                   n_local[d]-nb*j_local2, jb, j,
                                   mzone, dlAT(d, nb*j_local2, 0), lddat,
                                          dlpanel,                 ldpanel,
                                   zone,  dlAT(d, nb*j_local2, j), lddat);
                        }
                  }

                  /* factor the diagonal */
              cudaSetDevice(id);
                  cudaStreamSynchronize(stream[id][0]);
                  lapackf77_zpotrf(MagmaLowerStr, &jb, work, &jb, info);
                  if (*info != 0) {
                     *info = *info + j;
                     break;
                  }

                  /* send the diagonal to gpus */
              if ( (j+jb) < n) {
                    for( d=0; d<num_gpus; d++ ) {
                  cudaSetDevice(d);
                          lddat = lddat_local[d];
                          if( d == id ) {
                                  dlpanel = dlAT(d, nb*j_local, j);
                                  ldpanel = lddat;
                          } else {
                                  dlpanel = dlPT(d, 0, 0);
                                  ldpanel = nb;
                          }
                      cudaMemcpy2DAsync(dlpanel, ldpanel*sizeof(cuDoubleComplex),
                                        work,    jb     *sizeof(cuDoubleComplex),
                                        sizeof(cuDoubleComplex)*jb, jb,
                                        cudaMemcpyHostToDevice,stream[d][1]);
                        }
                  } else {
                    cudaMemcpy2DAsync(dlAT(id, nb*j_local, j), lddat*sizeof(cuDoubleComplex),
                                      work,                    jb   *sizeof(cuDoubleComplex),
                                      sizeof(cuDoubleComplex)*jb, jb,
                                      cudaMemcpyHostToDevice,stream[id][1]);
                  }
              if ( (j+jb) < n) {
                    for( d=0; d<num_gpus; d++ ) {
                  cudaSetDevice(d);
                          lddat = lddat_local[d];

                          /* next column */
                          j_local2 = j_local+1;
                          if( d > id ) j_local2--;
                          if( d == id ) {
                                  dlpanel = dlAT(d, nb*j_local, j);
                                  ldpanel = lddat;
                          } else {         
                                  dlpanel = dlPT(d, 0, 0);
                                  ldpanel = nb;
                          }
                          nb0 = min(nb, n_local[d]-nb*j_local2 );
                          nb0 = min(nb, lddat-nb*j_local2 );

                  cudaStreamSynchronize(stream[d][1]);
                          if( d == (j/nb+1)%num_gpus ) {
                            /* owns the next column, look-ahead the column */
                    cublasZtrsm( MagmaRight, MagmaLower, MagmaConjTrans, MagmaNonUnit, 
                                 nb0, jb, zone,
                                 dlpanel,                 ldpanel, 
                                 dlAT(d, nb*j_local2, j), lddat);

                            /* send the column to cpu */
                                if( j+jb+nb0 < n ) {
                      cudaMemcpy2DAsync(&work[nb0*nb0],          nb0  *sizeof(cuDoubleComplex), 
                                        dlAT(d, nb*j_local2, 0), lddat*sizeof(cuDoubleComplex), 
                                        nb0*sizeof(cuDoubleComplex), j+jb,
                                        cudaMemcpyDeviceToHost,stream[d][2]);
                                }

                                /* update the remaining blocks */
                        nb2 = n_local[d] - j_local2*nb - nb0;
                    cublasZtrsm( MagmaRight, MagmaLower, MagmaConjTrans, MagmaNonUnit, 
                                 nb2, jb, zone,
                                 dlpanel,                     ldpanel, 
                                 dlAT(d, nb*j_local2+nb0, j), lddat);
                          } else {
                        /* update the entire trailing matrix */
                        nb2 = n_local[d] - j_local2*nb;
                    cublasZtrsm( MagmaRight, MagmaLower, MagmaConjTrans, MagmaNonUnit, 
                                 nb2, jb, zone,
                                 dlpanel,                 ldpanel, 
                                 dlAT(d, nb*j_local2, j), lddat);
                          }
                        }
                  }
            }
          } /* end of else not upper */

          /* clean up */
          for( d=0; d<num_gpus; d++ ) {
            cublasFree( d_lP[d] );
            cudaStreamDestroy(stream[d][0]);
            cudaStreamDestroy(stream[d][1]);
            cudaStreamDestroy(stream[d][2]);
            cudaStreamDestroy(stream[d][3]);
          }

    } /* end of not lapack */

        /* free workspace */
        cudaFreeHost(work);

        return MAGMA_SUCCESS;
} /* magma_zpotrf_mgpu */
