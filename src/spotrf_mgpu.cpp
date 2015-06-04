/*
    -- MAGMA (version 1.2.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       June 2012

       @generated s Thu Jun 28 12:30:33 2012

*/
#include "common_magma.h"

/* === Define what BLAS to use ============================================ */
#define PRECISION_s
#if (defined(PRECISION_s) || defined(PRECISION_d)) 
  #define magma_sgemm magmablas_sgemm
  #define magma_strsm magmablas_strsm
#endif

#if (GPUSHMEM >= 200)
  #if (defined(PRECISION_s))
     #undef  magma_sgemm
     #define magma_sgemm magmablas_sgemm_fermi80
  #endif
#endif
/* === End defining what BLAS to use ======================================= */

#define dlA(id, i, j)  (d_lA[id] + (j)*ldda + (i))
#define dlP(id, i, j)  (d_lP[id] + (j)*ldda + (i))

#define dlAT(id, i, j)  (d_lA[id] + (j)*lddat + (i))
#define dlPT(id, i, j)  (d_lP[id] + (j)*nb    + (i))

extern "C" magma_int_t
magma_spotrf_mgpu(int num_gpus, char uplo, magma_int_t n, 
                 float **d_lA, magma_int_t ldda, magma_int_t *info)
{
/*  -- MAGMA (version 1.2.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       June 2012

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


    magma_int_t     j, jb, nb, nb0, nb2, d, id, j_local, j_local2;
    char            uplo_[2] = {uplo, 0};
    float c_one     = MAGMA_S_ONE;
    float c_neg_one = MAGMA_S_NEG_ONE;
    float *work;
    float          d_one     =  1.0;
    float          d_neg_one = -1.0;
    int upper = lapackf77_lsame(uplo_, "U");
    float *d_lP[MagmaMaxGPUs], *dlpanel;
    magma_int_t n_local[MagmaMaxGPUs], lddat_local[MagmaMaxGPUs], lddat, ldpanel;
    cudaStream_t stream[MagmaMaxGPUs][4];

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
        return *info;
    }
    nb = magma_get_spotrf_nb(n);

    if (MAGMA_SUCCESS != magma_smalloc_pinned( &work, n*nb )) {
        *info = MAGMA_ERR_HOST_ALLOC;
        return *info;
    }

    if ((nb <= 1) || (nb >= n)) {
      /*  Use unblocked code. */
          magma_setdevice(0);
      magma_sgetmatrix( n, n, dlA(0,0,0), ldda, work, n );
      lapackf77_spotrf(uplo_, &n, work, &n, info);
      magma_ssetmatrix( n, n, work, n, dlA(0,0,0), ldda );
    } else {

      for( d=0; d<num_gpus; d++ ) {
        magma_setdevice(d);

        /* local-n and local-ld */
        n_local[d] = ((n/nb)/num_gpus)*nb;
        if (d < (n/nb)%num_gpus)
          n_local[d] += nb;
        else if (d == (n/nb)%num_gpus)
          n_local[d] += n%nb;
        lddat_local[d] = ((n_local[d]+31)/32)*32;
        
        if (MAGMA_SUCCESS != magma_smalloc( &d_lP[d], nb*ldda )) {
          for( j=0; j<d; j++ ) {
            magma_setdevice(j);
            magma_free( d_lP[d] );
          }
          *info = MAGMA_ERR_DEVICE_ALLOC;
          return *info;
        }
        magma_queue_create( &stream[d][0] );
        magma_queue_create( &stream[d][1] );
        magma_queue_create( &stream[d][2] );
        magma_queue_create( &stream[d][3] );
      }

      /* Use blocked code. */
      if (upper) 
        {     
          /* Compute the Cholesky factorization A = U'*U. */
          for (j=0; j<n; j+=nb) {

            /* Set the GPU number that holds the current panel */
            id = (j/nb)%num_gpus;
            magma_setdevice(id);

            /* Set the local index where the current panel is */
            j_local = j/(nb*num_gpus);
            jb = min(nb, (n-j));
              
            if( j>0 && (j+jb)<n) {
              /* wait for the off-diagonal column off the current diagonal *
               * and send it to gpus                                       */
              magma_queue_sync( stream[id][2] );
              for( d=0; d<num_gpus; d++ ) {
                magma_setdevice(d);
                if( d != id )
                  magma_ssetmatrix_async( j, jb,
                                          &work[jb],   n,
                                          dlP(d,jb,0), ldda, stream[d][3] );
              }
            }
            
            /* Update the current diagonal block */
            magma_setdevice(id);
            magma_ssyrk(MagmaUpper, MagmaTrans, jb, j, 
                        d_neg_one, dlA(id, 0, nb*j_local), ldda, 
                        d_one,     dlA(id, j, nb*j_local), ldda);

            /* send the diagonal to cpu */
            magma_sgetmatrix_async( jb, jb,
                                    dlA(id, j, nb*j_local), ldda,
                                    work,                   n, stream[id][0] );

            if ( j>0 && (j+jb)<n) {
              /* Compute the local block column of the panel. */
              for( d=0; d<num_gpus; d++ ) {
                magma_setdevice(d);
                j_local2 = j_local+1;
                if( d > id ) j_local2 --;

                /* wait for the off-diagonal */
                if( d != id ) {
                  //magma_queue_sync( stream[id][3] );
                  dlpanel = dlP(d, jb, 0);
                } else {
                  dlpanel = dlA(d, 0, nb*j_local);
                }
                
                /* update the panel */
                magma_sgemm(MagmaTrans, MagmaNoTrans, 
                            jb, (n_local[d]-nb*(j_local2-1)-jb), j, 
                            c_neg_one, dlpanel,                ldda, 
                                       dlA(d, 0, nb*j_local2), ldda,
                            c_one,     dlA(d, j, nb*j_local2), ldda);
              }
            }
          
            /* factor the diagonal */
            magma_setdevice(id);
            magma_queue_sync( stream[id][0] );
            lapackf77_spotrf(MagmaUpperStr, &jb, work, &n, info);
            if (*info != 0) {
              *info = *info + j;
              break;
            }

            /* send the diagonal to gpus */
            if ( (j+jb) < n) {
              for( d=0; d<num_gpus; d++ ) {
                magma_setdevice(d);
                if( d == id ) dlpanel = dlA(d, j, nb*j_local);
                else          dlpanel = dlP(d, 0, 0);
                
                magma_ssetmatrix_async( jb, jb,
                                        work,    n,
                                        dlpanel, ldda, stream[d][1] );
              }
            } else {
              magma_setdevice(id);
              magma_ssetmatrix_async( jb, jb,
                                      work,                   n,
                                      dlA(id, j, nb*j_local), ldda, stream[id][1] );
            }

            /* panel-factorize the off-diagonal */
            if ( (j+jb) < n) {
              for( d=0; d<num_gpus; d++ ) {
                magma_setdevice(d);
                
                /* next column */
                j_local2 = j_local+1;
                if( d > id ) j_local2--;
                if( d == id ) dlpanel = dlA(d, j, nb*j_local);
                else          dlpanel = dlP(d, 0, 0);
                nb0 = min(nb, n_local[d]-nb*j_local2 );
                
                magma_queue_sync( stream[d][1] );
                if( d == (j/nb+1)%num_gpus ) {
                  /* owns the next column, look-ahead the column */
                  /*
                  magma_strsm( MagmaLeft, MagmaUpper, MagmaTrans, MagmaNonUnit, 
                               jb, nb0, c_one,
                               dlpanel,                ldda, 
                               dlA(d, j, nb*j_local2), ldda);
                  */
                  nb2 = n_local[d] - j_local2*nb;
                  magma_strsm( MagmaLeft, MagmaUpper, MagmaTrans, MagmaNonUnit,
                               jb, nb2, c_one,
                               dlpanel,                    ldda,
                               dlA(d, j, nb*j_local2), ldda);

                  
                  /* send the column to cpu */
                  if( j+jb+nb0 < n ) {
                    magma_sgetmatrix_async( (j+jb), nb0,
                                            dlA(d, 0, nb*j_local2), ldda,
                                            &work[nb0],             n, stream[d][2] );
                  }

                  /* update the remaining blocks */
                  /*
                  nb2 = n_local[d] - j_local2*nb - nb0;
                  magma_strsm( MagmaLeft, MagmaUpper, MagmaTrans, MagmaNonUnit, 
                               jb, nb2, c_one,
                               dlpanel,                    ldda, 
                               dlA(d, j, nb*j_local2+nb0), ldda);
                  */
                } else {
                  /* update the entire trailing matrix */
                  nb2 = n_local[d] - j_local2*nb;
                  magma_strsm( MagmaLeft, MagmaUpper, MagmaTrans, MagmaNonUnit, 
                               jb, nb2, c_one,
                               dlpanel,                ldda, 
                               dlA(d, j, nb*j_local2), ldda);
                }
              }
            } /* end of strsm */
          } /* end of for j=1, .., n */
        } else { 
            /* Compute the Cholesky factorization A = L*L'. */

            for (j=0; j<n; j+=nb) {

              /* Set the GPU number that holds the current panel */
              id = (j/nb)%num_gpus;
              magma_setdevice(id);

              /* Set the local index where the current panel is */
              j_local = j/(nb*num_gpus);
          jb = min(nb, (n-j));

                  if( j>0 && (j+jb)<n) {
                    /* wait for the off-diagonal row off the current diagonal *
                         * and send it to gpus                                    */
                magma_queue_sync( stream[id][2] );
            for( d=0; d<num_gpus; d++ ) {
                  magma_setdevice(d);
                      lddat = lddat_local[d];
                          if( d != id ) 
                  magma_ssetmatrix_async( jb, j,
                                          &work[jb*jb], jb,
                                          dlPT(d,0,jb), nb, stream[d][3] );
                        }
                  }

              /* Update the current diagonal block */
              magma_setdevice(id);
                  lddat = lddat_local[id];
                  magma_ssyrk(MagmaLower, MagmaNoTrans, jb, j,
                              d_neg_one, dlAT(id, nb*j_local, 0), lddat,
                              d_one,     dlAT(id, nb*j_local, j), lddat);

                  /* send the diagonal to cpu */
              magma_sgetmatrix_async( jb, jb,
                                      dlAT(id, nb*j_local, j), lddat,
                                      work,                    jb, stream[id][0] );

                  if ( j > 0 && (j+jb) < n) {
                        /* compute the block-rows of the panel */
                        for( d=0; d<num_gpus; d++ ) {
                  magma_setdevice(d);
                          lddat = lddat_local[d];
                          j_local2 = j_local+1;
                          if( d > id ) j_local2 --;

                          /* wait for the off-diagonal */
                  if( d != id ) {
                                  magma_queue_sync( stream[id][3] );
                                  dlpanel = dlPT(d, 0, jb);
                                  ldpanel = nb;
                          } else {
                                  dlpanel = dlAT(d, nb*j_local, 0);
                                  ldpanel = lddat;
                          }

                          /* update the panel */
                      magma_sgemm( MagmaNoTrans, MagmaTrans,
                                   n_local[d]-nb*j_local2, jb, j,
                                   c_neg_one, dlAT(d, nb*j_local2, 0), lddat,
                                              dlpanel,                 ldpanel,
                                   c_one,     dlAT(d, nb*j_local2, j), lddat);
                        }
                  }

                  /* factor the diagonal */
              magma_setdevice(id);
                  magma_queue_sync( stream[id][0] );
                  lapackf77_spotrf(MagmaLowerStr, &jb, work, &jb, info);
                  if (*info != 0) {
                     *info = *info + j;
                     break;
                  }

                  /* send the diagonal to gpus */
              if ( (j+jb) < n) {
                    for( d=0; d<num_gpus; d++ ) {
                  magma_setdevice(d);
                          lddat = lddat_local[d];
                          if( d == id ) {
                                  dlpanel = dlAT(d, nb*j_local, j);
                                  ldpanel = lddat;
                          } else {
                                  dlpanel = dlPT(d, 0, 0);
                                  ldpanel = nb;
                          }
                      magma_ssetmatrix_async( jb, jb,
                                              work,    jb,
                                              dlpanel, ldpanel, stream[d][1] );
                        }
                  } else {
                    magma_ssetmatrix_async( jb, jb,
                                            work,                    jb,
                                            dlAT(id, nb*j_local, j), lddat, stream[id][1] );
                  }
              if ( (j+jb) < n) {
                    for( d=0; d<num_gpus; d++ ) {
                  magma_setdevice(d);
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

                  magma_queue_sync( stream[d][1] );
                          if( d == (j/nb+1)%num_gpus ) {
                            /* owns the next column, look-ahead the column */
                    magma_strsm( MagmaRight, MagmaLower, MagmaTrans, MagmaNonUnit, 
                                 nb0, jb, c_one,
                                 dlpanel,                 ldpanel, 
                                 dlAT(d, nb*j_local2, j), lddat);

                            /* send the column to cpu */
                                if( j+jb+nb0 < n ) {
                      magma_sgetmatrix_async( nb0, j+jb,
                                              dlAT(d, nb*j_local2, 0), lddat,
                                              &work[nb0*nb0],          nb0, stream[d][2] );
                                }

                                /* update the remaining blocks */
                        nb2 = n_local[d] - j_local2*nb - nb0;
                    magma_strsm( MagmaRight, MagmaLower, MagmaTrans, MagmaNonUnit, 
                                 nb2, jb, c_one,
                                 dlpanel,                     ldpanel, 
                                 dlAT(d, nb*j_local2+nb0, j), lddat);
                          } else {
                        /* update the entire trailing matrix */
                        nb2 = n_local[d] - j_local2*nb;
                    magma_strsm( MagmaRight, MagmaLower, MagmaTrans, MagmaNonUnit, 
                                 nb2, jb, c_one,
                                 dlpanel,                 ldpanel, 
                                 dlAT(d, nb*j_local2, j), lddat);
                          }
                        }
                  }
            }
          } /* end of else not upper */

          /* clean up */
          for( d=0; d<num_gpus; d++ ) {
            magma_free( d_lP[d] );
            magma_queue_destroy( stream[d][0] );
            magma_queue_destroy( stream[d][1] );
            magma_queue_destroy( stream[d][2] );
            magma_queue_destroy( stream[d][3] );
          }

    } /* end of not lapack */

        /* free workspace */
        magma_free_pinned( work );

        return *info;
} /* magma_spotrf_mgpu */
