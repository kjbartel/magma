/*
    -- MAGMA (version 1.2.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       May 2012

       @precisions normal z -> s d c

*/
#include "common_magma.h"

/* === Define what BLAS to use ============================================ */
#define PRECISION_z
#if (defined(PRECISION_s) || defined(PRECISION_d)) 
  #define magma_zgemm magmablas_zgemm
  #define magma_ztrsm magmablas_ztrsm
#endif

#if (GPUSHMEM >= 200)
  #if (defined(PRECISION_s))
     #undef  magma_sgemm
     #define magma_sgemm magmablas_sgemm_fermi80
  #endif
#endif
/* === End defining what BLAS to use ======================================= */
#define A(i, j)  (a   +((j)+off_j)*lda  + (i)+off_i)

#define dlA(id, i, j)  (d_lA[(id)] + (j)*ldda + (i))
#define dlP(id, i, j)  (d_lP[(id)] + (j)*ldda + (i))

#define dlAT(id, i, j)  (d_lA[(id)] + (j)*ldda + (i))
#define dlPT(id, i, j)  (d_lP[(id)] + (j)*nb   + (i))

#define VERSION2
extern "C" magma_int_t
magma_zpotrf2_mgpu(int num_gpus, char uplo, magma_int_t m, magma_int_t n, magma_int_t off_i, magma_int_t off_j, magma_int_t nb,
                   cuDoubleComplex **d_lA, magma_int_t ldda, cuDoubleComplex **d_lP, magma_int_t lddp, 
                   cuDoubleComplex *a, magma_int_t lda, cudaStream_t stream[][4], magma_int_t *info ) 
{
/*  -- MAGMA (version 1.2.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       May 2012

    Purpose   
    =======   
    ZPOTRF computes the Cholesky factorization of a complex Hermitian   
    positive definite matrix dA.   
    Auxiliary subroutine for zpotrf2_ooc. It is multiple gpu interface to compute 
    Cholesky of a "rectangular" matrix.

    The factorization has the form   
       dA = U**H * U,  if UPLO = 'U', or   
       dA = L  * L**H,  if UPLO = 'L',   
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
            factorization dA = U**H * U or dA = L * L**H.   

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


    magma_int_t     j, jb, nb0, nb2, d, id, j_local, j_local2;
    char            uplo_[2] = {uplo, 0};
    cuDoubleComplex c_one     = MAGMA_Z_ONE;
    cuDoubleComplex c_neg_one = MAGMA_Z_NEG_ONE;
    double          d_one     =  1.0;
    double          d_neg_one = -1.0;
    long int        upper = lapackf77_lsame(uplo_, "U");
    cuDoubleComplex *dlpanel;
    magma_int_t n_local[4], ldpanel;

    *info = 0;
    if ( (! upper) && (! lapackf77_lsame(uplo_, "L")) ) {
        *info = -1;
    } else if (n < 0) {
        *info = -2;
    } else if (!upper && num_gpus*ldda < max(1,n)) {
        *info = -4;
    } else if (upper && ldda < max(1,m)) {
        *info = -4;
    }
    if (*info != 0) {
        magma_xerbla( __func__, -(*info) );
        return *info;
    }
    //nb = magma_get_zpotrf_nb(n);
    {

      /* Use blocked code. */
      for( d=0; d<num_gpus; d++ ) {
        magma_setdevice(d);
        /* local-n and local-ld */
        if (upper) {
          n_local[d] = ((n/nb)/num_gpus)*nb;
          if (d < (n/nb)%num_gpus)
            n_local[d] += nb;
          else if (d == (n/nb)%num_gpus)
            n_local[d] += n%nb;
        } else {
          n_local[d] = ((m/nb)/num_gpus)*nb;
          if (d < (m/nb)%num_gpus)
            n_local[d] += nb;
          else if (d == (m/nb)%num_gpus)
            n_local[d] += m%nb;
        }
      }

      if (upper) 
      {     
      /* Compute the Cholesky factorization A = U'*U. */
      for (j=0; j<m; j+=nb) {

        /* Set the GPU number that holds the current panel */
        id = (j/nb)%num_gpus;
        magma_setdevice(id);

        /* Set the local index where the current panel is */
        j_local = j/(nb*num_gpus);
        jb = min(nb, (m-j));
          
#if defined(VERSION1) || defined(VERSION2)
        if( j>0 && (j+jb)<n && num_gpus > 1 ) {
          /* wait for the off-diagonal column off the current diagonal *
           * and send it to gpus                                       */
          magma_queue_sync( stream[id][0] );
          for( d=0; d<num_gpus; d++ ) {
            magma_setdevice(d);
            if( d != id )
              magma_zsetmatrix_async( j, jb,
                                      A(0,j),      lda,
                                      dlP(d,jb,0), lddp, stream[d][3] );
          }
        }
#endif

        /* Update the current diagonal block */
        magma_setdevice(id);
        magma_zherk(MagmaUpper, MagmaConjTrans, jb, j, 
                    d_neg_one, dlA(id, 0, nb*j_local), ldda, 
                    d_one,     dlA(id, j, nb*j_local), ldda);

#if defined(VERSION1) || defined(VERSION2)
        /* send the diagonal to cpu */
        magma_zgetmatrix_async( jb, jb,
                                dlA(id, j, nb*j_local), ldda,
                                A(j,j),                 lda, stream[id][0] );

        if ( j>0 && (j+jb)<n) {
          /* Compute the local block column of the panel. */
          for( d=0; d<num_gpus; d++ ) {
            magma_setdevice(d);
            j_local2 = j_local+1;
            if( d > id ) j_local2 --;

            /* wait for the off-diagonal */
            if( d != id ) {
              dlpanel = dlP(d, jb, 0);
              ldpanel = lddp;
            } else {
              dlpanel = dlA(d, 0, nb*j_local);
              ldpanel = ldda;
            }
        
            /* update the panel */
            magma_zgemm(MagmaConjTrans, MagmaNoTrans, 
                        jb, (n_local[d]-nb*(j_local2-1)-jb), j, 
                        c_neg_one, dlpanel,                ldpanel, 
                                   dlA(d, 0, nb*j_local2), ldda,
                        c_one,     dlA(d, j, nb*j_local2), ldda);
          }
        }
#elif defined(VERSION3)
        /* send the whole column to cpu */
        magma_zgetmatrix_async( (j+jb), jb,
                                dlA(id, 0, nb*j_local), ldda,
                                A(0, j),                lda, stream[id][0] );
#endif

        /* wait for panel at cpu */
        magma_setdevice(id);
        magma_queue_sync( stream[id][0] );
#ifdef VERSION3
        if ( j>0 && (j+jb)<n) {
          /* send off-diagonals to gpus */
          for( d=0; d<num_gpus; d++ ) {
            magma_setdevice(d);
            if( d != id )
              magma_zsetmatrix_async( j, jb,
                                      A(0,j),      lda,
                                      dlP(d,jb,0), lddp, stream[d][3] );
          }

          /* Compute the local block column of the panel. */
          for( d=0; d<num_gpus; d++ ) {
            magma_setdevice(d);
            j_local2 = j_local+1;
            if( d > id ) j_local2 --;

            /* wait for the off-diagonal */
            if( d != id ) {
              magma_queue_sync( stream[id][3] );
              dlpanel = dlP(d, jb, 0);
              ldpanel = lddp;
            } else {
              dlpanel = dlA(d, 0, nb*j_local);
              ldpanel = ldda;
            }
        
            /* update the panel */
            magma_zgemm(MagmaConjTrans, MagmaNoTrans, 
                        jb, (n_local[d]-nb*(j_local2-1)-jb), j, 
                        c_neg_one, dlpanel,                ldpanel, 
                                   dlA(d, 0, nb*j_local2), ldda,
                        c_one,     dlA(d, j, nb*j_local2), ldda);
          }
        }
#endif
        /* factor the diagonal */
        lapackf77_zpotrf(MagmaUpperStr, &jb, A(j,j), &lda, info);
        if (*info != 0) {
          *info = *info + j;
          break;
        }

        /* send the diagonal to gpus */
        if ( (j+jb) < n) {
          for( d=0; d<num_gpus; d++ ) {
            magma_setdevice(d);
            if( d == id ) {
                dlpanel = dlA(d, j, nb*j_local);
                ldpanel = ldda;
            } else {
                dlpanel = dlP(d, 0, 0);
                ldpanel = lddp;
            }
            magma_zsetmatrix_async( jb, jb,
                                    A(j,j),  lda,
                                    dlpanel, ldpanel, stream[d][1] );
          }
        } else {
          magma_setdevice(id);
          magma_zsetmatrix_async( jb, jb,
                                  A(j,j),                 lda,
                                  dlA(id, j, nb*j_local), ldda, stream[id][1] );
        }

        /* panel-factorize the off-diagonal */
        if ( (j+jb) < n) {
          for( d=0; d<num_gpus; d++ ) {
            magma_setdevice(d);
        
            /* next column */
            j_local2 = j_local+1;
            if( d > id ) j_local2--;
            if( d == id ) {
                dlpanel = dlA(d, j, nb*j_local);
                ldpanel = ldda;
            } else {
                dlpanel = dlP(d, 0, 0);
                ldpanel = lddp;
            }
            nb0 = min(nb, n_local[d]-nb*j_local2 );
        
            magma_queue_sync( stream[d][1] );
            if( d == (j/nb+1)%num_gpus ) {
              /* owns the next column, look-ahead the column */
#ifdef  VERSION1
              magma_ztrsm( MagmaLeft, MagmaUpper, MagmaConjTrans, MagmaNonUnit, 
                           jb, nb0, c_one,
                           dlpanel,                ldpanel, 
                           dlA(d, j, nb*j_local2), ldda);

              /* send the column to cpu */
              if( j+jb < m ) 
              {
                magma_zgetmatrix_async( (j+jb), nb0,
                                        dlA(d, 0, nb*j_local2), ldda,
                                        A(0,j+jb),              lda, stream[d][0] );
              }

              /* update the remaining blocks */
              nb2 = n_local[d] - j_local2*nb - nb0;
              magma_ztrsm( MagmaLeft, MagmaUpper, MagmaConjTrans, MagmaNonUnit, 
                           jb, nb2, c_one,
                           dlpanel,                    ldpanel, 
                           dlA(d, j, nb*j_local2+nb0), ldda);
#elif defined (VERSION2)
              nb2 = n_local[d] - j_local2*nb;
              magma_ztrsm( MagmaLeft, MagmaUpper, MagmaConjTrans, MagmaNonUnit,
                           jb, nb2, c_one,
                           dlpanel,                ldpanel,
                           dlA(d, j, nb*j_local2), ldda);

              /* send the column to cpu */
              if( j+jb < m ) 
              {
                magma_zgetmatrix_async( (j+jb), nb0,
                                        dlA(d, 0, nb*j_local2), ldda,
                                        A(0,j+jb),              lda, stream[d][0] );
              }
#elif defined (VERSION3)
              nb2 = n_local[d] - j_local2*nb;
              magma_ztrsm( MagmaLeft, MagmaUpper, MagmaConjTrans, MagmaNonUnit,
                           jb, nb2, c_one,
                           dlpanel,                ldpanel,
                           dlA(d, j, nb*j_local2), ldda);
#endif
          
            } else {
              /* update the entire trailing matrix */
              nb2 = n_local[d] - j_local2*nb;
              magma_ztrsm( MagmaLeft, MagmaUpper, MagmaConjTrans, MagmaNonUnit, 
                           jb, nb2, c_one,
                           dlpanel,                ldpanel, 
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
          magma_setdevice(id);

          /* Set the local index where the current panel is */
          j_local = j/(nb*num_gpus);
          jb = min(nb, (n-j));

#if defined(VERSION1) || defined(VERSION2)
          if( j>0 ) {
            /* wait for the off-diagonal row off the current diagonal *
             * and send it to gpus                                    */
            magma_queue_sync( stream[id][0] );
            if( (j+jb)<m && num_gpus > 1 ) {
              for( d=0; d<num_gpus; d++ ) {
                magma_setdevice(d);
                if( d != id ) 
                magma_zsetmatrix_async( jb, j,
                                        A(j,0),       lda,
                                        dlPT(d,0,jb), nb, stream[d][3] );
              }
            }
          }
#endif

          /* Update the current diagonal block */
          magma_setdevice(id);
          magma_zherk(MagmaLower, MagmaNoTrans, jb, j,
                      d_neg_one, dlAT(id, nb*j_local, 0), ldda,
                      d_one,     dlAT(id, nb*j_local, j), ldda);

#if defined(VERSION1) || defined(VERSION2)
          /* send the diagonal to cpu */
          magma_zgetmatrix_async( jb, jb,
                                  dlAT(id, nb*j_local, j), ldda,
                                  A(j,j),                  lda, stream[id][0] );

          if ( j > 0 && (j+jb) < m) {
            /* compute the block-rows of the panel */
            for( d=0; d<num_gpus; d++ ) {
              magma_setdevice(d);
              j_local2 = j_local+1;
              if( d > id ) j_local2 --;

              /* wait for the off-diagonal */
              if( d != id ) {
                  magma_queue_sync( stream[id][3] );
                  dlpanel = dlPT(d, 0, jb);
                  ldpanel = nb;
              } else {
                  dlpanel = dlAT(d, nb*j_local, 0);
                  ldpanel = ldda;
              }

              /* update the panel */
              magma_zgemm( MagmaNoTrans, MagmaConjTrans,
                           n_local[d]-nb*j_local2, jb, j,
                           c_neg_one, dlAT(d, nb*j_local2, 0), ldda,
                                      dlpanel,                 ldpanel,
                           c_one,     dlAT(d, nb*j_local2, j), ldda);
            }
          }
#elif defined(VERSION3)
          /* send the whole row to cpu */
          magma_zgetmatrix_async( jb, (j+jb),
                                  dlAT(id, nb*j_local, 0), ldda,
                                  A(j,0),                  lda, stream[id][0] );
#endif

          /* wait for the panel at cpu */
          magma_setdevice(id);
          magma_queue_sync( stream[id][0] );
#ifdef VERSION3
          if( j>0 && (j+jb)<m ) {
            /* send off-diagonals to gpus */
            for( d=0; d<num_gpus; d++ ) {
              magma_setdevice(d);
              if( d != id ) 
              magma_zsetmatrix_async( jb, j,
                                      A(j,0),       lda,
                                      dlPT(d,0,jb), nb, stream[d][3] );
            }

            /* compute the block-rows of the panel */
            for( d=0; d<num_gpus; d++ ) {
              magma_setdevice(d);
              j_local2 = j_local+1;
              if( d > id ) j_local2 --;

              /* wait for the off-diagonal */
              if( d != id ) {
                  magma_queue_sync( stream[id][3] );
                  dlpanel = dlPT(d, 0, jb);
                  ldpanel = nb;
              } else {
                  dlpanel = dlAT(d, nb*j_local, 0);
                  ldpanel = ldda;
              }

              /* update the panel */
              magma_zgemm( MagmaNoTrans, MagmaConjTrans,
                           n_local[d]-nb*j_local2, jb, j,
                           c_neg_one, dlAT(d, nb*j_local2, 0), ldda,
                                      dlpanel,                 ldpanel,
                           c_one,     dlAT(d, nb*j_local2, j), ldda);
            }
          }
#endif
          /* factor the diagonal */
          lapackf77_zpotrf(MagmaLowerStr, &jb, A(j,j), &lda, info);
          if (*info != 0) {
             *info = *info + j;
             break;
          }

          /* send the diagonal to gpus */
          if ( (j+jb) < m) {
            for( d=0; d<num_gpus; d++ ) {
              magma_setdevice(d);
              if( d == id ) {
                  dlpanel = dlAT(d, nb*j_local, j);
                  ldpanel = ldda;
              } else {
                  dlpanel = dlPT(d, 0, 0);
                  ldpanel = nb;
              }
              magma_zsetmatrix_async( jb, jb,
                                      A(j,j),  lda,
                                      dlpanel, ldpanel, stream[d][1] );
            }
          } else {
            magma_setdevice(id);
            magma_zsetmatrix_async( jb, jb,
                                    A(j,j),                  lda,
                                    dlAT(id, nb*j_local, j), ldda, stream[id][1] );
          }
          if ( (j+jb) < m) {
            for( d=0; d<num_gpus; d++ ) {
              magma_setdevice(d);

              /* next column */
              j_local2 = j_local+1;
              if( d > id ) j_local2--;
              if( d == id ) {
                  dlpanel = dlAT(d, nb*j_local, j);
                  ldpanel = ldda;
              } else {         
                  dlpanel = dlPT(d, 0, 0);
                  ldpanel = nb;
              }
              nb0 = min(nb, n_local[d]-nb*j_local2 );
              //nb0 = min(nb, ldda-nb*j_local2 );

              magma_queue_sync( stream[d][1] );
              if( d == (j/nb+1)%num_gpus ) {
#ifdef VERSION1
                /* owns the next column, look-ahead the column */
                magma_ztrsm( MagmaRight, MagmaLower, MagmaConjTrans, MagmaNonUnit, 
                             nb0, jb, c_one,
                             dlpanel,                 ldpanel, 
                             dlAT(d, nb*j_local2, j), ldda);

                /* send the column to cpu */
                //if( j+jb+nb0 < n && num_gpus > 1 ) {
                if( j+jb < n ) {
                  magma_zgetmatrix_async( nb0, j+jb,
                                          dlAT(d, nb*j_local2, 0), ldda,
                                          A(j+jb,0),               lda, stream[d][0] );
                }

                /* update the remaining blocks */
                nb2 = n_local[d] - j_local2*nb - nb0;
                magma_ztrsm( MagmaRight, MagmaLower, MagmaConjTrans, MagmaNonUnit, 
                             nb2, jb, c_one,
                             dlpanel,                     ldpanel, 
                             dlAT(d, nb*j_local2+nb0, j), ldda);
#elif defined (VERSION2)
                nb2 = n_local[d] - j_local2*nb;
                magma_ztrsm( MagmaRight, MagmaLower, MagmaConjTrans, MagmaNonUnit, 
                             nb2, jb, c_one,
                             dlpanel,                 ldpanel, 
                             dlAT(d, nb*j_local2, j), ldda);

                /* send the column to cpu */
                //if( j+jb+nb0 < n && num_gpus > 1 ) {
                if( j+jb < n ) {
                  magma_zgetmatrix_async( nb0, j+jb,
                                          dlAT(d, nb*j_local2, 0), ldda,
                                          A(j+jb,0),               lda, stream[d][0] );
                }
#elif defined (VERSION3)
                nb2 = n_local[d] - j_local2*nb;
                magma_ztrsm( MagmaRight, MagmaLower, MagmaConjTrans, MagmaNonUnit, 
                             nb2, jb, c_one,
                             dlpanel,                 ldpanel, 
                             dlAT(d, nb*j_local2, j), ldda);
#endif
              } else {
                /* update the entire trailing matrix */
                nb2 = n_local[d] - j_local2*nb;
                magma_ztrsm( MagmaRight, MagmaLower, MagmaConjTrans, MagmaNonUnit, 
                             nb2, jb, c_one,
                             dlpanel,                 ldpanel, 
                             dlAT(d, nb*j_local2, j), ldda);
              }
            }
          }
        }
      } /* end of else not upper */

    } /* end of not lapack */

    return *info;
} /* magma_zpotrf_mgpu */

#undef A
#define A(i, j)  (a +(j)*lda  + (i))
#define dA(d, i, j) (dwork[(d)]+(j)*ldda + (i))

extern "C" magma_int_t
magma_zhtodpo(int num_gpus, char *uplo, magma_int_t m, magma_int_t n, magma_int_t off_i, magma_int_t off_j, magma_int_t nb,
              cuDoubleComplex *a, magma_int_t lda, cuDoubleComplex **dwork, magma_int_t ldda, cudaStream_t stream[][4],
              magma_int_t *info) {

      magma_int_t k;

      if( lapackf77_lsame(uplo, "U") ) {

        /* go through each column */
        magma_int_t j, jj, jb, mj;
        for (j=off_j; j<n; j+=nb) {
          jj = (j-off_j)/(nb*num_gpus);
          k  = ((j-off_j)/nb)%num_gpus;
          magma_setdevice(k);
          jb = min(nb, (n-j));
          if(j+jb < off_j+m) mj = (j-off_i)+jb;
          else mj = m;
          magma_zsetmatrix_async( mj, jb,
                                  A(off_i, j),     lda,
                                  dA(k, 0, jj*nb), ldda, stream[k][0] );
        }
      } else {
        magma_int_t i, ii, ib, ni;

        /* go through each row */
        for(i=off_i; i<m; i+=nb){
          ii = (i-off_i)/(nb*num_gpus);
          k  = ((i-off_i)/nb)%num_gpus;
          magma_setdevice(k);

          ib = min(nb, (m-i));
          if(i+ib < off_i+n) ni = (i-off_i)+ib;
          else ni = n;

          magma_zsetmatrix_async( ib, ni,
                                  A(i, off_j),     lda,
                                  dA(k, ii*nb, 0), ldda, stream[k][2] );
        }
      }
      for( k=0; k<num_gpus; k++ ) {
        magma_setdevice(k);
        magma_queue_sync( stream[k][0] );
      }
      magma_setdevice(0);

      return *info;
}

extern "C" magma_int_t
magma_zdtohpo(int num_gpus, char *uplo, magma_int_t m, magma_int_t n, magma_int_t off_i, magma_int_t off_j, magma_int_t nb, magma_int_t NB,
              cuDoubleComplex *a, magma_int_t lda, cuDoubleComplex **dwork, magma_int_t ldda, cudaStream_t stream[][4],
              magma_int_t *info) {

      magma_int_t k, nk;

      if( lapackf77_lsame(uplo, "U") ) {
        magma_int_t j, jj, jb, mj;

        for (j=off_j+NB; j<n; j+=nb) {
          jj = (j-off_j)/(nb*num_gpus);
          k  = ((j-off_j)/nb)%num_gpus;
          magma_setdevice(k);

          jb = min(nb, (n-j));
          if(j+jb < off_j+m) mj = (j-off_i)+jb;
          else mj = m;
          magma_zgetmatrix_async( mj, jb,
                                  dA(k, 0, jj*nb), ldda,
                                  A(off_i, j),     lda, stream[k][0] );
        }
      } else {
        magma_int_t i, ii, ib, ni;

        /* go through each row */
        for(i=off_i+NB; i<m; i+=nb){
          ii = (i-off_i)/(nb*num_gpus);
          k  = ((i-off_i)/nb)%num_gpus;
          magma_setdevice(k);

          ib = min(nb, (m-i));
          if(i+ib < off_i+n) ni = (i-off_i)+ib;
          else ni = n;

          magma_zgetmatrix_async( ib, ni,
                                  dA(k, ii*nb, 0), ldda,
                                  A(i, off_j),     lda, stream[k][2] );
        }
      }
      for( k=0; k<num_gpus; k++ ) {
        magma_setdevice(k);
        magma_queue_sync( stream[k][0] );
      }
      magma_setdevice(0);

      return *info;
}

#undef A
#undef dA
#undef dlA
#undef dlP
#undef dlAT
#undef dlPT
