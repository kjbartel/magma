/*
    -- MAGMA (version 1.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2011

       @precisions normal z -> s d c

*/
#include <math.h>
#include "common_magma.h"

// === Define what BLAS to use ============================================
#define PRECISION_z
#if (defined(PRECISION_s) || defined(PRECISION_d))
  #define cublasZgemm magmablas_zgemm
  #define cublasZtrsm magmablas_ztrsm
#endif
// === End defining what BLAS to use =======================================

extern "C" magma_int_t
magma_zgetrf3(magma_int_t num_gpus, 
              magma_int_t m, magma_int_t n, 
              cuDoubleComplex *a, magma_int_t lda,
              magma_int_t *ipiv, magma_int_t *info)
{
/*  -- MAGMA (version 1.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2011

    Purpose
    =======
    ZGETRF3 computes an LU factorization of a general M-by-N matrix A
    using partial pivoting with row interchanges.

    The factorization has the form
       A = P * L * U
    where P is a permutation matrix, L is lower triangular with unit
    diagonal elements (lower trapezoidal if m > n), and U is upper
    triangular (upper trapezoidal if m < n).

    This is the right-looking Level 3 BLAS version of the algorithm.

    Arguments
    =========
    NUM_GPUS 
            (input) INTEGER
            The number of GPUS to be used for the factorization.

    M       (input) INTEGER
            The number of rows of the matrix A.  M >= 0.

    N       (input) INTEGER
            The number of columns of the matrix A.  N >= 0.

    A       (input/output) COMPLEX_16 array, dimension (LDA,N).
            On entry, the M-by-N matrix to be factored.
            On exit, the factors L and U from the factorization
            A = P*L*U; the unit diagonal elements of L are not stored.

    LDA     (input) INTEGER
            The leading dimension of the array A.  LDA >= max(1,M).

    IPIV    (output) INTEGER array, dimension (min(M,N))
            The pivot indices; for 1 <= i <= min(M,N), row i of the
            matrix was interchanged with row IPIV(i).

    INFO    (output) INTEGER
            = 0:  successful exit
            < 0:  if INFO = -i, the i-th argument had an illegal value
                  if INFO = -7, internal GPU memory allocation failed.
            > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
                  has been completed, but the factor U is exactly
                  singular, and division by zero will occur if it is used
                  to solve a system of equations.
    =====================================================================    */

#define inAT(id,i,j) (d_lAT[(id)] + (i)*nb*lddat + (j)*nb)

    cuDoubleComplex c_one     = MAGMA_Z_ONE;
    cuDoubleComplex c_neg_one = MAGMA_Z_NEG_ONE;

    magma_int_t iinfo, nb, n_local[4], ldat_local[4];
    magma_int_t maxm, mindim;
    magma_int_t i, j, d, rows, cols, s, lddat, lddwork, ldpan[4];
        magma_int_t id, i_local, i_local2, nb0, nb1;
    cuDoubleComplex *d_lAT[4], *d_lAP[4], *d_panel[4], *panel_local[4], *work;
    static cudaStream_t streaml[4][2];

    /* Check arguments */
    *info = 0;
    if (m < 0)
        *info = -2;
    else if (n < 0)
        *info = -3;
    else if (lda < max(1,m))
        *info = -5;

    if (*info != 0) {
        magma_xerbla( __func__, -(*info) );
        return MAGMA_ERR_ILLEGAL_VALUE;
    }

    /* Quick return if possible */
    if (m == 0 || n == 0)
        return MAGMA_SUCCESS;

    mindim = min(m, n);
    nb     = magma_get_zgetrf_nb(m);

    if (nb <= 1 || nb >= n) 
      {
        /* Use CPU code. */
        lapackf77_zgetrf(&m, &n, a, &lda, ipiv, info);
      } 
    else 
      {
          /* Use hybrid blocked code. */
          maxm = ((m + 31)/32)*32;
          if( num_gpus > ceil((double)n/nb) ) {
            printf( " * too many GPUs for the matrix size, using %d GPUs\n",num_gpus );
            *info = -1;
            return MAGMA_ERR_ILLEGAL_VALUE;
          }

          cuDoubleComplex *d_lA[4];

          /* allocate workspace for each GPU */
          for(i=0; i<num_gpus; i++)
            {
              cudaSetDevice(i);

              /* local-n and local-ld */
              n_local[i] = ((n/nb)/num_gpus)*nb;
              if (i < (n/nb)%num_gpus)
                n_local[i] += nb;
              else if (i == (n/nb)%num_gpus)
                n_local[i] += n%nb;
              ldat_local[i] = ((n_local[i]+31)/32)*32;

              /* workspaces */
              if ( CUBLAS_STATUS_SUCCESS != cublasAlloc((2*nb + ldat_local[i])*maxm,
                                                        sizeof(cuDoubleComplex), (void**)&d_lA[i]))
                {
                  for( j=0; j<i; j++ ) {
                    cudaSetDevice(j);
                    cublasFree( d_lA[j] );
                  }
                  return MAGMA_ERR_CUBLASALLOC;
                }

              d_lAP[i]   = d_lA[i];              /* a panel workspace of size nb * maxm */
              d_panel[i] = d_lA[i] +   nb*maxm;  /* another panel of size nb * maxm     */

              lddat = ldat_local[i];
              d_lAT[i]   = d_lA[i] + 2*nb*maxm;  /* local-matrix storage (transpose) of
                                                    size lddat * maxm                   */
              
              panel_local[i] = inAT(i,0,0);
              ldpan[i] = lddat;
              
              /* create the streams */
              if( cudaSuccess != cudaStreamCreate(&streaml[i][0]) ||
                  cudaSuccess != cudaStreamCreate(&streaml[i][1])   ) 
                {
                  for( j=0; j<=i; j++ ) {
                    cudaSetDevice(j);
                    cublasFree( d_lA[j] );
                  }
                  return cudaErrorInvalidValue;
                }
            }
      
          /* Read incrementally from 'a' and transpose in d_lAT; d_lA is work space */
          magmablas_zsetmatrix_transpose2( m, n, a, lda, d_lAT, ldat_local,
              d_lA, maxm, nb, num_gpus, streaml );
           
          /* use 'a' as work space */
          /*
          lddwork = maxm;
          if ( cudaSuccess != cudaMallocHost( (void**)&work, 
                                              lddwork*nb*sizeof(cuDoubleComplex) ) ) 
            {
              for(i=0; i<num_gpus; i++ ) {
                cudaSetDevice(i);
                cublasFree( d_lA[i] );
              }
              return MAGMA_ERR_HOSTALLOC;
            }
          */
          work = a;
          lddwork = lda;

#ifdef ROW_MAJOR_PROFILE
          magma_timestr_t       start, end;
          start = get_current_time();
#endif
          s = mindim / nb;
          for( i=0; i<s; i++ )
            {
              /* Set the GPU number that holds the current panel */
              id = i%num_gpus;
              cudaSetDevice(id);

              /* Set the local index where the current panel is */
              i_local = i/num_gpus;
              cols  = maxm - i*nb;
              rows  = m - i*nb;
              lddat = ldat_local[id];

              if ( i>0 ){
              /* start sending the panel to cpu */
              magmablas_ztranspose( d_lAP[id], cols, inAT(id,i,i_local), lddat, nb, cols );
              //cublasGetMatrix( m-i*nb, nb, sizeof(cuDoubleComplex), d_lAP[id], cols, work, lddwork);
              cudaMemcpy2DAsync( work,     lddwork*sizeof(cuDoubleComplex),
                                 d_lAP[id], cols  *sizeof(cuDoubleComplex),
                                 sizeof(cuDoubleComplex)*rows, nb,
                                 cudaMemcpyDeviceToHost, streaml[id][1]);

              /* make sure that gpu queue is empty */
              //cuCtxSynchronize();

              /* the remaining updates */
              //if ( i>0 ){
                /* id-th gpu update the remaining matrix */
                cublasZtrsm( MagmaRight, MagmaUpper, MagmaNoTrans, MagmaUnit, 
                             n_local[id] - (i_local+1)*nb, nb, 
                             c_one, panel_local[id],        ldpan[id], 
                             inAT(id,i-1,i_local+1), lddat );
                cublasZgemm( MagmaNoTrans, MagmaNoTrans, 
                             n_local[id]-(i_local+1)*nb, rows, nb, 
                             c_neg_one, inAT(id,i-1,i_local+1),           lddat, 
                             &(panel_local[id][nb*ldpan[id]]), ldpan[id], 
                             c_one,     inAT(id,i,  i_local+1),           lddat );
              }

              /* synchronize i-th panel from id-th gpu into work */
              cudaStreamSynchronize(streaml[id][1]);
                
              /* i-th panel factorization */
              lapackf77_zgetrf( &rows, &nb, work, &lddwork, ipiv+i*nb, &iinfo);
              if ( (*info == 0) && (iinfo > 0) ) {
                *info = iinfo + i*nb;
                break;
              }

              /* start sending the panel to all the gpus */
              for( d=0; d<num_gpus; d++ ) {
                cudaSetDevice(d);
                lddat = ldat_local[d];
                //cublasSetMatrix(rows, nb, sizeof(cuDoubleComplex), work, lddwork, d_lAP[d], maxm);
                cudaMemcpy2DAsync(d_lAP[d], maxm   *sizeof(cuDoubleComplex),
                                  work,     lddwork*sizeof(cuDoubleComplex),
                                  sizeof(cuDoubleComplex)*rows, nb,
                                  cudaMemcpyHostToDevice, streaml[d][0]);
              }
              
              for( d=0; d<num_gpus; d++ ) {
                cudaSetDevice(d);
                lddat = ldat_local[d];
                /* apply the pivoting */
                if( d == 0 ) 
                  magmablas_zpermute_long2( d_lAT[d], lddat, ipiv, nb, i*nb );
                else
                  magmablas_zpermute_long3( d_lAT[d], lddat, ipiv, nb, i*nb );
                
                /* storage for panel */
                if( d == id ) {
                  /* the panel belond to this gpu */
                  panel_local[d] = inAT(d,i,i_local);
                  ldpan[d] = lddat;
                  /* next column */
                  i_local2 = i_local+1;
                } else {
                  /* the panel belong to another gpu */
                  panel_local[d] = d_panel[d];
                  ldpan[d] = nb;
                  /* next column */
                  i_local2 = i_local;
                  if( d < id ) i_local2 ++;
                }
                /* the size of the next column */
                if ( s > (i+1) ) {
                  nb0 = nb;
                } else {
                  nb0 = n_local[d]-nb*(s/num_gpus);
                  if( d < s%num_gpus ) nb0 -= nb;
                }
                if( d == (i+1)%num_gpus) {
                  /* owns the next column, look-ahead the column */
                  nb1 = nb0;
                } else {
                  /* update the entire trailing matrix */
                  nb1 = n_local[d] - i_local2*nb;
                }
                
                /* synchronization */
                //cudaStreamSynchronize(streaml[d][0]);
                //cublasSetMatrix(rows, nb, sizeof(cuDoubleComplex), work, lddwork, d_lAP[d], maxm);
                //magmablas_ztranspose2(panel_local[d], ldpan[d], d_lAP[d], maxm, cols, nb);
                magmablas_ztranspose2s(panel_local[d], ldpan[d], d_lAP[d], maxm, cols, nb, 
                                       &streaml[d][0]);
                /* gpu updating the trailing matrix */
                cublasZtrsm( MagmaRight, MagmaUpper, MagmaNoTrans, MagmaUnit, 
                             nb1, nb, c_one,
                             panel_local[d],       ldpan[d],
                             inAT(d, i, i_local2), lddat);
                cublasZgemm( MagmaNoTrans, MagmaNoTrans, 
                             nb1, m-(i+1)*nb, nb, 
                             c_neg_one, inAT(d, i,   i_local2),         lddat,
                             &(panel_local[d][nb*ldpan[d]]), ldpan[d], 
                             c_one,     inAT(d, i+1, i_local2),         lddat );
                
              } /* end of gpu updates */
            } /* end of for i=1..s */
          
          /* Set the GPU number that holds the last panel */
          id = s%num_gpus;
          
          /* Set the local index where the last panel is */
          i_local = s/num_gpus;
          
          /* size of the last diagonal-block */
          nb0 = min(m - s*nb, n - s*nb);
          rows = m    - s*nb;
          cols = maxm - s*nb;
          lddat = ldat_local[id];
          
          if( nb0 > 0 ) {
            /* send the last panel to cpu */
            cudaSetDevice(id);
            //cudaStreamSynchronize(streaml[id][1]);
            magmablas_ztranspose2( d_lAP[id], maxm, inAT(id,s,i_local), lddat, nb0, rows);
            cublasGetMatrix(rows, nb0, sizeof(cuDoubleComplex), d_lAP[id], maxm, work, lddwork);
            
            /* make sure that gpu queue is empty */
            //cuCtxSynchronize();

            /* factor on cpu */
            lapackf77_zgetrf( &rows, &nb0, work, &lddwork, ipiv+s*nb, &iinfo);
            if ( (*info == 0) && (iinfo > 0) )
              *info = iinfo + s*nb;

            /* send the factor to gpus */
            for( d=0; d<num_gpus; d++ ) {
              cudaSetDevice(d);
              lddat = ldat_local[d];
              i_local2 = i_local;
              if( d < id ) i_local2 ++;
              
              if( d == id || n_local[d] > i_local2*nb ) 
                {
                  //cublasSetMatrix(rows, nb0, sizeof(cuDoubleComplex), work, lddwork, d_lAP[d], maxm);
                  cudaMemcpy2DAsync(d_lAP[d], maxm   *sizeof(cuDoubleComplex),
                                    work,     lddwork*sizeof(cuDoubleComplex),
                                    sizeof(cuDoubleComplex)*rows, nb0,
                                    cudaMemcpyHostToDevice, streaml[d][0]);
                }
            }
          }
          
          /* clean up */
          for( d=0; d<num_gpus; d++ ) {
            cudaSetDevice(d);
            lddat = ldat_local[d];
            
            if( nb0 > 0 ) {
              if( d == 0 ) 
                magmablas_zpermute_long2( d_lAT[d], lddat, ipiv, nb0, s*nb );
              else
                magmablas_zpermute_long3( d_lAT[d], lddat, ipiv, nb0, s*nb );
              
              
              i_local2 = i_local;
              if( d < id ) i_local2++;
              if( d == id ) {
                /* the panel belond to this gpu */
                panel_local[d] = inAT(d,s,i_local);
                
                /* next column */
                nb1 = n_local[d] - i_local*nb-nb0;
                
                //cublasSetMatrix(rows, nb0, sizeof(cuDoubleComplex), work, lddwork, d_lAP[d], maxm);
                cudaStreamSynchronize(streaml[d][0]);
                magmablas_ztranspose2( panel_local[d], lddat, d_lAP[d], maxm, rows, nb0);
                
                if( nb1 > 0 )
                  cublasZtrsm( MagmaRight, MagmaUpper, MagmaNoTrans, MagmaUnit, 
                               nb1, nb0, c_one,
                               panel_local[d],        lddat, 
                               inAT(d,s,i_local)+nb0, lddat);
              } else if( n_local[d] > i_local2*nb ) {
                /* the panel belong to another gpu */
                panel_local[d] = d_panel[d];
                
                /* next column */
                nb1 = n_local[d] - i_local2*nb;
                
                //cublasSetMatrix(rows, nb0, sizeof(cuDoubleComplex), work, lddwork, d_lAP[d], maxm);
                cudaStreamSynchronize(streaml[d][0]);
                
                magmablas_ztranspose2( panel_local[d], nb0, d_lAP[d], maxm, rows, nb0);
                cublasZtrsm( MagmaRight, MagmaUpper, MagmaNoTrans, MagmaUnit, 
                             nb1, nb0, c_one,
                             panel_local[d],     nb0, 
                             inAT(d,s,i_local2), lddat);
              }
            }
          } /* end of for d=1,..,num_gpus */

#ifdef ROW_MAJOR_PROFILE
          end = get_current_time();
          printf("\n Performance %f GFlop/s\n", (2./3.*n*n*n /1000000.) / GetTimerValue(start, end));
#endif
          /* save on output */
          magmablas_zgetmatrix_transpose2( m, n, d_lAT, ldat_local, a, lda,
              d_lA, maxm, nb, num_gpus, streaml );
          
          for( d=0; d<num_gpus; d++ ) 
            {
              cublasFree(d_lA[d]);
              cudaStreamDestroy(streaml[d][0]);
              cudaStreamDestroy(streaml[d][1]);
            } /* end of for d=1,..,num_gpus */
          //cudaFreeHost(work);
      }
    
    return MAGMA_SUCCESS;
    
    /* End of MAGMA_ZGETRF3 */
}

#undef inAT
