/*
    -- MAGMA (version 1.2.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2010

       @generated c Tue May 15 18:17:28 2012

*/
#include <math.h>
#include "common_magma.h"

// === Define what BLAS to use ============================================
#define PRECISION_c
#if (defined(PRECISION_s) || defined(PRECISION_d))
  #define magma_cgemm magmablas_cgemm
  /* note: magma_blas_dtrsm seems to have a bug with not-multiple-of-32 N *
   * and rhs > N on Fermi (Pluto)?                                        */
  //#define magma_ctrsm magmablas_ctrsm
#else
  #define magmablas_ctrsm magma_ctrsm
#endif
// === End defining what BLAS to use =======================================

/* to appy pivoting from the previous big-panel: need some index-adjusting */
extern "C" void
magmablas_cpermute_long3( cuFloatComplex *dAT, int lda, int *ipiv, int nb, int ind );


extern "C" magma_int_t
magma_cgetrf1_mgpu(magma_int_t num_gpus, 
         magma_int_t m, magma_int_t n, magma_int_t nb, magma_int_t offset,
         cuFloatComplex **d_lAT, magma_int_t lddat, magma_int_t *ipiv, 
         cuFloatComplex **d_lAP, cuFloatComplex *work, magma_int_t lddwork, 
         cudaStream_t **streaml0, magma_int_t *info)
{
/*  -- MAGMA (version 1.2.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2010

    Purpose
    =======

    CGETRF computes an LU factorization of a general M-by-N matrix A
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

    A       (input/output) COMPLEX array on the GPU, dimension (LDDA,N).
            On entry, the M-by-N matrix to be factored.
            On exit, the factors L and U from the factorization
            A = P*L*U; the unit diagonal elements of L are not stored.

    LDDA     (input) INTEGER
            The leading dimension of the array A.  LDDA >= max(1,M).

    IPIV    (output) INTEGER array, dimension (min(M,N))
            The pivot indices; for 1 <= i <= min(M,N), row i of the
            matrix was interchanged with row IPIV(i).

    INFO    (output) INTEGER
            = 0:  successful exit
            < 0:  if INFO = -i, the i-th argument had an illegal value
                  or another error occured, such as memory allocation failed.
            > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
                  has been completed, but the factor U is exactly
                  singular, and division by zero will occur if it is used
                  to solve a system of equations.
    =====================================================================    */

#define inAT(id,i,j) (d_lAT[(id)] + ((offset)+(i)*nb)*lddat + (j)*nb)

    cuFloatComplex c_one     = MAGMA_C_ONE;
    cuFloatComplex c_neg_one = MAGMA_C_NEG_ONE;

    magma_int_t iinfo, n_local[4];
    magma_int_t maxm, mindim;
    magma_int_t i, d, rows, cols, s, ldpan[4];
    magma_int_t id, i_local, i_local2, nb0, nb1;
    cuFloatComplex *d_panel[4], *panel_local[4];
    static cudaStream_t streaml[4][2];

    /* Check arguments */
    *info = 0;
    if (m < 0)
    *info = -2;
    else if (n < 0)
    *info = -3;
    else if (num_gpus*lddat < max(1,n))
    *info = -5;

    if (*info != 0) {
        magma_xerbla( __func__, -(*info) );
        return *info;
    }

    /* Quick return if possible */
    if (m == 0 || n == 0)
        return *info;

    /* Function Body */
    mindim = min(m, n);
    //nb     = magma_get_cgetrf_nb(m);
    if( num_gpus > ceil((float)n/nb) ) {
      printf( " * too many GPUs for the matrix size, using %d GPUs\n",num_gpus );
      *info = -1;
      return *info;
    }

    {
      /* Use hybrid blocked code. */
      maxm = ((m + 31)/32)*32;

      /* allocate workspace for each GPU */
      for(i=0; i<num_gpus; i++){
        magma_setdevice(i);

        /* local-n and local-ld */
        n_local[i] = ((n/nb)/num_gpus)*nb;
        if (i < (n/nb)%num_gpus)
           n_local[i] += nb;
        else if (i == (n/nb)%num_gpus)
           n_local[i] += n%nb;

        d_panel[i] = &(d_lAP[i][nb*maxm]);

        /* streams */
        magma_queue_create( &streaml[i][0] );
        magma_queue_create( &streaml[i][1] );
      }

      s = mindim / nb;
      for( i=0; i<s; i++ )
            {
                /* Set the GPU number that holds the current panel */
                id = i%num_gpus;
                magma_setdevice(id);

                /* Set the local index where the current panel is */
                i_local = i/num_gpus;
                cols  = maxm - i*nb;
                rows  = m - i*nb;

                /* start sending the panel to cpu */
                magmablas_ctranspose( d_lAP[id], cols, inAT(id,i,i_local), lddat, nb, cols );
                magma_cgetmatrix_async( rows, nb,
                                        d_lAP[id], cols,
                                        work,      lddwork, streaml[id][1] );

                /* make sure that gpu queue is empty */
                magma_device_sync();

                /* the remaining updates */
                if ( i>0 ){
                    /* id-th gpu update the remaining matrix */
                    magmablas_ctrsm( MagmaRight, MagmaUpper, MagmaNoTrans, MagmaUnit, 
                                     n_local[id] - (i_local+1)*nb, nb, 
                                     c_one, panel_local[id],        ldpan[id], 
                                     inAT(id,i-1,i_local+1), lddat );
                    magma_cgemm( MagmaNoTrans, MagmaNoTrans, 
                                 n_local[id]-(i_local+1)*nb, rows, nb, 
                                 c_neg_one, inAT(id,i-1,i_local+1),           lddat, 
                                            &(panel_local[id][nb*ldpan[id]]), ldpan[id], 
                                 c_one,     inAT(id,i,  i_local+1),           lddat );
                }

                /* stnchrnoize i-th panel from id-th gpu into work */
                magma_queue_sync( streaml[id][1] );
                
                /* i-th panel factorization */
                lapackf77_cgetrf( &rows, &nb, work, &lddwork, ipiv+i*nb, &iinfo);
                if ( (*info == 0) && (iinfo > 0) ) {
                    *info = iinfo + i*nb;
                    //break;
                }

                /* start sending the panel to all the gpus */
                for( d=0; d<num_gpus; d++ ) {
                  magma_setdevice(d);
                  magma_csetmatrix_async( rows, nb,
                                          work,     lddwork,
                                          d_lAP[d], maxm, streaml[d][0] );
                }

                for( d=0; d<num_gpus; d++ ) {
                  magma_setdevice(d);
                  /* apply the pivoting */
                  if( d == 0 ) 
                    magmablas_cpermute_long2( inAT(d,0,0), lddat, ipiv, nb, i*nb );
                  else
                    magmablas_cpermute_long3( inAT(d,0,0), lddat, ipiv, nb, i*nb );

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
                  } else { /* no look-ahead for the remaining columns for now */
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
                  magma_queue_sync( streaml[d][0] );
                  magmablas_ctranspose2(panel_local[d], ldpan[d], d_lAP[d], maxm, cols, nb);
                  /* gpu updating the trailing matrix */
                  magmablas_ctrsm( MagmaRight, MagmaUpper, MagmaNoTrans, MagmaUnit, 
                                   nb1, nb, c_one,
                                   panel_local[d],       ldpan[d],
                                   inAT(d, i, i_local2), lddat);
                  magma_cgemm( MagmaNoTrans, MagmaNoTrans, 
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
            if( nb0 > 0 ) {
              magma_setdevice(id);

              /* send the last panel to cpu (no look-ahead for the remaining for remaining columns) */
              magmablas_ctranspose2( d_lAP[id], maxm, inAT(id,s,i_local), lddat, nb0, rows);
              magma_cgetmatrix( rows, nb0, d_lAP[id], maxm, work, lddwork );

              /* make sure that gpu queue is empty */
              magma_device_sync();

              /* factor on cpu */
              lapackf77_cgetrf( &rows, &nb0, work, &lddwork, ipiv+s*nb, &iinfo);
              if ( (*info == 0) && (iinfo > 0) )
              *info = iinfo + s*nb;

              /* start sending the factor to gpus */
              for( d=0; d<num_gpus; d++ ) {
                magma_setdevice(d);
                i_local2 = i_local;
                if( d < id ) i_local2 ++;

                if( d == id || n_local[d] > i_local2*nb ) 
                {
                  magma_csetmatrix_async( rows, nb0,
                                          work,     lddwork,
                                          d_lAP[d], maxm, streaml[d][0] );
                }
              }
            }


            /* clean up */
            for( d=0; d<num_gpus; d++ ) {
              magma_setdevice(d);
              if( nb0 > 0 ) {
                if( d == 0 ) 
                  magmablas_cpermute_long2( inAT(d,0,0), lddat, ipiv, nb0, s*nb );
                else
                  magmablas_cpermute_long3( inAT(d,0,0), lddat, ipiv, nb0, s*nb );

                i_local2 = i_local;
                if( d < id ) i_local2++;
                if( d == id ) {
                  /* the panel belond to this gpu */
                  panel_local[d] = inAT(d,s,i_local);

                  /* next column */
                  nb1 = n_local[d] - i_local*nb-nb0;

                  magma_queue_sync( streaml[d][0] );
                  magmablas_ctranspose2( panel_local[d], lddat, d_lAP[d], maxm, rows, nb0);

                  if( nb1 > 0 )
                  magma_ctrsm( MagmaRight, MagmaUpper, MagmaNoTrans, MagmaUnit, 
                               nb1, nb0, c_one,
                               panel_local[d],        lddat, 
                               inAT(d,s,i_local)+nb0, lddat);
                } else if( n_local[d] > i_local2*nb ) {
                  /* the panel belong to another gpu */
                  panel_local[d] = d_panel[d];

                  /* next column */
                  nb1 = n_local[d] - i_local2*nb;

                  magma_queue_sync( streaml[d][0] );

                  magmablas_ctranspose2( panel_local[d], nb0, d_lAP[d], maxm, rows, nb0);
                  magma_ctrsm( MagmaRight, MagmaUpper, MagmaNoTrans, MagmaUnit, 
                               nb1, nb0, c_one,
                               panel_local[d],     nb0, 
                               inAT(d,s,i_local2), lddat);
                }
              }
              //magma_device_sync();
              magma_queue_destroy( streaml[d][0] );
              magma_queue_destroy( streaml[d][1] );
              
            } /* end of for d=1,..,num_gpus */
    }

    return *info;

    /* End of MAGMA_CGETRF_MGPU */
}

#undef inAT
