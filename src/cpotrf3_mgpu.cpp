/*
    -- MAGMA (version 1.2.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       June 2012

       @generated c Thu Jun 28 12:30:35 2012

*/
#include "common_magma.h"
//#include "trace.h"

/* === Define what BLAS to use ============================================ */
#define PRECISION_c
#if (defined(PRECISION_s) || defined(PRECISION_d)) 
  #define magma_cgemm magmablas_cgemm
  //#define magma_ctrsm magmablas_ctrsm
#endif

#if (GPUSHMEM >= 200)
  #if (defined(PRECISION_s))
     #undef  magma_sgemm
     #define magma_sgemm magmablas_sgemm_fermi80
  #endif
#endif
/* === End defining what BLAS to use ======================================= */
#define A(i, j)  (a   +((j)+off_j)*lda  + (i)+off_i)

#define dlA(id, i, j)     (d_lA[(id)] +               (j)*ldda + (i))
#define dlP(id, i, j, k)  (d_lP[(id)] + (k)*nb*lddp + (j)*lddp + (i))

#define dlAT(id, i, j)    (d_lA[(id)] + (j)*ldda + (i))
#define dlPT(id, i, j, k) (d_lP[(id)] + (k)*nb*lddp + (j)*nb   + (i))

#define VERSION1
extern "C" magma_int_t
magma_cpotrf3_mgpu(int num_gpus, char uplo, magma_int_t m, magma_int_t n, magma_int_t off_i, magma_int_t off_j, magma_int_t nb,
                   cuFloatComplex **d_lA, magma_int_t ldda, cuFloatComplex **d_lP, magma_int_t lddp, 
                   cuFloatComplex *a, magma_int_t lda, cudaStream_t stream[][3], magma_int_t *info ) 
{
/*  -- MAGMA (version 1.2.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       June 2012

    Purpose   
    =======   
    CPOTRF computes the Cholesky factorization of a complex Hermitian   
    positive definite matrix dA.   
    Auxiliary subroutine for cpotrf2_ooc. It is multiple gpu interface to compute 
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

    dA      (input/output) COMPLEX array on the GPU, dimension (LDDA,N)   
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


    magma_int_t     j, jb, nb0, nb2, d, id, j_local, j_local2, buf;
    char            uplo_[2] = {uplo, 0};
    cuFloatComplex c_one     = MAGMA_C_ONE;
    cuFloatComplex c_neg_one = MAGMA_C_NEG_ONE;
    float          d_one     =  1.0;
    float          d_neg_one = -1.0;
    int upper = lapackf77_lsame(uplo_, "U");
    cuFloatComplex *dlpanel;
    magma_int_t n_local[4], ldpanel;
    cudaEvent_t event0[4], /* compute next block -> cherk     */
                event1[4], /* send panel to GPU -> update     */
                event2[4]; /* send block-row to GPU -> update */
    magma_int_t stream1 = 0, stream2 = 1, stream3 = 2;

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

    /* initialization */
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
      magma_event_create( &event0[d] );
      magma_event_create( &event1[d] );
      magma_event_create( &event2[d] );
    }

    /* == initialize the trace */
    //trace_init( 1, num_gpus, 2, (CUstream_st**)stream );

    if (upper) 
    {     
      /* Compute the Cholesky factorization A = U'*U. */
      for (j=0; j<m; j+=nb) {

        /* Set the GPU number that holds the current panel */
        id  = (j/nb)%num_gpus;
        buf = (j/nb)%2;

        /* Set the local index where the current panel is */
        j_local = j/(nb*num_gpus);
        jb = min(nb, (m-j));
          
        /* Update the current diagonal block on stream1 */
        magma_setdevice(id);
        //trace_gpu_start( id, 0, stream[id][0], "syrk", "syrk" );
        magmablasSetKernelStream(stream[id][stream1]);
        magma_cherk(MagmaUpper, MagmaConjTrans, jb, j, 
                    d_neg_one, dlA(id, 0, nb*j_local), ldda, 
                    d_one,     dlA(id, j, nb*j_local), ldda);
        //trace_gpu_end( id, 0, stream[id][0] );

        /* send the diagonal to cpu on stream1 (GPU interface: fist block is on gpu)*/
        //trace_gpu_start( id, 0, stream[id][0], "comm", "D to CPU" );
        magma_cgetmatrix_async( jb, jb,
                                dlA(id, j, nb*j_local), ldda,
                                A(j,j),                 lda, stream[id][stream1] );
        //trace_gpu_end( id, 0, stream[id][0] );

        /* updating the remaing blocks in this column */
        if( j>0 && (j+jb)<n ) {
          if( num_gpus > 1 ) {
            /* wait for the off-diagonal on cpu */
            magma_queue_sync( stream[id][stream3] );
            /* broadcast rows to gpus on stream2 */
            for( d=0; d<num_gpus; d++ ) {
              if( d != id ) {
                magma_setdevice(d);
                //trace_gpu_start( d, 0, stream[d][0], "comm", "row to GPUs" );
                magma_csetmatrix_async( j, jb,
                                        A(0,j),          lda,
                                        dlP(d,jb,0,buf), lddp, stream[d][stream2] );
                //trace_gpu_end( d, 0, stream[d][0] );
                //magma_event_record( event2[d], stream[d][stream2] );
              }
            }
          }

          /* update the remaining blocks of the panel. on stream2 */
          for( d=0; d<num_gpus; d++ ) {

            j_local2 = j_local+1;
            if( d > id ) j_local2 --;
            if( d != id ) {
              dlpanel = dlP(d,jb,0,buf);
              ldpanel = lddp;
            } else {
              dlpanel = dlA(d, 0, nb*j_local);
              ldpanel = ldda;
            }
        
            /* update the panel */
            magma_setdevice(d);
            magmablasSetKernelStream(stream[d][stream2]);
            //trace_gpu_start( d, 0, stream[d][0], "gemm", "gemm" );
            magma_cgemm(MagmaConjTrans, MagmaNoTrans, 
                        jb, (n_local[d]-nb*(j_local2-1)-jb), j, 
                        c_neg_one, dlpanel,                ldpanel, 
                                   dlA(d, 0, nb*j_local2), ldda,
                        c_one,     dlA(d, j, nb*j_local2), ldda);
            magma_event_record( event2[d], stream[d][stream2] );
            //trace_gpu_end( d, 0, stream[d][0] );
          }
        }

        /* wait for panel at cpu */
        magma_setdevice(id);
        magma_queue_sync( stream[id][stream1] );
        /* factor the diagonal */
        //trace_cpu_start( 0, "getrf", "getrf" );
        lapackf77_cpotrf(MagmaUpperStr, &jb, A(j,j), &lda, info);
        //trace_cpu_end( 0 );
        if (*info != 0) {
          *info = *info + j;
          break;
        }

        /* send the diagonal to gpus on stream1 */
        if ( (j+jb) < n) {
          for( d=0; d<num_gpus; d++ ) {
            if( d == id ) {
                dlpanel = dlA(d, j, nb*j_local);
                ldpanel = ldda;
            } else {
                dlpanel = dlP(d,0,0,buf);
                ldpanel = lddp;
            }
            magma_setdevice(d);
            //trace_gpu_start( d, 0, stream[d][0], "comm", "comm" );
            magma_csetmatrix_async( jb, jb,
                                    A(j,j),  lda,
                                    dlpanel, ldpanel, stream[d][stream1] );
            //trace_gpu_end( d, 0, stream[d][0] );
            magma_event_record( event1[d], stream[d][stream1] );
          }
        } else {
          magma_setdevice(id);
          //trace_gpu_start( id, 0, stream[id][0], "comm", "comm" );
          magma_csetmatrix_async( jb, jb,
                                  A(j,j),                 lda,
                                  dlA(id, j, nb*j_local), ldda, stream[id][stream1] );
          //trace_gpu_end( id, 0, stream[id][0] );
        }

        /* panel-factorize the off-diagonal */
        if ( (j+jb) < n) {
          for( d=0; d<num_gpus; d++ ) {
        
            /* next column */
            j_local2 = j_local+1;
            if( d > id ) j_local2--;
            if( d == id ) {
                dlpanel = dlA(d,j,nb*j_local);
                ldpanel = ldda;
            } else {
                dlpanel = dlP(d,0,0,buf);
                ldpanel = lddp;
            }
            nb0 = min(nb, n_local[d]-nb*j_local2 );
        
            magma_setdevice(d);
            if( j+jb < m && d == (j/nb+1)%num_gpus ) { 
              /* owns the next column, look-ahead next block on stream1 */
              magmablasSetKernelStream(stream[d][stream1]);
              magma_queue_wait_event( stream[d][stream1], event2[d] ); //update done
              //trace_gpu_start( d, 1, stream[d][1], "trsm", "trsm" );
              magma_ctrsm( MagmaLeft, MagmaUpper, MagmaConjTrans, MagmaNonUnit, 
                           jb, nb0, c_one,
                           dlpanel,                ldpanel, 
                           dlA(d, j, nb*j_local2), ldda);
              magma_event_record( event0[d], stream[d][stream1] );
              //trace_gpu_end( d, 1, stream[d][1] );

              //if( j+jb < m ) 
              {
                /* send the column to cpu on stream 3                                      */
                /* note: cpotrf2_ooc copies only off-diagonal submatrices to cpu           */
                /* so even on 1gpu, we need to copy to cpu, but don't have to wait for it. */     
                magma_queue_wait_event( stream[d][stream3], event0[d] ); // look-ahead done
                //trace_gpu_start( d, 1, stream[d][1], "comm", "row to CPU" );
                magma_cgetmatrix_async( (j+jb), nb0,
                                        dlA(d, 0, nb*j_local2), ldda,
                                        A(0,j+jb),              lda, stream[d][stream3] );
                //trace_gpu_end( d, 1, stream[d][1] );
              }
            } else {
              /* update all the blocks on stream2 */
              nb2 = n_local[d] - j_local2*nb;
              magmablasSetKernelStream(stream[d][stream2]);
              magma_queue_wait_event( stream[d][stream2], event1[d] ); // panel received
              //trace_gpu_start( d, 0, stream[d][0], "trsm", "trsm" );
              magma_ctrsm( MagmaLeft, MagmaUpper, MagmaConjTrans, MagmaNonUnit, 
                           jb, nb2, c_one,
                           dlpanel,                ldpanel, 
                           dlA(d, j, nb*j_local2), ldda);
              //trace_gpu_end( d, 0, stream[d][0] );
            }
          } /* end of for */

          /* gpu owning the next column                    */
          /* after look ahead, update the remaining blocks */
          if( j+jb < m ) {
            d = (j/nb+1)%num_gpus;
            /* next column */
            j_local2 = j_local+1;
            if( d > id ) j_local2--;
            if( d == id ) {
              dlpanel = dlA(d, j, nb*j_local);
              ldpanel = ldda;
            } else {
              dlpanel = dlP(d,0,0,buf);
              ldpanel = lddp;
            }
            nb0 = min(nb, n_local[d]-nb*j_local2 );
            nb2 = n_local[d] - j_local2*nb - nb0;
        
            magma_setdevice(d);
            /* update the remaining blocks */
            magmablasSetKernelStream(stream[d][stream2]);  
            magma_queue_wait_event( stream[d][stream2], event1[d] ); // panel received
            //trace_gpu_start( d, 0, stream[d][0], "trsm", "trsm" );
            magma_ctrsm( MagmaLeft, MagmaUpper, MagmaConjTrans, MagmaNonUnit, 
                         jb, nb2, c_one,
                         dlpanel,                    ldpanel, 
                         dlA(d, j, nb*j_local2+nb0), ldda);
            //trace_gpu_end( d, 0, stream[d][0] );
          }
        } /* end of ctrsm */
      } /* end of for j=1, .., n */
    } else { 
      /* Compute the Cholesky factorization A = L*L'. */
      for (j=0; j<n; j+=nb) {

        /* Set the GPU number that holds the current panel */
        id = (j/nb)%num_gpus;
        buf = (j/nb)%2;

        /* Set the local index where the current panel is */
        j_local = j/(nb*num_gpus);
        jb = min(nb, (n-j));

        /* Update the current diagonal block on stream1 */
        magma_setdevice(id);
        magmablasSetKernelStream(stream[id][stream1]);
        //trace_gpu_start( id, 0, stream[id][0], "syrk", "syrk" );
        magma_cherk(MagmaLower, MagmaNoTrans, jb, j,
                    d_neg_one, dlAT(id, nb*j_local, 0), ldda,
                    d_one,     dlAT(id, nb*j_local, j), ldda);
        //trace_gpu_end( id, 0, stream[id][0] );

        /* send the diagonal to cpu on stream1 */
        //trace_gpu_start( id, 1, stream[id][1], "comm", "D to CPU" );
        magma_cgetmatrix_async( jb, jb,
                                dlAT(id, nb*j_local, j), ldda,
                                A(j,j),                  lda, stream[id][stream1] );
        //trace_gpu_end( id, 1, stream[id][0] );

        if( j>0 && (j+jb)<m ) {
          if( num_gpus > 1 ) {
            /* wait for off-diagonal row on CPU *
             * and send it to gpus on stream2   */
            magma_queue_sync( stream[id][stream3] );
            for( d=0; d<num_gpus; d++ ) {
              if( d != id ) {
                magma_setdevice(d);
                //trace_gpu_start( d, 0, stream[d][0], "comm", "row to GPU" );
                magma_csetmatrix_async( jb, j,
                                        A(j,0),           lda,
                                        dlPT(d,0,jb,buf), nb, stream[d][stream2] );
                //trace_gpu_end( d, 0, stream[d][0] );
                //magma_event_record( event2[d], stream[d][stream2] );
              }
            }
          }

          /* update the remaining block-rows of the panel */
          for( d=0; d<num_gpus; d++ ) {
            magma_setdevice(d);
            magmablasSetKernelStream(stream[d][stream2]);

            j_local2 = j_local+1;
            if( d > id ) j_local2 --;
            if( d != id ) {
                dlpanel = dlPT(d,0,jb,buf);
                ldpanel = nb;
            } else {
                dlpanel = dlAT(d, nb*j_local, 0);
                ldpanel = ldda;
            }

            /* update the panel */
            //trace_gpu_start( d, 0, stream[d][0], "gemm", "gemm" );
            magma_cgemm( MagmaNoTrans, MagmaConjTrans,
                         n_local[d]-nb*j_local2, jb, j,
                         c_neg_one, dlAT(d, nb*j_local2, 0), ldda,
                                    dlpanel,                 ldpanel,
                         c_one,     dlAT(d, nb*j_local2, j), ldda);
            //trace_gpu_end( d, 0, stream[d][0] );
            magma_event_record( event2[d], stream[d][stream2] );
          }
        }

        /* wait for the panel on cpu */
        magma_setdevice(id);
        magma_queue_sync( stream[id][stream1] );
        /* factor the diagonal */
        //trace_cpu_start( 0, "getrf", "getrf" );
        lapackf77_cpotrf(MagmaLowerStr, &jb, A(j,j), &lda, info);
        //trace_cpu_end( 0 );
        if (*info != 0) {
           *info = *info + j;
           break;
        }

        /* send the diagonal to gpus on stream1 */
        if ( (j+jb) < m) {
          for( d=0; d<num_gpus; d++ ) {
            if( d == id ) {
                dlpanel = dlAT(d, nb*j_local, j);
                ldpanel = ldda;
            } else {
                dlpanel = dlPT(d,0,0,buf);
                ldpanel = nb;
            }
            magma_setdevice(d);
            //trace_gpu_start( d, 0, stream[d][0], "comm", "D to GPU" );
            magma_csetmatrix_async( jb, jb,
                                    A(j,j),  lda,
                                    dlpanel, ldpanel, stream[d][stream1] );
            //trace_gpu_end( d, 0, stream[d][0] );
            magma_event_record( event1[d], stream[d][stream1] );
          }
        } else {
          magma_setdevice(id);
          //trace_gpu_start( id, 0, stream[id][0], "comm", "D to GPU" );
          magma_csetmatrix_async( jb, jb,
                                  A(j,j),                  lda,
                                  dlAT(id, nb*j_local, j), ldda, stream[id][stream1] );
          //trace_gpu_end( id, 0, stream[id][0] );
        }
        if ( (j+jb) < m) {
          for( d=0; d<num_gpus; d++ ) {
            /* next column */
            j_local2 = j_local+1;
            if( d > id ) j_local2--;
            if( d == id ) {
                dlpanel = dlAT(d, nb*j_local, j);
                ldpanel = ldda;
            } else {         
                dlpanel = dlPT(d,0,0,buf);
                ldpanel = nb;
            }
            nb0 = min(nb, n_local[d]-nb*j_local2 );

            magma_setdevice(d);
            if( j+nb < n && d == (j/nb+1)%num_gpus ) { /* owns next column, look-ahead next block on stream1 */
              magmablasSetKernelStream(stream[d][stream1]);
              magma_queue_wait_event( stream[d][stream1], event2[d] ); // update done
              magma_ctrsm( MagmaRight, MagmaLower, MagmaConjTrans, MagmaNonUnit, 
                           nb0, jb, c_one,
                           dlpanel,                 ldpanel, 
                           dlAT(d, nb*j_local2, j), ldda);
              //trace_gpu_end( d, 0, stream[d][0] );
              magma_event_record( event0[d], stream[d][stream1] );

              //if( j+jb < m ) 
              {
                /* send the column to cpu on stream3                                      */
                /* note: cpotrf2_ooc copies only off-diagonal submatrices to cpu,         */
                /* so even on 1gpu, we need to copy to cpu, but don't have to wait for it.*/     
                magma_queue_wait_event( stream[d][stream3], event0[d] ); // lookahead done
                //trace_gpu_start( d, 0, stream[d][0], "comm", "row to CPU" );
                magma_cgetmatrix_async( nb0, j+jb,
                                        dlAT(d, nb*j_local2, 0), ldda,
                                        A(j+jb,0),               lda, stream[d][stream3] );
                //trace_gpu_end( d, 0, stream[d][0] );
              }
            } else { /* other gpus updating all the blocks on stream2 */
              /* update the entire column */
              nb2 = n_local[d] - j_local2*nb;
              magmablasSetKernelStream(stream[d][stream2]);
              magma_queue_wait_event( stream[d][stream2], event1[d] ); // panel received
              //trace_gpu_start( d, 0, stream[d][0], "trsm", "trsm" );
              magma_ctrsm( MagmaRight, MagmaLower, MagmaConjTrans, MagmaNonUnit, 
                           nb2, jb, c_one,
                           dlpanel,                 ldpanel, 
                           dlAT(d, nb*j_local2, j), ldda);
              //trace_gpu_end( d, 0, stream[d][0] );
            }
          } /* end for d */

          /* gpu owing the next column updates remaining blocks on stream2 */
          if( j+nb < n ) {
            d = (j/nb+1)%num_gpus;
            magma_setdevice(d);
            magmablasSetKernelStream(stream[d][stream2]);

            /* next column */
            j_local2 = j_local+1;
            if( d > id ) j_local2--;
            if( d == id ) {
              dlpanel = dlAT(d, nb*j_local, j);
              ldpanel = ldda;
            } else {         
              dlpanel = dlPT(d,0,0,buf);
              ldpanel = nb;
            }
            nb0 = min(nb, n_local[d]-nb*j_local2 );
            nb2 = n_local[d] - j_local2*nb - nb0;

            /* update the remaining blocks in the column */
            magma_queue_wait_event( stream[d][stream2], event1[d] ); // panel received
            //trace_gpu_start( d, 1, stream[d][1], "trsm", "trsm" );
            magma_ctrsm( MagmaRight, MagmaLower, MagmaConjTrans, MagmaNonUnit, 
                         nb2, jb, c_one,
                         dlpanel,                     ldpanel, 
                         dlAT(d, nb*j_local2+nb0, j), ldda);
            //trace_gpu_end( d, 1, stream[d][1] );
          }
        }
      }
    } /* end of else not upper */

    /* == finalize the trace == */
    //trace_finalize( "cpotrf.svg","trace.css" );
    for( d=0; d<num_gpus; d++ ) {
      magma_setdevice(d);
      magma_event_destroy( event0[d] );
      magma_event_destroy( event1[d] );
      magma_event_destroy( event2[d] );
      magmablasSetKernelStream(NULL);
    }
    magma_setdevice(0);

    return *info;
} /* magma_cpotrf_mgpu */

#undef A                     
#define A(i, j)  (a +(j)*lda  + (i))
#define dA(d, i, j) (dwork[(d)]+(j)*ldda + (i))
      
extern "C" magma_int_t
magma_chtodpo(int num_gpus, char *uplo, magma_int_t m, magma_int_t n, magma_int_t off_i, magma_int_t off_j, magma_int_t nb,
              cuFloatComplex *a, magma_int_t lda, cuFloatComplex **dwork, magma_int_t ldda, cudaStream_t stream[][3],
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
          magma_csetmatrix_async( mj, jb,
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
          
          magma_csetmatrix_async( ib, ni,
                                  A(i, off_j),     lda,
                                  dA(k, ii*nb, 0), ldda, stream[k][0] );
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
magma_cdtohpo(int num_gpus, char *uplo, magma_int_t m, magma_int_t n, magma_int_t off_i, magma_int_t off_j, magma_int_t nb, magma_int_t NB,
              cuFloatComplex *a, magma_int_t lda, cuFloatComplex **dwork, magma_int_t ldda, cudaStream_t stream[][3],
              magma_int_t *info) {

      magma_int_t k, stream_id = 1;
      if( lapackf77_lsame(uplo, "U") ) {
        magma_int_t j, jj, jb, mj;

        /* go through each column */
        for (j=off_j+NB; j<n; j+=nb) {
          jj = (j-off_j)/(nb*num_gpus);
          k  = ((j-off_j)/nb)%num_gpus;
          magma_setdevice(k);

          jb = min(nb, (n-j));
          if(j+jb < off_j+m) mj = (j-off_i)+jb;
          else mj = m;
          magma_cgetmatrix_async( mj, jb,
                                  dA(k, 0, jj*nb), ldda,
                                  A(off_i, j),     lda, stream[k][stream_id] );
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

          magma_cgetmatrix_async( ib, ni,
                                  dA(k, ii*nb, 0), ldda,
                                  A(i, off_j),     lda, stream[k][stream_id] );
        }
      }
      for( k=0; k<num_gpus; k++ ) {
        magma_setdevice(k);
        magma_queue_sync( stream[k][stream_id] );
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
