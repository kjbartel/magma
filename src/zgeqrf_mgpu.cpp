/*
    -- MAGMA (version 1.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2011

       @precisions normal z -> s d c

*/
#include "common_magma.h"

#define MultiGPUs

extern "C" magma_int_t
magma_zgeqrf2_mgpu( int num_gpus, magma_int_t m, magma_int_t n,
                    cuDoubleComplex **dlA, magma_int_t ldda,
                    cuDoubleComplex *tau, 
                    magma_int_t *info )
{
/*  -- MAGMA (version 1.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2011

    Purpose
    =======
    ZGEQRF2_MGPU computes a QR factorization of a complex M-by-N matrix A:
    A = Q * R. This is a GPU interface of the routine.

    Arguments
    =========
    M       (input) INTEGER
            The number of rows of the matrix A.  M >= 0.

    N       (input) INTEGER
            The number of columns of the matrix A.  N >= 0.

    dA      (input/output) COMPLEX_16 array on the GPU, dimension (LDDA,N)
            On entry, the M-by-N matrix dA.
            On exit, the elements on and above the diagonal of the array
            contain the min(M,N)-by-N upper trapezoidal matrix R (R is
            upper triangular if m >= n); the elements below the diagonal,
            with the array TAU, represent the orthogonal matrix Q as a
            product of min(m,n) elementary reflectors (see Further
            Details).

    LDDA    (input) INTEGER
            The leading dimension of the array dA.  LDDA >= max(1,M).
            To benefit from coalescent memory accesses LDDA must be
            dividable by 16.

    TAU     (output) COMPLEX_16 array, dimension (min(M,N))
            The scalar factors of the elementary reflectors (see Further
            Details).

    INFO    (output) INTEGER
            = 0:  successful exit
            < 0:  if INFO = -i, the i-th argument had an illegal value
                  if INFO = -9, internal GPU memory allocation failed.

    Further Details
    ===============

    The matrix Q is represented as a product of elementary reflectors

       Q = H(1) H(2) . . . H(k), where k = min(m,n).

    Each H(i) has the form

       H(i) = I - tau * v * v'

    where tau is a complex scalar, and v is a complex vector with
    v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i),
    and tau in TAU(i).
    =====================================================================    */

    #define dlA(gpu,a_1,a_2) ( dlA[gpu]+(a_2)*(ldda) + (a_1))
    #define work_ref(a_1)    ( work + (a_1))
    #define hwork            ( work + (nb)*(m))

    #define hwrk_ref(a_1)    ( local_work + (a_1))
    #define lhwrk            ( local_work + (nb)*(m))

    cuDoubleComplex *dwork[4], *panel[4], *local_work;

    magma_int_t i, j, k, ldwork, lddwork, old_i, old_ib, rows;
    magma_int_t nbmin, nx, ib, nb;
    magma_int_t lhwork, lwork;

    magma_int_t cdevice;
    cudaGetDevice(&cdevice);

    int panel_gpunum, i_local, n_local[4], la_gpu, displacement; 

    *info = 0;
    if (m < 0) {
        *info = -1;
    } else if (n < 0) {
        *info = -2;
    } else if (ldda < max(1,m)) {
        *info = -4;
    }
    if (*info != 0) {
        magma_xerbla( __func__, -(*info) );
        return MAGMA_ERR_ILLEGAL_VALUE;
    }

    k = min(m,n);
    if (k == 0)
        return MAGMA_SUCCESS;

    nb = magma_get_zgeqrf_nb(m);

    displacement = n * nb;
    lwork  = (m+n+64) * nb;
    lhwork = lwork - (m)*nb;

    for(i=0; i<num_gpus; i++){
      #ifdef  MultiGPUs
         cudaSetDevice(i);
      #endif
         if ( CUBLAS_STATUS_SUCCESS != cublasAlloc((n+ldda)*nb,
                                                sizeof(cuDoubleComplex),
                                                (void**)&(dwork[i])) ) {
        *info = -9;
        return MAGMA_ERR_CUBLASALLOC;
      }
    }

    /* Set the number of local n for each GPU */
    for(i=0; i<num_gpus; i++){
      n_local[i] = ((n/nb)/num_gpus)*nb;
      if (i < (n/nb)%num_gpus)
        n_local[i] += nb;
      else if (i == (n/nb)%num_gpus)
        n_local[i] += n%nb;
    }

    if ( cudaSuccess != cudaMallocHost( (void**)&local_work, lwork*sizeof(cuDoubleComplex)) ) {
      *info = -9;
      for(i=0; i<num_gpus; i++){
        #ifdef  MultiGPUs
          cudaSetDevice(i);
        #endif
        cublasFree( dwork[i] );
      }

      return MAGMA_ERR_HOSTALLOC;
    }

    static cudaStream_t streaml[4][2];
    for(i=0; i<num_gpus; i++){
      #ifdef  MultiGPUs
         cudaSetDevice(i);
      #endif
      cudaStreamCreate(&streaml[i][0]);
      cudaStreamCreate(&streaml[i][1]);
    }  

    nbmin = 2;
    nx    = nb;
    ldwork = m;
    lddwork= n;

    if (nb >= nbmin && nb < k && nx < k) {
        /* Use blocked code initially */
        old_i = 0; old_ib = nb;
        for (i = 0; i < k-nx; i += nb) 
          {
            /* Set the GPU number that holds the current panel */
            panel_gpunum = (i/nb)%num_gpus;
            
            /* Set the local index where the current panel is */
            i_local = i/(nb*num_gpus)*nb;
            
            ib = min(k-i, nb);
            rows = m -i;
            /* Send current panel to the CPU */
            #ifdef  MultiGPUs
               cudaSetDevice(panel_gpunum);
            #endif
            cudaMemcpy2DAsync( hwrk_ref(i), ldwork*sizeof(cuDoubleComplex),
                               dlA(panel_gpunum, i, i_local), ldda*sizeof(cuDoubleComplex),
                               sizeof(cuDoubleComplex)*rows, ib,
                               cudaMemcpyDeviceToHost, streaml[panel_gpunum][1]);

            if (i>0){
                /* Apply H' to A(i:m,i+2*ib:n) from the left; this is the look-ahead
                   application to the trailing matrix                                     */
                la_gpu = panel_gpunum;

                /* only the GPU that has next panel is done look-ahead */
                #ifdef  MultiGPUs
                     cudaSetDevice(la_gpu);
                #endif
                   
                magma_zlarfb_gpu( MagmaLeft, MagmaConjTrans, MagmaForward, MagmaColumnwise,
                                  m-old_i, n_local[la_gpu]-i_local-old_ib, old_ib,
                                  panel[la_gpu], ldda, dwork[la_gpu],      lddwork,
                                  dlA(la_gpu, old_i, i_local+old_ib), ldda, 
                                  dwork[la_gpu]+old_ib, lddwork);
                  
                la_gpu = ((i-nb)/nb)%num_gpus;
                #ifdef  MultiGPUs
                cudaSetDevice(la_gpu);
                #endif
                cudaMemcpy2DAsync( panel[la_gpu], ldda  *sizeof(cuDoubleComplex),
                                   hwrk_ref(old_i),  ldwork*sizeof(cuDoubleComplex),
                                   sizeof(cuDoubleComplex)*old_ib, old_ib,
                                   cudaMemcpyHostToDevice, streaml[la_gpu][0]);
            }
            
            #ifdef  MultiGPUs
               cudaSetDevice(panel_gpunum);
            #endif
            cudaStreamSynchronize(streaml[panel_gpunum][1]);

            lapackf77_zgeqrf(&rows, &ib, hwrk_ref(i), &ldwork, tau+i, lhwrk, &lhwork, info);

            // Form the triangular factor of the block reflector
            // H = H(i) H(i+1) . . . H(i+ib-1) 
            lapackf77_zlarft( MagmaForwardStr, MagmaColumnwiseStr,
                              &rows, &ib,
                              hwrk_ref(i), &ldwork, tau+i, lhwrk, &ib);

            zpanel_to_q( MagmaUpper, ib, hwrk_ref(i), ldwork, lhwrk+ib*ib );
            // Send the current panel back to the GPUs 
            // Has to be done with asynchronous copies
            for(j=0; j<num_gpus; j++)
              {  
                #ifdef  MultiGPUs
                   cudaSetDevice(j);
                #endif
                if (j == panel_gpunum)
                  panel[j] = dlA(j, i, i_local);
                else
                  panel[j] = dwork[j]+displacement;
                cudaMemcpy2DAsync(panel[j],    ldda  *sizeof(cuDoubleComplex),
                                  hwrk_ref(i), ldwork*sizeof(cuDoubleComplex),
                                  sizeof(cuDoubleComplex)*rows, ib,
                                  cudaMemcpyHostToDevice, streaml[j][0]);
              }
            for(j=0; j<num_gpus; j++)
              {
                #ifdef  MultiGPUs
                cudaSetDevice(j);
                #endif
                cudaStreamSynchronize(streaml[j][0]);
              }

            /* Restore the panel */
            zq_to_panel( MagmaUpper, ib, hwrk_ref(i), ldwork, lhwrk+ib*ib );

            if (i + ib < n) 
              {
                /* Send the T matrix to the GPU. 
                   Has to be done with asynchronous copies */
                for(j=0; j<num_gpus; j++)
                  {
                    #ifdef  MultiGPUs
                       cudaSetDevice(j);
                    #endif
                       cudaMemcpy2DAsync(dwork[j], lddwork *sizeof(cuDoubleComplex),
                                         lhwrk,    ib      *sizeof(cuDoubleComplex),
                                         sizeof(cuDoubleComplex)*ib, ib,
                                         cudaMemcpyHostToDevice, streaml[j][0]);
                  }

                if (i+nb < k-nx)
                  {
                    /* Apply H' to A(i:m,i+ib:i+2*ib) from the left;
                       This is update for the next panel; part of the look-ahead    */
                    la_gpu = (panel_gpunum+1)%num_gpus;
                    int i_loc = (i+nb)/(nb*num_gpus)*nb;
                    for(j=0; j<num_gpus; j++){
                      #ifdef  MultiGPUs
                      cudaSetDevice(j);
                      #endif
                      //cudaStreamSynchronize(streaml[j][0]);
                      if (j==la_gpu)
                        magma_zlarfb_gpu( MagmaLeft, MagmaConjTrans, MagmaForward, MagmaColumnwise,
                                          rows, ib, ib,
                                          panel[j], ldda, dwork[j],    lddwork,
                                          dlA(j, i, i_loc), ldda, dwork[j]+ib, lddwork);
                      else if (j<=panel_gpunum)
                        magma_zlarfb_gpu( MagmaLeft, MagmaConjTrans, MagmaForward, MagmaColumnwise,
                                          rows, n_local[j]-i_local-ib, ib,
                                          panel[j], ldda, dwork[j],    lddwork,
                                          dlA(j, i, i_local+ib), ldda, dwork[j]+ib, lddwork);
                      else
                        magma_zlarfb_gpu( MagmaLeft, MagmaConjTrans, MagmaForward, MagmaColumnwise,
                                          rows, n_local[j]-i_local, ib,
                                          panel[j], ldda, dwork[j],    lddwork,
                                          dlA(j, i, i_local), ldda, dwork[j]+ib, lddwork);
                    }     
                  }
                else {
                  /* do the entire update as we exit and there would be no lookahead */
                  la_gpu = (panel_gpunum+1)%num_gpus;
                  int i_loc = (i+nb)/(nb*num_gpus)*nb;

                  #ifdef  MultiGPUs
                     cudaSetDevice(la_gpu);
                  #endif
                  magma_zlarfb_gpu( MagmaLeft, MagmaConjTrans, MagmaForward, MagmaColumnwise,
                                    rows, n_local[la_gpu]-i_loc, ib,
                                    panel[la_gpu], ldda, dwork[la_gpu],    lddwork,
                                    dlA(la_gpu, i, i_loc), ldda, dwork[la_gpu]+ib, lddwork);
                  #ifdef  MultiGPUs
                     cudaSetDevice(panel_gpunum);
                  #endif
                  cublasSetMatrix(ib, ib, sizeof(cuDoubleComplex),
                                  hwrk_ref(i), ldwork,
                                  dlA(panel_gpunum, i, i_local),     ldda);
                }
                old_i  = i;
                old_ib = ib;
              }
          }
    } else {
      i = 0;
    }
    
    for(j=0; j<num_gpus; j++){
      #ifdef  MultiGPUs
      cudaSetDevice(j);
      #endif
      cublasFree(dwork[j]);
    }
    
    /* Use unblocked code to factor the last or only block. */
    if (i < k) {
        ib   = n-i;
        rows = m-i;
        lhwork = lwork - rows*ib;

        panel_gpunum = (panel_gpunum+1)%num_gpus;
        int i_loc = (i)/(nb*num_gpus)*nb;

        #ifdef  MultiGPUs
           cudaSetDevice(panel_gpunum);
        #endif
        cublasGetMatrix(rows, ib, sizeof(cuDoubleComplex),
                        dlA(panel_gpunum, i, i_loc), ldda,
                        lhwrk, rows);

        lhwork = lwork - rows*ib;
        lapackf77_zgeqrf(&rows, &ib, lhwrk, &rows, tau+i, lhwrk+ib*rows, &lhwork, info);

        cublasSetMatrix(rows, ib, sizeof(cuDoubleComplex),
                        lhwrk,     rows,
                        dlA(panel_gpunum, i, i_loc), ldda);
    }

    for(i=0; i<num_gpus; i++){
      #ifdef  MultiGPUs
         cudaSetDevice(i);
      #endif
      cudaStreamDestroy(streaml[i][0]);
      cudaStreamDestroy(streaml[i][1]);
    }

    cudaSetDevice(cdevice);

    return MAGMA_SUCCESS;
} /* magma_zgeqrf2_mgpu */
