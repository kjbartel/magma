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

//===========================================================================
#include <sys/time.h>
#include <assert.h>

double get_current_cpu_time(void)
{
  struct timeval  time_val;

  gettimeofday(&time_val, NULL);

  return (double)(time_val.tv_sec) + (double)(time_val.tv_usec) / 1000000.0;
}

#define MAX_THREADS 5

//#define MAX_EVENTS 163840
#define MAX_EVENTS 1048576

int    event_num        [MAX_THREADS]               __attribute__ ((aligned (128)));
double event_start_time [MAX_THREADS]               __attribute__ ((aligned (128)));
double event_end_time   [MAX_THREADS]               __attribute__ ((aligned (128)));
double event_log        [MAX_THREADS][MAX_EVENTS]   __attribute__ ((aligned (128)));
int log_events = 1;

#define core_cpu_event_start(my_core_id)                      \
  event_start_time[my_core_id] = get_current_cpu_time();      \

#define core_cpu_event_end(my_core_id)                        \
  event_end_time[my_core_id] = get_current_cpu_time();        \

#define core_gpu_event_start(my_core_id, e1, e2)              \
  cudaEventElapsedTime(&ctime, e1, e2);                       \
  event_start_time[my_core_id] = ctime/1000.+dtime;           \

#define core_gpu_event_end(my_core_id, e1, e2)                \
  cudaEventElapsedTime(&ctime, e1, e2);                       \
  event_end_time[my_core_id] = ctime/1000.+dtime;             \

#define core_log_event(event, my_core_id)             \
  event_log[my_core_id][event_num[my_core_id]+0] = my_core_id;\
  event_log[my_core_id][event_num[my_core_id]+1] = event_start_time[my_core_id];\
  event_log[my_core_id][event_num[my_core_id]+2] = event_end_time[my_core_id];\
  event_log[my_core_id][event_num[my_core_id]+3] = (event);\
  event_num[my_core_id] += (log_events << 2);         \
  event_num[my_core_id] &= (MAX_EVENTS-1);

void dump_trace(int cores_num)
{
  char trace_file_name[32];
  FILE *trace_file;
  int event;
  int core;

  double scale = 100000.0;
  double large = 100.0;

  sprintf(trace_file_name, "trace.svg");
  trace_file = fopen(trace_file_name, "w");
  assert(trace_file != NULL);

  fprintf(trace_file,
          "<?xml version=\"1.0\" standalone=\"no\"?>"
          "<svg version=\"1.1\" baseProfile=\"full\" xmlns=\"http://www.w3.org/2000/svg\" "
          "xmlns:xlink=\"http://www.w3.org/1999/xlink\" xmlns:ev=\"http://www.w3.org/2001/xml-events\"  "
          ">\n"
          "  <g font-size=\"20\">\n");   
  
  for (core = 0; core < cores_num; core++)
    for (event = 0; event < event_num[core]; event += 4)
      {
        int    thread = event_log[core][event+0];
        double start  = event_log[core][event+1];
        double end    = event_log[core][event+2];
        int    color  = event_log[core][event+3];

        start -= event_log[core][2];
        end   -= event_log[core][2];
        /*
            fprintf(trace_file,
                "    "
                "<rect x=\"%.2lf\" y=\"%.0lf\" width=\"%.2lf\" height=\"%.0lf\" "
                "fill=\"#%06x\" stroke=\"#000000\" stroke-width=\"1\"/>\n",
                start * scale,
                thread * 100.0,
                (end - start) * scale,
                90.0,
                color);
        */
        fprintf(trace_file,
            "    "
            "<rect x=\"%.2lf\" y=\"%.0lf\" width=\"%.2lf\" height=\"%.0lf\" "
                //            "fill=\"#%06x\" />\n",
                "fill=\"#%06x\" stroke=\"#000000\" stroke-width=\"1\"/>\n",
                start * scale,
                thread * (large+20.0),
                (end - start) * scale,
                large,
                color);
      }

  fprintf(trace_file,
          "  </g>\n"
          "</svg>\n");

  fclose(trace_file);
}

//===========================================================================

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

    float ctime, dtime;

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
    cudaEvent_t start[4], stop[4][10];
    for(i=0; i<num_gpus; i++){
      #ifdef  MultiGPUs
         cudaSetDevice(i);
      #endif
      cudaStreamCreate(&streaml[i][0]);
      cudaStreamCreate(&streaml[i][1]);

      cudaEventCreate( &start[i] );
      for(j=0; j<10; j++)
        cudaEventCreate( &stop[i][j]  );

      cudaEventRecord( start[i], 0 );
    }  

    core_cpu_event_start(num_gpus);
    core_cpu_event_end(num_gpus);
    core_log_event(0x666666, num_gpus);
    
    dtime = get_current_cpu_time();

    for(j=0; j<num_gpus; j++){
      cudaSetDevice(j);
      cudaEventRecord(stop[j][0], 0);
      cudaEventRecord(stop[j][1], 0);
      cudaThreadSynchronize();
      core_gpu_event_start(j, start[j], stop[j][0]);
      core_gpu_event_end(j, start[j], stop[j][1]);
      core_log_event(0x666666, j);
    }
      
    nbmin = 2;
    nx    = nb;
    ldwork = m;
    lddwork= n;

    // cublasSetKernelStream(streaml[0][0]);

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
                cudaEventRecord(stop[la_gpu][0], 0);
                magma_zlarfb_gpu( MagmaLeft, MagmaConjTrans, MagmaForward, MagmaColumnwise,
                                  m-old_i, n_local[la_gpu]-i_local-old_ib, old_ib,
                                  panel[la_gpu], ldda, dwork[la_gpu],      lddwork,
                                  dlA(la_gpu, old_i, i_local+old_ib), ldda, 
                                  dwork[la_gpu]+old_ib, lddwork);
                cudaEventRecord(stop[la_gpu][1], 0);

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
            core_cpu_event_start(num_gpus);
            lapackf77_zgeqrf(&rows, &ib, hwrk_ref(i), &ldwork, tau+i, lhwrk, &lhwork, info);

            // Form the triangular factor of the block reflector
            // H = H(i) H(i+1) . . . H(i+ib-1) 
            lapackf77_zlarft( MagmaForwardStr, MagmaColumnwiseStr,
                              &rows, &ib,
                              hwrk_ref(i), &ldwork, tau+i, lhwrk, &ib);

            zpanel_to_q( MagmaUpper, ib, hwrk_ref(i), ldwork, lhwrk+ib*ib );
            core_cpu_event_end(num_gpus);
            core_log_event(0x006680, num_gpus);

            core_cpu_event_start(num_gpus);
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
            core_cpu_event_end(num_gpus);
            core_log_event(0xDD0000, num_gpus);

            //=================== take the values of all counters
            /* Restore the panel */
            core_cpu_event_start(num_gpus);
            zq_to_panel( MagmaUpper, ib, hwrk_ref(i), ldwork, lhwrk+ib*ib );
            core_cpu_event_end(num_gpus);
            core_log_event(0x006680, num_gpus);
            
            if (i>0){
              for(j=0; j<num_gpus; j++){
                cudaSetDevice(j);
                core_gpu_event_start(j, start[j], stop[j][2]);
                core_gpu_event_end(j, start[j], stop[j][3]);
                core_log_event(0x880000, j);

                core_gpu_event_start(j, start[j], stop[j][3]);
                core_gpu_event_end(j, start[j], stop[j][0]);
                core_log_event(0xDD0000, j);
              }
            }
            cudaSetDevice(panel_gpunum);
            core_gpu_event_start(panel_gpunum, start[panel_gpunum], stop[panel_gpunum][0]);
            core_gpu_event_end(panel_gpunum, start[panel_gpunum], stop[panel_gpunum][1]);
            core_log_event(0x660000, panel_gpunum);
            

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
                      cuCtxSynchronize();

                      cudaEventRecord(stop[j][2], 0);
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
                      cudaEventRecord(stop[j][3], 0);
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
    dump_trace(num_gpus+1);

    return MAGMA_SUCCESS;
} /* magma_zgeqrf2_mgpu */
