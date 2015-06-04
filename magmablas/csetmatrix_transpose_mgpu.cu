/*
    -- MAGMA (version 1.2.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       June 2012

       @generated c Thu Jun 28 12:31:24 2012
       @author Ichitaro Yamazaki
*/
#include "common_magma.h"
#define PRECISION_c
#include "commonblas.h"

//
//    m, n - dimensions in the source (input) matrix.
//             This routine copies the ha matrix from the CPU
//             to dat on the GPU. In addition, the output matrix
//             is transposed. The routine uses a buffer of size
//             2*lddb*nb pointed to by dB (lddb > m) on the GPU. 
//             Note that lda >= m and lddat >= n.
//
extern "C" void 
magmablas_csetmatrix_transpose_mgpu(
                  magma_int_t num_gpus, cudaStream_t **stream0,
                  cuFloatComplex  *ha,  magma_int_t lda, 
                  cuFloatComplex **dat, magma_int_t ldda, magma_int_t starti,
                  cuFloatComplex **dB,  magma_int_t lddb,
                  magma_int_t m, magma_int_t n, magma_int_t nb)
{
    int d, i = 0, offset = 0, j = 0, j_local, ib;
    cudaStream_t stream[4][2];

    /* Quick return */
    if ( (m == 0) || (n == 0) )
        return;

    if (lda < m || num_gpus*ldda < n || lddb < m){
        printf( "Wrong arguments in zhtodt2 (%d<%d), (%d*%d<%d), or (%d<%d).\n",
                (int) lda, (int) m, (int) num_gpus, (int) ldda, (int) n, (int) lddb, (int) m );
        return;
    }
    
    for( d=0; d<num_gpus; d++ ) {
      magma_setdevice(d);
      magma_queue_create( &stream[d][0] );
      magma_queue_create( &stream[d][1] );
    }

    /* Move data from CPU to GPU in the first panel in the dB buffer */
    for( d=0; d<num_gpus; d++ ) {
      magma_setdevice(d);
      i  = d * nb;
      ib = min(n-i, nb);

      //printf( " sent A(:,%d:%d) to %d-th panel of %d gpu\n",i,i+ib-1,0,d );
      magma_csetmatrix_async( m, ib,
                              ha + i*lda, lda,
                              dB[d],      lddb, stream[d][0] );
      j++;
    }

    for(i=num_gpus*nb; i<n; i+=nb){
       /* Move data from CPU to GPU in the second panel in the dB buffer */
       d       = j%num_gpus;
       j_local = j/num_gpus;
       magma_setdevice(d);


       ib = min(n-i, nb);
       //printf( " sent A(:,%d:%d) to %d-th panel of %d gpu (%d,%d)\n",i,i+ib-1,j_local%2,d,m,ib );
       magma_csetmatrix_async( m, ib,
                               ha+i*lda,                        lda,
                               dB[d] + (j_local%2) * nb * lddb, lddb, stream[d][j_local%2] );
  
       /* Make sure that the previous panel (i.e., j%2) has arrived 
          and transpose it directly into the dat matrix                  */
       j_local = (j-num_gpus)/num_gpus;
       //printf( " wait for %d-th panel of %d gpu and store in %d:%d\n\n",j_local%2,d,nb*(starti+j_local),nb*(starti+j_local)+nb-1 );
       magma_queue_sync( stream[d][j_local%2] );
       magmablas_ctranspose2(dat[d]+nb*(starti+j_local), ldda, 
                             dB[d] +(j_local%2)*nb*lddb, lddb, 
                             m, nb);
       j++;
       offset += nb;
    }

    /* Transpose the last part of the matrix.                            */
    for( i=0; i<num_gpus; i++ ) {
      d       = j%num_gpus;
      j_local = (j-num_gpus)/num_gpus;
      ib      = min(n-offset, nb);

      magma_setdevice(d);

      //printf( " wait for %d-th panel of %d gpu and store in %d:%d ((%d+%d)*%d)\n\n",j_local%2,d,nb*(starti+j_local),nb*(starti+j_local)+ib-1,starti,j_local,nb );
      magma_queue_sync( stream[d][j_local%2] );
      magmablas_ctranspose2(dat[d]+nb*(starti+j_local), ldda, 
                            dB[d] +(j_local%2)*nb*lddb, lddb, 
                            m, ib);
      j++;
      offset += nb;
    }

    
    for( d=0; d<num_gpus; d++ ) {
      magma_setdevice(d);
      magma_queue_destroy( stream[d][0] );
      magma_queue_destroy( stream[d][1] );
    }
    
}

