/*
    -- MAGMA (version 1.2.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       May 2012

       @generated d Tue May 15 18:18:06 2012

*/
#include "common_magma.h"
#define PRECISION_d
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
magmablas_dsetmatrix_transpose3(
                  magma_int_t num_gpus, cudaStream_t **stream0,
                  double  *ha,  int lda, 
                  double **dat, int ldda, int starti,
                  double **dB,  int lddb,
                  int m, int n , int nb)
{
    int d, i = 0, offset = 0, j = 0, j_local, ib;
    static cudaStream_t stream[4][2];

    /* Quick return */
    if ( (m == 0) || (n == 0) )
        return;

    if (lda < m || num_gpus*ldda < n || lddb < m){
        printf("Wrong arguments in zhtodt2 (%d<%d), (%d*%d<%d), or (%d<%d).\n",lda,m,num_gpus,ldda,n,lddb,m);
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
      magma_dsetmatrix_async( m, ib,
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
       magma_dsetmatrix_async( m, ib,
                               ha+i*lda,                        lda,
                               dB[d] + (j_local%2) * nb * lddb, lddb, stream[d][j_local%2] );
  
       /* Make sure that the previous panel (i.e., j%2) has arrived 
          and transpose it directly into the dat matrix                  */
       j_local = (j-num_gpus)/num_gpus;
       //printf( " wait for %d-th panel of %d gpu and store in %d:%d\n\n",j_local%2,d,nb*(starti+j_local),nb*(starti+j_local)+nb-1 );
       magma_queue_sync( stream[d][j_local%2] );
       magmablas_dtranspose2(dat[d]+nb*(starti+j_local), ldda, 
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
      magmablas_dtranspose2(dat[d]+nb*(starti+j_local), ldda, 
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

