/*
    -- MAGMA (version 1.2.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2010

       @precisions normal z -> s d c

*/
#include "common_magma.h"
#define PRECISION_z
#include "commonblas.h"

//
//    m, n - dimensions in the output (ha) matrix.
//             This routine copies the dat matrix from the GPU
//             to ha on the CPU. In addition, the output matrix
//             is transposed. The routine uses a buffer of size
//             2*lddb*nb pointed to by dB (lddb > m) on the GPU. 
//             Note that lda >= m and lddat >= n.
//
extern "C" void 
magmablas_zgetmatrix_transpose3(
                  magma_int_t num_gpus, cudaStream_t **stream0,
                  cuDoubleComplex **dat, int ldda,
                  cuDoubleComplex   *ha, int lda,
                  cuDoubleComplex  **dB, int lddb,
                  int m, int n , int nb)
{
    int i = 0, j = 0, j_local, d, ib;
    static cudaStream_t stream[4][2];

    /* Quick return */
    if ( (m == 0) || (n == 0) )
        return;

    if (lda < m || num_gpus*ldda < n || lddb < m){
        printf("Wrong arguments in zdtoht (%d<%d), (%d*%d<%d), or (%d<%d).\n",lda,m,num_gpus,ldda,n,lddb,m);
        return;
    }
    
    for( d=0; d<num_gpus; d++ ) {
      magma_setdevice(d);
      magma_queue_create( &stream[d][0] );
      magma_queue_create( &stream[d][1] );
    }
    

    for(d=0; d<num_gpus; d++ ) {
       magma_setdevice(d);
       i  = nb*d;
       ib = min(n-i, nb);

       /* transpose and send the first tile */
       magmablas_ztranspose2( dB[d], lddb, dat[d], ldda, ib, m);
       //printf( " (%dx%d) (%d,%d)->(%d,%d) (%d,%d)\n",m,ib,0,i,0,0,lda,lddb );
       magma_zgetmatrix_async( m, ib,
                               dB[d],    lddb,
                               ha+i*lda, lda, stream[d][0] );
       j++;
    }

    for(i=num_gpus*nb; i<n; i+=nb){
       d = j%num_gpus;
       magma_setdevice(d);

       ib = min(n-i, nb);

       /* Move data from GPU to CPU using two buffers; first transpose the data on the GPU */
       j_local = j/num_gpus;
       magmablas_ztranspose2( dB[d] +(j_local%2)*nb*lddb, lddb, 
                              dat[d]+nb*j_local,          ldda, 
                              ib, m);
       //printf( " (%dx%d) (%d,%d)->(%d,%d) (%d,%d)\n",m,ib,0,i,0,(j_local%2)*nb,lda,lddb );
       magma_zgetmatrix_async( m, ib,
                               dB[d] + (j_local%2)*nb*lddb, lddb,
                               ha+i*lda,                    lda, stream[d][j_local%2] );

       /* wait for the previous tile */
       j_local = (j-num_gpus)/num_gpus;
       magma_queue_sync( stream[d][j_local%2] );
       j++;
    }

    for( i=0; i<num_gpus; i++ ) {
       d = j%num_gpus;
       magma_setdevice(d);

       j_local = (j-num_gpus)/num_gpus;
       magma_queue_sync( stream[d][j_local%2] );
       j++;
    }

    
    for( d=0; d<num_gpus; d++ ) {
      magma_setdevice(d);
      magma_queue_destroy( stream[d][0] );
      magma_queue_destroy( stream[d][1] );
    }
    
}



