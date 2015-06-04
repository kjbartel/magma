/*
    -- MAGMA (version 1.3.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2012

       @precisions normal z -> s d c

*/
#include "common_magma.h"
#define PRECISION_z
#include "commonblas.h"


//
//      m, n - dimensions in the source (input) matrix.
//             This routine copies the ha matrix from the CPU
//             to dat on the GPU. In addition, the output matrix
//             is transposed. The routine uses a buffer of size
//             2*lddb*nb pointed to by dB (lddb > m) on the GPU. 
//             Note that lda >= m and lddat >= n.
//
extern "C" void 
magmablas_zsetmatrix_transpose( magma_int_t m, magma_int_t n,
                                const cuDoubleComplex  *ha, magma_int_t lda, 
                                cuDoubleComplex       *dat, magma_int_t ldda,
                                cuDoubleComplex        *dB, magma_int_t lddb, magma_int_t nb )
{
    magma_int_t i = 0, j = 0, ib;

    /* Quick return */
    if ( (m == 0) || (n == 0) )
        return;

    if (lda < m || ldda < n || lddb < m){
        printf("Wrong arguments in zhtodt.\n");
        return;
    }

    cudaStream_t stream[2];
    magma_queue_create( &stream[0] );
    magma_queue_create( &stream[1] );
   
    /* Move data from CPU to GPU in the first panel in the dB buffer */
    ib   = min(n-i, nb);
    magma_zsetmatrix_async( m, ib,
                            ha + i*lda,             lda,
                            dB + (j%2) * nb * lddb, lddb, stream[j%2] );
    j++;

    for(i=nb; i<n; i+=nb){
       /* Move data from CPU to GPU in the second panel in the dB buffer */
       ib   = min(n-i, nb);
       magma_zsetmatrix_async( m, ib,
                               ha+i*lda,               lda,
                               dB + (j%2) * nb * lddb, lddb, stream[j%2] );
       j++;
  
       /* Note that the previous panel (i.e., j%2) comes through the stream
          for the kernel so there is no need to synchronize.             */
       // magmablas_ztranspose2( dat+i-nb, ldda, dB + (j%2)*nb*lddb, lddb, m, nb);
       magmablas_ztranspose2s( dat+i-nb, ldda, dB + (j%2)*nb*lddb, lddb, m, nb, &stream[j%2]);
    }

    /* Transpose the last part of the matrix.                            */
    j++;
    // magmablas_ztranspose2( dat+i-nb, ldda, dB + (j%2)*nb*lddb, lddb, m, ib);
    magmablas_ztranspose2s( dat+i-nb, ldda, dB + (j%2)*nb*lddb, lddb, m, ib, &stream[j%2]);

    magma_queue_destroy( stream[0] );
    magma_queue_destroy( stream[1] );
}

//===========================================================================
//  This version is similar to the above but for multiGPUs. The distribution
//  is 1D block cyclic. The input arrays are pointers for the corresponding 
//  GPUs. The streams are passed as argument, in contrast to the single GPU
//  routine.
//  NOTE: see magmablas_zsetmatrix_transpose_mgpu.
//===========================================================================
extern "C" void 
magmablas_zsetmatrix_transpose2( magma_int_t m, magma_int_t n,
                                 const cuDoubleComplex  *ha,  magma_int_t  lda, 
                                 cuDoubleComplex       **dat, magma_int_t *ldda,
                                 cuDoubleComplex       **dB,  magma_int_t  lddb, magma_int_t nb,
                                 magma_int_t num_gpus, cudaStream_t stream[][2] )
{
    magma_int_t i = 0, j[4] = {0, 0, 0, 0}, ib, k = 0;

    /* Quick return */
    if ( (m == 0) || (n == 0) )
        return;

    if (lda < m || lddb < m){
        printf("Wrong arguments in zhtodt2.\n");
        return;
    }

    if (n<num_gpus*nb){
       for(i=0; i<n; i+=nb){
          k = (i/nb)%num_gpus;
          magma_setdevice(k);

          ib = min(n-i, nb);
          magma_zsetmatrix_async( m, ib,
                                  ha+i*lda, lda,
                                  dB[k],    lddb, stream[k][0] );
       }
       for(i=0; i<n; i+=nb){
          k = (i/nb)%num_gpus;
          magma_setdevice(k);

          ib = min(n-i, nb);
          //magma_queue_sync( stream[k][0] );
          //magmablas_ztranspose2( dat[k]+ i/(nb*num_gpus)*nb, ldda[k],
          //                       dB[k], lddb, m, ib);
          magmablas_ztranspose2s( dat[k]+ i/(nb*num_gpus)*nb, ldda[k],
                                 dB[k], lddb, m, ib, &stream[k][0]);
       }
    } 
    else
    {
      for(i=0; i<(n + num_gpus*nb); i+=nb){
         k = (i/nb)%num_gpus;
         magma_setdevice(k);

         if (i<n){
            /* Move data from CPU to GPU in the second panel in the dB buffer */
            ib = min(n-i, nb);
            magma_zsetmatrix_async( m, ib,
                                    ha+i*lda,                 lda,
                                    dB[k] + (j[k]%2)*nb*lddb, lddb, stream[k][j[k]%2] );
         }
         j[k]++;
  
         if (i> (num_gpus-1)*nb){
            /* Make sure that the previous panel (i.e., j[k]%2) has arrived 
               and transpose it directly into the dat matrix                  */
            //magma_queue_sync( stream[k][ j[k]%2 ] );
            ib = min(n - i + num_gpus*nb, nb);
            //magmablas_ztranspose2( dat[k]+ i/(nb*num_gpus)*nb -nb, ldda[k],
            //                       dB[k] +(j[k]%2)*nb*lddb, lddb, m, ib);
            magmablas_ztranspose2s( dat[k]+ i/(nb*num_gpus)*nb -nb, ldda[k],
                                   dB[k] +(j[k]%2)*nb*lddb, lddb, m, ib, &stream[k][j[k]%2]);

         }
      }
    }
}
