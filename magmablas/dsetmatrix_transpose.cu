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

extern "C" void
magmablas_dtranspose2s(double *odata, int ldo,
                       double *idata, int ldi,
                       int m, int n, cudaStream_t *stream );


//
//      m, n - dimensions in the source (input) matrix.
//             This routine copies the ha matrix from the CPU
//             to dat on the GPU. In addition, the output matrix
//             is transposed. The routine uses a buffer of size
//             2*lddb*nb pointed to by dB (lddb > m) on the GPU. 
//             Note that lda >= m and lddat >= n.
//
extern "C" void 
magmablas_dsetmatrix_transpose( int m, int n,
                                double  *ha, int lda, 
                                double *dat, int ldda,
                                double  *dB, int lddb, int nb )
{
    int i = 0, j = 0, ib;

    /* Quick return */
    if ( (m == 0) || (n == 0) )
        return;

    if (lda < m || ldda < n || lddb < m){
        printf("Wrong arguments in zhtodt.\n");
        return;
    }

    static cudaStream_t stream[2];
    magma_queue_create( &stream[0] );
    magma_queue_create( &stream[1] );
   
    /* Move data from CPU to GPU in the first panel in the dB buffer */
    ib   = min(n-i, nb);
    magma_dsetmatrix_async( m, ib,
                            ha + i*lda,             lda,
                            dB + (j%2) * nb * lddb, lddb, stream[j%2] );
    j++;

    for(i=nb; i<n; i+=nb){
       /* Move data from CPU to GPU in the second panel in the dB buffer */
       ib   = min(n-i, nb);
       magma_dsetmatrix_async( m, ib,
                               ha+i*lda,               lda,
                               dB + (j%2) * nb * lddb, lddb, stream[j%2] );
       j++;
  
       /* Note that the previous panel (i.e., j%2) comes through the stream
          for the kernel so there is no need to synchronize.             */
       // magmablas_dtranspose2( dat+i-nb, ldda, dB + (j%2)*nb*lddb, lddb, m, nb);
       magmablas_dtranspose2s( dat+i-nb, ldda, dB + (j%2)*nb*lddb, lddb, m, nb, &stream[j%2]);
    }

    /* Transpose the last part of the matrix.                            */
    j++;
    // magmablas_dtranspose2( dat+i-nb, ldda, dB + (j%2)*nb*lddb, lddb, m, ib);
    magmablas_dtranspose2s( dat+i-nb, ldda, dB + (j%2)*nb*lddb, lddb, m, ib, &stream[j%2]);

    magma_queue_destroy( stream[0] );
    magma_queue_destroy( stream[1] );
}

//===========================================================================
//  This version is similar to the above but for multiGPUs. The distribution
//  is 1D block cyclic. The input arrays are pointers for the corresponding 
//  GPUs. The streams are passed as argument, in contrast to the single GPU
//  routine.
//===========================================================================
extern "C" void 
magmablas_dsetmatrix_transpose2( int m, int n,
                                 double  *ha,  int  lda, 
                                 double **dat, int *ldda,
                                 double **dB,  int  lddb, int nb,
                                 int num_gpus, cudaStream_t stream[][2] )
{
    int i = 0, j[4] = {0, 0, 0, 0}, ib, k = 0;

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
          magma_dsetmatrix_async( m, ib,
                                  ha+i*lda, lda,
                                  dB[k],    lddb, stream[k][0] );
       }
       for(i=0; i<n; i+=nb){
          k = (i/nb)%num_gpus;
          magma_setdevice(k);

          ib = min(n-i, nb);
          //magma_queue_sync( stream[k][0] );
          //magmablas_dtranspose2( dat[k]+ i/(nb*num_gpus)*nb, ldda[k],
          //                       dB[k], lddb, m, ib);
          magmablas_dtranspose2s( dat[k]+ i/(nb*num_gpus)*nb, ldda[k],
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
            magma_dsetmatrix_async( m, ib,
                                    ha+i*lda,                 lda,
                                    dB[k] + (j[k]%2)*nb*lddb, lddb, stream[k][j[k]%2] );
         }
         j[k]++;
  
         if (i> (num_gpus-1)*nb){
            /* Make sure that the previous panel (i.e., j[k]%2) has arrived 
               and transpose it directly into the dat matrix                  */
            //magma_queue_sync( stream[k][ j[k]%2 ] );
            ib = min(n - i + num_gpus*nb, nb);
            //magmablas_dtranspose2( dat[k]+ i/(nb*num_gpus)*nb -nb, ldda[k],
            //                       dB[k] +(j[k]%2)*nb*lddb, lddb, m, ib);
            magmablas_dtranspose2s( dat[k]+ i/(nb*num_gpus)*nb -nb, ldda[k],
                                   dB[k] +(j[k]%2)*nb*lddb, lddb, m, ib, &stream[k][j[k]%2]);

         }
      }
    }
}
