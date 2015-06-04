/*
    -- MAGMA (version 1.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2011

       @precisions normal z -> s d c

*/
#include "common_magma.h"
#define PRECISION_z
#include "commonblas.h"

extern "C" void
magmablas_ztranspose2s(cuDoubleComplex *odata, int ldo,
                       cuDoubleComplex *idata, int ldi,
                       int m, int n, cudaStream_t *stream );


//
//      m, n - dimensions in the output (ha) matrix.
//             This routine copies the dat matrix from the GPU
//             to ha on the CPU. In addition, the output matrix
//             is transposed. The routine uses a buffer of size
//             2*lddb*nb pointed to by dB (lddb > m) on the GPU. 
//             Note that lda >= m and lddat >= n.
//
extern "C" void 
magmablas_zgetmatrix_transpose( int m, int n,
                                cuDoubleComplex *dat, int ldda,
                                cuDoubleComplex  *ha, int lda,
                                cuDoubleComplex  *dB, int lddb, int nb )
{
    int i = 0, j = 0, ib;

    /* Quick return */
    if ( (m == 0) || (n == 0) )
        return;

    if (lda < m || ldda < n || lddb < m){
        printf("Wrong arguments in zdtoht.\n");
        return;
    }

    static cudaStream_t stream[2];
    cudaStreamCreate(&stream[0]);
    cudaStreamCreate(&stream[1]);

    for(i=0; i<n; i+=nb){
       /* Move data from GPU to CPU using 2 buffers; 1st transpose the data on the GPU */
       ib   = min(n-i, nb);

       //magmablas_ztranspose2 ( dB + (j%2)*nb*lddb, lddb, dat+i, ldda, ib, m);
       magmablas_ztranspose2s( dB + (j%2)*nb*lddb, lddb, dat+i, ldda, ib, m, &stream[j%2]);
       cudaMemcpy2DAsync(ha+i*lda, lda*sizeof(cuDoubleComplex),
                         dB + (j%2) * nb * lddb, lddb*sizeof(cuDoubleComplex),
                         sizeof(cuDoubleComplex)*m, ib, 
                         cudaMemcpyDeviceToHost, stream[j%2]);
       j++;
    }

    cudaStreamDestroy( stream[0] );
    cudaStreamDestroy( stream[1] );
}

//===========================================================================
//  This version is similar to the above but for multiGPUs. The distribution
//  is 1D block cyclic. The input arrays are pointers for the corresponding
//  GPUs. The streams are passed as argument, in contrast to the single GPU
//  routine.
//===========================================================================
extern "C" void
magmablas_zgetmatrix_transpose2( int m, int n,
                                 cuDoubleComplex **dat, int *ldda,
                                 cuDoubleComplex  *ha,  int  lda,
                                 cuDoubleComplex **dB,  int  lddb, int nb,
                                 int num_gpus, cudaStream_t stream[][2] )
{
    int i = 0, j[4] = {0, 0, 0, 0}, ib, k;

    /* Quick return */
    if ( (m == 0) || (n == 0) )
        return;

    if (lda < m || lddb < m){
        printf("Wrong arguments in zdtoht.\n");
        return;
    }

    for(i=0; i<n; i+=nb){
       /* Move data from GPU to CPU using 2 buffers; 1st transpose the data on the GPU */
       k = (i/nb)%num_gpus;
       ib   = min(n-i, nb);
       cudaSetDevice(k);

       //cudaStreamSynchronize(stream[k][j[k]%2]);
       //magmablas_ztranspose2( dB[k] + (j[k]%2)*nb*lddb, lddb, 
       //                       dat[k]+i/(nb*num_gpus)*nb, ldda[k], ib, m);
       magmablas_ztranspose2s(dB[k] + (j[k]%2)*nb*lddb, lddb,
                              dat[k]+i/(nb*num_gpus)*nb, ldda[k], 
                              ib, m, &stream[k][j[k]%2]);
       cudaMemcpy2DAsync(ha+i*lda, lda*sizeof(cuDoubleComplex),
                         dB[k] + (j[k]%2) * nb * lddb, lddb*sizeof(cuDoubleComplex),
                         sizeof(cuDoubleComplex)*m, ib,
                         cudaMemcpyDeviceToHost,
                         stream[k][j[k]%2]);
       j[k]++;
    }
}

