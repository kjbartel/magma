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

//===========================================================================
//  Set a matrix from CPU to multi-GPUs is 1D block cyclic distribution. 
//  The da arrays are pointers to the matrix data for the corresponding GPUs. 
//===========================================================================
extern "C" void 
magmablas_zsetmatrix_1D_bcyclic( int m, int n,
                                 cuDoubleComplex  *ha, int lda, 
                                 cuDoubleComplex  *da[], int ldda, 
                                 int num_gpus, int nb )
{
    int i, k, nk, cdevice;

    cudaGetDevice(&cdevice);

    for(i=0; i<n; i+=nb){
       k = (i/nb)%num_gpus;
       cudaSetDevice(k);
         
       nk = min(nb, n-i);
       //cublasSetMatrix( m, nk, sizeof(cuDoubleComplex), ha+i*lda, lda,
       //                 da[k]+i/(nb*num_gpus)*nb*ldda, ldda);
       cudaMemcpy2DAsync(da[k]+i/(nb*num_gpus)*nb*ldda, ldda*sizeof(cuDoubleComplex),
                         ha + i*lda, lda*sizeof(cuDoubleComplex),
                         sizeof(cuDoubleComplex)*m, nk,
                         cudaMemcpyHostToDevice, NULL);
    }

    cudaSetDevice(cdevice);
}


//===========================================================================
//  Get a matrix with 1D block cyclic distribution on multiGPUs to the CPU.
//  The da arrays are pointers to the matrix data for the corresponding GPUs.
//===========================================================================
extern "C" void
magmablas_zgetmatrix_1D_bcyclic( int m, int n,
                                 cuDoubleComplex  *da[], int ldda,
                                 cuDoubleComplex  *ha, int lda,
                                 int num_gpus, int nb )
{
    int i, k, nk, cdevice;

    cudaGetDevice(&cdevice);

    for(i=0; i<n; i+=nb){
       k = (i/nb)%num_gpus;
       cudaSetDevice(k);

       nk = min(nb, n-i);
       //cublasGetMatrix( m, nk, sizeof(cuDoubleComplex),
       //                 da[k]+i/(nb*num_gpus)*nb*ldda, ldda,
       //                 ha+i*lda, lda);
       cudaMemcpy2DAsync(ha + i*lda, lda*sizeof(cuDoubleComplex),
                         da[k]+i/(nb*num_gpus)*nb*ldda, ldda*sizeof(cuDoubleComplex),
                         sizeof(cuDoubleComplex)*m, nk,
                         cudaMemcpyDeviceToHost, NULL);
    }
        
    cudaSetDevice(cdevice);
}

