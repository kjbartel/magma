/*
    -- MAGMA (version 1.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2011

       @generated c Sun Nov 13 20:48:41 2011

*/
#include "common_magma.h"
#define PRECISION_c
#include "commonblas.h"

//===========================================================================
//  Set a matrix from CPU to multi-GPUs is 1D block cyclic distribution. 
//  The da arrays are pointers to the matrix data for the corresponding GPUs. 
//===========================================================================
extern "C" void 
magmablas_csetmatrix_1D_bcyclic( int m, int n,
                                 cuFloatComplex  *ha, int lda, 
                                 cuFloatComplex  *da[], int ldda, 
                                 int num_gpus, int nb )
{
    int i, k, nk, cdevice;

    cudaGetDevice(&cdevice);

    for(i=0; i<n; i+=nb){
       k = (i/nb)%num_gpus;
       cudaSetDevice(k);
         
       nk = min(nb, n-i);
       //cublasSetMatrix( m, nk, sizeof(cuFloatComplex), ha+i*lda, lda,
       //                 da[k]+i/(nb*num_gpus)*nb*ldda, ldda);
       cudaMemcpy2DAsync(da[k]+i/(nb*num_gpus)*nb*ldda, ldda*sizeof(cuFloatComplex),
                         ha + i*lda, lda*sizeof(cuFloatComplex),
                         sizeof(cuFloatComplex)*m, nk,
                         cudaMemcpyHostToDevice, NULL);
    }

    cudaSetDevice(cdevice);
}


//===========================================================================
//  Get a matrix with 1D block cyclic distribution on multiGPUs to the CPU.
//  The da arrays are pointers to the matrix data for the corresponding GPUs.
//===========================================================================
extern "C" void
magmablas_cgetmatrix_1D_bcyclic( int m, int n,
                                 cuFloatComplex  *da[], int ldda,
                                 cuFloatComplex  *ha, int lda,
                                 int num_gpus, int nb )
{
    int i, k, nk, cdevice;

    cudaGetDevice(&cdevice);

    for(i=0; i<n; i+=nb){
       k = (i/nb)%num_gpus;
       cudaSetDevice(k);

       nk = min(nb, n-i);
       //cublasGetMatrix( m, nk, sizeof(cuFloatComplex),
       //                 da[k]+i/(nb*num_gpus)*nb*ldda, ldda,
       //                 ha+i*lda, lda);
       cudaMemcpy2DAsync(ha + i*lda, lda*sizeof(cuFloatComplex),
                         da[k]+i/(nb*num_gpus)*nb*ldda, ldda*sizeof(cuFloatComplex),
                         sizeof(cuFloatComplex)*m, nk,
                         cudaMemcpyDeviceToHost, NULL);
    }
        
    cudaSetDevice(cdevice);
}

