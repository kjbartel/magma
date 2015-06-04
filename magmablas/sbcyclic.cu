/*
    -- MAGMA (version 1.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2011

       @generated s Sun Nov 13 20:48:40 2011

*/
#include "common_magma.h"
#define PRECISION_s
#include "commonblas.h"

//===========================================================================
//  Set a matrix from CPU to multi-GPUs is 1D block cyclic distribution. 
//  The da arrays are pointers to the matrix data for the corresponding GPUs. 
//===========================================================================
extern "C" void 
magmablas_ssetmatrix_1D_bcyclic( int m, int n,
                                 float  *ha, int lda, 
                                 float  *da[], int ldda, 
                                 int num_gpus, int nb )
{
    int i, k, nk, cdevice;

    cudaGetDevice(&cdevice);

    for(i=0; i<n; i+=nb){
       k = (i/nb)%num_gpus;
       cudaSetDevice(k);
         
       nk = min(nb, n-i);
       //cublasSetMatrix( m, nk, sizeof(float), ha+i*lda, lda,
       //                 da[k]+i/(nb*num_gpus)*nb*ldda, ldda);
       cudaMemcpy2DAsync(da[k]+i/(nb*num_gpus)*nb*ldda, ldda*sizeof(float),
                         ha + i*lda, lda*sizeof(float),
                         sizeof(float)*m, nk,
                         cudaMemcpyHostToDevice, NULL);
    }

    cudaSetDevice(cdevice);
}


//===========================================================================
//  Get a matrix with 1D block cyclic distribution on multiGPUs to the CPU.
//  The da arrays are pointers to the matrix data for the corresponding GPUs.
//===========================================================================
extern "C" void
magmablas_sgetmatrix_1D_bcyclic( int m, int n,
                                 float  *da[], int ldda,
                                 float  *ha, int lda,
                                 int num_gpus, int nb )
{
    int i, k, nk, cdevice;

    cudaGetDevice(&cdevice);

    for(i=0; i<n; i+=nb){
       k = (i/nb)%num_gpus;
       cudaSetDevice(k);

       nk = min(nb, n-i);
       //cublasGetMatrix( m, nk, sizeof(float),
       //                 da[k]+i/(nb*num_gpus)*nb*ldda, ldda,
       //                 ha+i*lda, lda);
       cudaMemcpy2DAsync(ha + i*lda, lda*sizeof(float),
                         da[k]+i/(nb*num_gpus)*nb*ldda, ldda*sizeof(float),
                         sizeof(float)*m, nk,
                         cudaMemcpyDeviceToHost, NULL);
    }
        
    cudaSetDevice(cdevice);
}

