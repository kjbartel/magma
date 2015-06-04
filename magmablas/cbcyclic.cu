/*
    -- MAGMA (version 1.3.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2012

       @author Stan Tomov
       @generated c Wed Nov 14 22:53:53 2012
*/
#include "common_magma.h"
#define PRECISION_c
#include "commonblas.h"

//===========================================================================
//  Set a matrix from CPU to multi-GPUs is 1D block cyclic distribution.
//  The dA arrays are pointers to the matrix data for the corresponding GPUs.
//===========================================================================
extern "C" void
magmablas_csetmatrix_1D_bcyclic( magma_int_t m, magma_int_t n,
                                 const cuFloatComplex *hA,   magma_int_t lda,
                                 cuFloatComplex       *dA[], magma_int_t ldda,
                                 magma_int_t num_gpus, magma_int_t nb )
{
    magma_int_t i, d, nk;
    magma_device_t cdevice;

    magma_getdevice( &cdevice );

    for( i = 0; i < n; i += nb ) {
        d = (i/nb) % num_gpus;
        magma_setdevice( d );
        nk = min(nb, n-i);
        magma_csetmatrix_async( m, nk,
                                hA + i*lda, lda,
                                dA[d] + i/(nb*num_gpus)*nb*ldda, ldda, NULL );
    }

    magma_setdevice( cdevice );
}


//===========================================================================
//  Get a matrix with 1D block cyclic distribution on multiGPUs to the CPU.
//  The dA arrays are pointers to the matrix data for the corresponding GPUs.
//===========================================================================
extern "C" void
magmablas_cgetmatrix_1D_bcyclic( magma_int_t m, magma_int_t n,
                                 cuFloatComplex  *dA[], magma_int_t ldda,
                                 cuFloatComplex  *hA,   magma_int_t lda,
                                 magma_int_t num_gpus, magma_int_t nb )
{
    magma_int_t i, d, nk;
    magma_device_t cdevice;

    magma_getdevice( &cdevice );

    for( i = 0; i < n; i += nb ) {
        d = (i/nb) % num_gpus;
        magma_setdevice( d );
        nk = min(nb, n-i);
        magma_cgetmatrix_async( m, nk,
                                dA[d] + i/(nb*num_gpus)*nb*ldda, ldda,
                                hA + i*lda, lda, NULL );
    }

    magma_setdevice( cdevice );
}
