/*
    -- MAGMA (version 1.2.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       May 2012

       @author Stan Tomov
       @generated s Tue May 15 18:18:06 2012
*/
#include "common_magma.h"
#define PRECISION_s
#include "commonblas.h"

//===========================================================================
//  Set a matrix from CPU to multi-GPUs is 1D block cyclic distribution.
//  The dA arrays are pointers to the matrix data for the corresponding GPUs.
//===========================================================================
extern "C" void
magmablas_ssetmatrix_1D_bcyclic( int m, int n,
                                 float  *hA,   int lda,
                                 float  *dA[], int ldda,
                                 int num_gpus, int nb )
{
    int i, d, nk, cdevice;

    magma_getdevice( &cdevice );

    for( i = 0; i < n; i += nb ) {
        d = (i/nb) % num_gpus;
        magma_setdevice( d );
        nk = min(nb, n-i);
        magma_ssetmatrix_async( m, nk,
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
magmablas_sgetmatrix_1D_bcyclic( int m, int n,
                                 float  *dA[], int ldda,
                                 float  *hA,   int lda,
                                 int num_gpus, int nb )
{
    int i, d, nk, cdevice;

    magma_getdevice( &cdevice );

    for( i = 0; i < n; i += nb ) {
        d = (i/nb) % num_gpus;
        magma_setdevice( d );
        nk = min(nb, n-i);
        magma_sgetmatrix_async( m, nk,
                                dA[d] + i/(nb*num_gpus)*nb*ldda, ldda,
                                hA + i*lda, lda, NULL );
    }

    magma_setdevice( cdevice );
}
