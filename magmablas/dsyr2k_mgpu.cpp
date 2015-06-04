/*
    -- MAGMA (version 1.2.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       June 2012

       @generated d Thu Jun 28 12:31:23 2012
       @author Mark Gates
*/
#include "common_magma.h"

extern "C"
void magmablas_dsyr2k_mgpu2(
    char uplo, char trans, magma_int_t n, magma_int_t k,
    double alpha, double *dA[], magma_int_t lda,
                           double *dB[], magma_int_t ldb,
    double beta,           double *dC[], magma_int_t ldc,  magma_int_t offset,
    magma_int_t ngpu, magma_int_t nb, cudaStream_t streams[][10], magma_int_t nstream )
{
    #define dA(dev, i, j) (dA[dev] + (i) + (j)*lda)
    #define dB(dev, i, j) (dB[dev] + (i) + (j)*ldb)
    #define dC(dev, i, j) (dC[dev] + (i) + (j)*ldc)
    
    double b = MAGMA_D_MAKE( beta, 0. );
    
    magma_device_t cdev;
    cudaGetDevice( &cdev );
    
    // loop over all blocks
    // Faster to have two loops: first does A*B', second does B*A'
    for( magma_int_t i = 0; i < n; i += nb ) {
        magma_int_t ib     = min( nb, n-i );      // block size
        magma_int_t ioff   = i + offset;          // start global index in parent matrix
        magma_int_t iblock = (ioff / nb) / ngpu;  // local block id
        magma_int_t dev    = (ioff / nb) % ngpu;
        magma_int_t di     = iblock*nb;           // local index in parent matrix
        
        cudaSetDevice( dev );
        magma_int_t s = iblock % nstream;
        magmablasSetKernelStream( streams[ dev ][ s ] );
        
        // C[i:n,i] += A[i:n,0] * B[i,0]'
        //printf( "dgemm  n=%4d, ib=%4d, k=%4d, i=%4d\n", n-i, ib, k, i );
        magma_dgemm( MagmaNoTrans, MagmaTrans, n-i, ib, k,
                     alpha, dA(dev,i,0), lda,
                            dB(dev,i,0), ldb,
                     b,     dC(dev,ioff,di), ldc );
    }
    for( magma_int_t i = 0; i < n; i += nb ) {
        magma_int_t ib     = min( nb, n-i );      // block size
        magma_int_t ioff   = i + offset;          // start global index in parent matrix
        magma_int_t iblock = (ioff / nb) / ngpu;  // local block id
        magma_int_t dev    = (ioff / nb) % ngpu;
        magma_int_t di     = iblock*nb;           // local index in parent matrix
        
        cudaSetDevice( dev );
        magma_int_t s = iblock % nstream;
        magmablasSetKernelStream( streams[ dev ][ s ] );
        
        // C[i:n,i] += B[i:n,0] * A[i,0]'
        //printf( "dgemm  n=%4d, ib=%4d, k=%4d, i=%4d\n", n-i, ib, k, i );
        magma_dgemm( MagmaNoTrans, MagmaTrans, n-i, ib, k,
                     alpha, dB(dev,i,0), ldb,
                            dA(dev,i,0), lda,
                     b,     dC(dev,ioff,di), ldc );
    }
    
    cudaSetDevice( cdev );
}
