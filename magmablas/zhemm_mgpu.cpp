/*
    -- MAGMA (version 1.2.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       June 2012

       @precisions normal z -> s d c
       @author Mark Gates
       
       This still has poor performance. Work in progress.
*/
#include "common_magma.h"
#include <assert.h>

extern "C" void
magmablas_zsymmetrize( char uplo, magma_int_t m, cuDoubleComplex *dA, magma_int_t ldda );

extern "C"
void magmablas_zhemm_mgpu(
    char side, char uplo, magma_int_t m, magma_int_t n,
    cuDoubleComplex alpha, cuDoubleComplex *dA[], magma_int_t ldda,  magma_int_t offset,
                           cuDoubleComplex *dB[], magma_int_t lddb,
    cuDoubleComplex beta,  cuDoubleComplex *dC[], magma_int_t lddc,
                           cuDoubleComplex *C,    magma_int_t ldc,
    magma_int_t ngpu, magma_int_t nb, cudaStream_t streams[][20], magma_int_t nstream )
{
    #define dA(dev, i, j) (dA[dev] + (i) + (j)*ldda)
    #define dB(dev, i, j) (dB[dev] + (i) + (j)*lddb)
    #define dC(dev, i, j) (dC[dev] + (i) + (j)*lddc)
    #define C(i, j) (C + (i) + (j)*ldc)
    
    assert( ldda >= m );
    assert( lddb >= m );
    assert( lddc >= m );
    
    cuDoubleComplex c_one  = MAGMA_Z_ONE;
    cuDoubleComplex c_zero = MAGMA_Z_ZERO;
    magma_int_t ione = 1;
    
    magma_device_t cdev;
    magma_getdevice( &cdev );

    // create events for sync
    magma_event_t event[ MagmaMaxGPUs ][ 20 ];
    for( magma_int_t d = 0; d < ngpu; ++d ) {
        magma_setdevice( d );
        cudaMemset( dC[d], 0, lddc*n*sizeof(cuDoubleComplex) );
        for( magma_int_t s = 0; s < nstream; ++s ) {
            cudaEventCreateWithFlags( &event[d][s], cudaEventDisableTiming );
        }
    }
    
    // loop over all blocks
    // Faster to have several loops:
    // first  symmetrizes A[i,i]
    // second does C[i]      = A[i:m,  i]'*B[i:m]
    // third  does C[i+1:m] += A[i+1:m,i] *B[i]
    for( magma_int_t i = 0; i < m; i += nb ) {
        magma_int_t ib     = min( nb, m-i );      // block size
        magma_int_t ioff   = i + offset;          // start global index in parent matrix
        magma_int_t iblock = (ioff / nb) / ngpu;  // local block id
        magma_int_t dev    = (ioff / nb) % ngpu;
        magma_int_t di     = iblock*nb;           // local index in parent matrix
        
        magma_setdevice( dev );
        magma_int_t s = iblock % nstream;
        magmablasSetKernelStream( streams[ dev ][ s ] );
        
        // make diagonal block symmetric
        magmablas_zsymmetrize( MagmaLower, ib, dA(dev,ioff,di), ldda );
    }
    for( magma_int_t i = 0; i < m; i += nb ) {
        magma_int_t ib     = min( nb, m-i );      // block size
        magma_int_t ioff   = i + offset;          // start global index in parent matrix
        magma_int_t iblock = (ioff / nb) / ngpu;  // local block id
        magma_int_t dev    = (ioff / nb) % ngpu;
        magma_int_t di     = iblock*nb;           // local index in parent matrix
        
        magma_setdevice( dev );
        magma_int_t s = iblock % nstream;
        magmablasSetKernelStream( streams[ dev ][ s ] );

        // C[i] = A[i:m,i]' * B[i:m0]
        //printf( "gemm1: A[%4d,%4d]*B[%4d] -> C[%4d] ib     %4d, n %4d, m-i %4d\n",
        //        ioff, di, i, i, ib, n, m-i );
        magma_zgemm( MagmaConjTrans, MagmaNoTrans, ib, n, m-i,
                     alpha,  dA(dev,ioff,di), ldda,
                             dB(dev,i,0),     lddb,
                     c_zero, dC(dev,i,0),     lddc );
        magma_event_record( event[dev][s], streams[dev][s] );
    }
    for( magma_int_t i = 0; i < m; i += nb ) {
        magma_int_t ib     = min( nb, m-i );      // block size
        magma_int_t ioff   = i + offset;          // start global index in parent matrix
        magma_int_t iblock = (ioff / nb) / ngpu;  // local block id
        magma_int_t dev    = (ioff / nb) % ngpu;
        magma_int_t di     = iblock*nb;           // local index in parent matrix
        
        // these have to be on same stream, since they write into same block,
        // unless we used separate C workspace for each stream
        magma_setdevice( dev );
        magmablasSetKernelStream( streams[dev][0] );
        for( magma_int_t s = 0; s < nstream; ++s ) {
            magma_queue_wait_event( streams[dev][0], event[dev][s] );
        }
        
        // C[i+1:n] += A[i+1:n,i] * B[i]
        //printf( "gemm2: A[%4d,%4d]*B[%4d] -> C[%4d] m-i-ib %4d, n %4d, ib  %4d\n",
        //        ioff+ib, di, i, i+ib, m-i-ib, n, ib );
        magma_zgemm( MagmaNoTrans, MagmaNoTrans, m-i-ib, n, ib,
                     alpha, dA(dev,ioff+ib,di), ldda,
                            dB(dev,i,0),        lddb,
                     c_one, dC(dev,i+ib,0),     lddc );
    }
    
    // meanwhile on CPU, scale C := beta*C
    for( magma_int_t j = 0; j < n; ++j ) {
        blasf77_zscal( &m, &beta, C(0,j), &ione );
    }
    
    // wait and reduce results
    magma_int_t size = ldc*n;
    cuDoubleComplex *Ctmp = C(0,n);
    for( magma_int_t d = 0; d < ngpu; ++d ) {
        magma_setdevice( d );
        magma_zgetmatrix( m, n, dC[d], lddc, Ctmp, ldc );
        blasf77_zaxpy( &size, &c_one, Ctmp, &ione, C, &ione );
    }
    
    magma_setdevice( cdev );
}
