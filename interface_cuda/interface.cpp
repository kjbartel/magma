/*
 *   -- MAGMA (version 1.2.1) --
 *      Univ. of Tennessee, Knoxville
 *      Univ. of California, Berkeley
 *      Univ. of Colorado, Denver
 *      June 2012
 *
 * @author Mark Gates
 */

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "common_magma.h"
#include "error.h"

#ifdef HAVE_CUBLAS

// ========================================
// initialization
// --------------------
extern "C"
void magma_init()
{
}

// --------------------
extern "C"
void magma_finalize()
{
}


// ========================================
// device support
// --------------------
extern "C"
void magma_getdevices(
    magma_device_t* devices,
    magma_int_t     size,
    magma_int_t*    numPtr )
{
    cudaError_t err;
    int cnt;
    err = cudaGetDeviceCount( &cnt );
    check_error( err );
    
    cnt = min( cnt, size );
    for( int i = 0; i < cnt; ++i ) {
        devices[i] = i;
    }
    *numPtr = cnt;
}

// --------------------
extern "C"
void magma_getdevice( magma_device_t* device )
{
    cudaError_t err;
    err = cudaGetDevice( device );
    check_error( err );
}

// --------------------
extern "C"
void magma_setdevice( magma_device_t device )
{
    cudaError_t err;
    err = cudaSetDevice( device );
    check_error( err );
}

// --------------------
extern "C"
void magma_device_sync()
{
    cudaError_t err;
    err = cudaDeviceSynchronize();
    check_error( err );
}


// ========================================
// queue support
// At the moment, MAGMA queue == CUDA stream.
// In the future, MAGMA queue may be CUBLAS handle.
// --------------------
extern "C"
void magma_queue_create( /*magma_device_t device,*/ magma_queue_t* queuePtr )
{
    //cudaStream_t   stream;
    //cublasStatus_t stat;
    cudaError_t    err;
    //err  = cudaSetDevice( device );
    //stat = cublasCreate( queuePtr );
    err  = cudaStreamCreate( queuePtr );  //&stream );
    //stat = cublasSetStream( *queuePtr, stream );
    check_error( err );
}

// --------------------
extern "C"
void magma_queue_destroy( magma_queue_t queue )
{
    //cudaStream_t   stream;
    //cublasStatus_t stat;
    cudaError_t    err;
    //stat = cublasGetStream( queue, &stream );
    err  = cudaStreamDestroy( queue );  //stream );
    //stat = cublasDestroy( queue );
    check_error( err );
}

// --------------------
extern "C"
void magma_queue_sync( magma_queue_t queue )
{
    //cudaStream_t   stream;
    //cublasStatus_t stat;
    cudaError_t    err;
    //stat = cublasGetStream( queue, &stream );
    err  = cudaStreamSynchronize( queue );  //stream );
    check_error( err );
}


// ========================================
// event support
// --------------------
extern "C"
void magma_event_create( magma_event_t* event )
{
    cudaError_t err;
    err = cudaEventCreate( event );
    check_error( err );
}

// --------------------
extern "C"
void magma_event_destroy( magma_event_t event )
{
    cudaError_t err;
    err = cudaEventDestroy( event );
    check_error( err );
}

// --------------------
extern "C"
void magma_event_record( magma_event_t event, magma_queue_t queue )
{
    cudaError_t err;
    err = cudaEventRecord( event, queue );
    check_error( err );
}

// --------------------
// blocks CPU until event occurs
extern "C"
void magma_event_sync( magma_event_t event )
{
    cudaError_t err;
    err = cudaEventSynchronize( event );
    check_error( err );
}

// --------------------
// blocks queue (but not CPU) until event occurs
extern "C"
void magma_queue_wait_event( magma_queue_t queue, magma_event_t event )
{
    cudaError_t err;
    err = cudaStreamWaitEvent( queue, event, 0 );
    check_error( err );
}

#endif // HAVE_CUBLAS
