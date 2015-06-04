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

#include "magma.h"
#include "error.h"

#ifdef HAVE_CUBLAS

// ========================================
// memory allocation
// Allocate size bytes on GPU, returning pointer in ptrPtr.
extern "C"
magma_err_t magma_malloc( magma_devptr* ptrPtr, size_t size )
{
    if ( cudaSuccess != cudaMalloc( ptrPtr, size )) {
        return MAGMA_ERR_DEVICE_ALLOC;
    }
    return MAGMA_SUCCESS;
}

// --------------------
// Free GPU memory allocated by magma_malloc.
extern "C"
magma_err_t magma_free( magma_devptr ptr )
{
    cudaError_t err = cudaFree( ptr );
    check_error( err );
    if ( err != cudaSuccess ) {
        return MAGMA_ERR_INVALID_PTR;
    }
    return MAGMA_SUCCESS;
}

// --------------------
// Allocate size bytes on CPU, returning pointer in ptrPtr.
// The purpose of using this instead of malloc() is to properly align arrays
// for vector (SSE) instructions. The default implementation uses
// posix_memalign to align memory to a 32 byte boundary.
// This memory can be freed by free().
extern "C"
magma_err_t magma_malloc_cpu( void** ptrPtr, size_t size )
{
#if 1
    int err = posix_memalign( ptrPtr, 32, size );
    if ( err != 0 ) {
        *ptrPtr = NULL;
        return MAGMA_ERR_HOST_ALLOC;
    }
#else
    *ptrPtr = malloc( size );
    if ( ptrPtr == NULL ) {
        return MAGMA_ERR_HOST_ALLOC;
    }
#endif
    return MAGMA_SUCCESS;
}

// --------------------
// Allocate size bytes on CPU in pinned memory, returning pointer in ptrPtr.
extern "C"
magma_err_t magma_malloc_pinned( void** ptrPtr, size_t size )
{
    if ( cudaSuccess != cudaMallocHost( ptrPtr, size )) {
        return MAGMA_ERR_HOST_ALLOC;
    }
    return MAGMA_SUCCESS;
}

// --------------------
// Free CPU pinned memory previously allocated by magma_malloc_pinned.
extern "C"
magma_err_t magma_free_pinned( void* ptr )
{
    cudaError_t err = cudaFreeHost( ptr );
    check_error( err );
    if ( cudaSuccess != err ) {
        return MAGMA_ERR_INVALID_PTR;
    }
    return MAGMA_SUCCESS;
}

#endif // HAVE_CUBLAS
