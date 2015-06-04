/*
 *   -- MAGMA (version 1.2.0) --
 *      Univ. of Tennessee, Knoxville
 *      Univ. of California, Berkeley
 *      Univ. of Colorado, Denver
 *      May 2012
 *
 * @author Mark Gates
 */

#include <stdlib.h>
#include <stdio.h>

#include "magma.h"

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
    if ( cudaSuccess != cudaFree( ptr )) {
        return MAGMA_ERR_INVALID_PTR;
    }
    return MAGMA_SUCCESS;
}

// --------------------
// Allocate size bytes on CPU in pinned memory, returning pointer in ptrPtr.
extern "C"
magma_err_t magma_malloc_host( void** ptrPtr, size_t size )
{
    if ( cudaSuccess != cudaMallocHost( ptrPtr, size )) {
        return MAGMA_ERR_HOST_ALLOC;
    }
    return MAGMA_SUCCESS;
}

// --------------------
// Free CPU pinned memory previously allocated by magma_malloc_host.
extern "C"
magma_err_t magma_free_host( void* ptr )
{
    if ( cudaSuccess != cudaFree( ptr )) {
        return MAGMA_ERR_INVALID_PTR;
    }
    return MAGMA_SUCCESS;
}

#endif // HAVE_CUBLAS
