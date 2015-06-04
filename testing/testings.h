#ifndef TESTINGS_H
#define TESTINGS_H

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>

#include "magma.h"

#ifndef min
#define min(a,b)  (((a)<(b))?(a):(b))
#endif

#ifndef max
#define max(a,b)  (((a)<(b))?(b):(a))
#endif


#define TESTING_CUDA_INIT()                                                \
    if( CUBLAS_STATUS_SUCCESS != cublasInit() ) {                          \
        fprintf(stderr, "ERROR: cublasInit failed\n");                     \
        exit(-1);                                                          \
    }                                                                      \
    printout_devices();


#define TESTING_CUDA_FINALIZE()                                            \
    cublasShutdown();


#define TESTING_CUDA_INIT_MGPU()                                           \
{                                                                          \
    int ndevices;                                                          \
    cudaGetDeviceCount( &ndevices );                                       \
    for( int idevice = 0; idevice < ndevices; ++idevice ) {                \
        cudaSetDevice(idevice);                                            \
        if( CUBLAS_STATUS_SUCCESS != cublasInit() ) {                      \
            fprintf(stderr, "ERROR: cublasInit failed\n");                 \
            exit(-1);                                                      \
        }                                                                  \
    }                                                                      \
    cudaSetDevice(0);                                                      \
    printout_devices();                                                    \
}


#define TESTING_CUDA_FINALIZE_MGPU()                                       \
{                                                                          \
    int ndevices;                                                          \
    cudaGetDeviceCount( &ndevices );                                       \
    for( int idevice = 0; idevice < ndevices; ++idevice ) {                \
        cudaSetDevice(idevice);                                            \
        cublasShutdown();                                                  \
    }                                                                      \
}


#define TESTING_MALLOC( ptr, type, size )                                  \
    if ( MAGMA_SUCCESS !=                                                  \
            magma_malloc_cpu( (void**) &ptr, (size)*sizeof(type) )) {      \
        fprintf( stderr, "!!!! malloc failed for: %s\n", #ptr );           \
        exit(-1);                                                          \
    }


#define TESTING_HOSTALLOC( ptr, type, size )                                  \
    if ( MAGMA_SUCCESS !=                                                     \
            magma_malloc_pinned( (void**) &ptr, (size)*sizeof(type) )) {      \
        fprintf( stderr, "!!!! magma_malloc_pinned failed for: %s\n", #ptr ); \
        exit(-1);                                                             \
    }


#define TESTING_DEVALLOC( ptr, type, size )                                \
    if ( MAGMA_SUCCESS !=                                                  \
            magma_malloc( (void**) &ptr, (size)*sizeof(type) )) {          \
        fprintf( stderr, "!!!! magma_malloc failed for: %s\n", #ptr );     \
        exit(-1);                                                          \
    }


#define TESTING_FREE(ptr)                                                  \
    magma_free_cpu(ptr);


#define TESTING_HOSTFREE(ptr)                                              \
    magma_free_pinned( ptr );


#define TESTING_DEVFREE(ptr)                                               \
    magma_free( ptr );


#ifdef __cplusplus
extern "C" {
#endif

void magma_zhermitian( magma_int_t N, cuDoubleComplex* A, magma_int_t lda );
void magma_chermitian( magma_int_t N, cuFloatComplex*  A, magma_int_t lda );
void magma_dhermitian( magma_int_t N, double*          A, magma_int_t lda );
void magma_shermitian( magma_int_t N, float*           A, magma_int_t lda );

void magma_zhpd( magma_int_t N, cuDoubleComplex* A, magma_int_t lda );
void magma_chpd( magma_int_t N, cuFloatComplex*  A, magma_int_t lda );
void magma_dhpd( magma_int_t N, double*          A, magma_int_t lda );
void magma_shpd( magma_int_t N, float*           A, magma_int_t lda );

void magma_assert( bool condition, const char* msg, ... );

#ifdef __cplusplus
}
#endif

#endif /* TESTINGS_H */
