#ifndef TESTINGS_H
#define TESTINGS_H

#ifndef min
#define min(a,b)  (((a)<(b))?(a):(b))
#endif

#ifndef max
#define max(a,b)  (((a)<(b))?(b):(a))
#endif


#define TESTING_CUDA_INIT()                                                \
    if( CUBLAS_STATUS_SUCCESS != cublasInit() ) {                          \
        fprintf(stderr, "CUBLAS: Not initialized\n");                      \
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
            fprintf(stderr, "CUBLAS: Not initialized\n");                  \
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
    free(ptr);


#define TESTING_HOSTFREE(ptr)                                              \
    magma_free_pinned( ptr );


#define TESTING_DEVFREE(ptr)                                               \
    magma_free( ptr );

#endif /* TESTINGS_H */
