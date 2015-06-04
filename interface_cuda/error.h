#ifndef ERROR_H
#define ERROR_H

#include "common_magma.h"

#include <cuda.h>

// overloaded C++ functions to deal with errors
void magma_xerror( cudaError_t    err, const char* func, const char* file, int line );
void magma_xerror( CUresult       err, const char* func, const char* file, int line );
void magma_xerror( cublasStatus_t err, const char* func, const char* file, int line );
void magma_xerror( magma_err_t    err, const char* func, const char* file, int line );

// cuda provides cudaGetErrorString,
// but not cuGetErrorString or cublasGetErrorString, so provide our own.
const char* cuGetErrorString( CUresult error );
const char* cublasGetErrorString( cublasStatus_t error );
const char* magmaGetErrorString( magma_err_t error );

#ifdef NDEBUG
#define check_error( err )                     ((void)0)
#define check_xerror( err, func, file, line )  ((void)0)
#else
#define check_error( err )                     magma_xerror( err, __func__, __FILE__, __LINE__ )
#define check_xerror( err, func, file, line )  magma_xerror( err, func, file, line )
#endif

#endif        //  #ifndef ERROR_H
