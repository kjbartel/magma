/*
    -- MAGMA (version 1.3.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2012
*/

#ifndef _MAGMABLAS_
#define _MAGMABLAS_

#include <cublas.h>

typedef int  magma_int_t;
typedef int  magma_err_t;
typedef void* magma_devptr;

// For now, make these compatible with old cublas v1 prototypes.
// In the future, we will redefine these data types and
// add queues (opencl queues, cublas handles).
typedef char magma_trans_t;
typedef char magma_side_t ;
typedef char magma_uplo_t ;
typedef char magma_diag_t ;

typedef cudaStream_t magma_stream_t;
typedef cudaStream_t magma_queue_t;
typedef cudaEvent_t  magma_event_t;
typedef int          magma_device_t;

// needed by magmablas*.h, but should eventually go in magma_types.h (see clMAGMA)
#define MagmaMaxGPUs       8

#include "magmablas_z.h"
#include "magmablas_c.h"
#include "magmablas_d.h"
#include "magmablas_s.h"
#include "magmablas_zc.h"
#include "magmablas_ds.h"

#if (GPUSHMEM < 200)  
#define magmablas_zgemm cublasZgemm
#endif
#define magmablas_cgemm cublasCgemm

/** ****************************************************************************
 *  Define magma streams
 */

extern cudaStream_t magma_stream;

#ifdef __cplusplus
extern "C" {
#endif

cublasStatus_t magmablasSetKernelStream( cudaStream_t stream );
cublasStatus_t magmablasGetKernelStream( cudaStream_t *stream );

#ifdef __cplusplus
}
#endif

#endif
