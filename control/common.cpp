/**
 *
 * @file common.cpp
 *
 *  MAGMA (version 1.1) --
 *  Univ. of Tennessee, Knoxville
 *  Univ. of California, Berkeley
 *  Univ. of Colorado, Denver
 *  November 2011
 *
 **/
#include "common_magma.h"

cudaStream_t magma_stream = 0;

cublasStatus_t magmablasSetKernelStream( cudaStream_t stream )
{
/*  -- MAGMA (version 1.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2011

    Purpose
    =======

    magmablasSetKernelStream sets the CUDA stream that all MAGMA BLAS and
    CUBLAS routines use.

    Arguments
    =========

    stream  (input) cudaStream_t
            The CUDA stream.

    =====================================================================   */
    magma_stream = stream;
    return cublasSetKernelStream( stream );
}


cublasStatus_t magmablasGetKernelStream( cudaStream_t *stream )
{
/*  -- MAGMA (version 1.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2011

    Purpose
    =======

    magmablasSetKernelStream gets the CUDA stream that all MAGMA BLAS
    routines use.

    Arguments
    =========

    stream  (output) cudaStream_t
            The CUDA stream.

    =====================================================================   */
    *stream = magma_stream;
    return CUBLAS_STATUS_SUCCESS;
}
