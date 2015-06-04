/**
 *
 * @file common.h
 *
 *  MAGMA codelets kernel
 *  MAGMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver,
 *  and INRIA Bordeaux Sud-Ouest
 *
 * @version 2.3.1
 * @author Mathieu Faverge
 * @date 2011-06-01
 *
 **/

/***************************************************************************//**
 *  MAGMA facilities of interest to both src and magmablas directories
 **/
#ifndef _COMMON_H_
#define _COMMON_H_

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <ctype.h>

#if defined( _WIN32 ) || defined( _WIN64 )
#include <io.h>
#else
#include <unistd.h>
#endif

#if defined(MORSE_USE_MPI)
#include <mpi.h>
#endif

#if defined(MORSE_USE_CUDA)
#include <cuda.h>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#include <cublas.h>
#include "magma.h"
#endif

#include <plasma.h>
#include <core_blas.h>

/* Line to avoid conflict with magma, because, we don't know why but lapacke provide a wrong interface of lapack in fortran */
#ifndef LAPACK_NAME
#define LAPACK_NAME(a, b) lapackef77_##a
#endif
#include <lapacke.h>

#include "magma_morse.h"
#include "context.h"
#include "async.h"
#include "descriptor.h"
#include "morse.h"

#include "compute_z.h"
#include "compute_c.h"
#include "compute_d.h"
#include "compute_s.h"

/***************************************************************************//**
 *  Global utilities
 **/
#ifndef max
#define max(a, b) ((a) > (b) ? (a) : (b))
#endif
#ifndef min
#define min(a, b) ((a) < (b) ? (a) : (b))
#endif
#ifndef roundup
#define roundup(a, b) (b <= 0) ? (a) : (((a) + (b)-1) & ~((b)-1))
#endif

/***************************************************************************//**
 *  Global shortcuts
 **/
#define MAGMA_NB          magma->nb
#define MAGMA_IB          magma->ib
#define MAGMA_RHBLK       magma->rhblock
#define MAGMA_TRANSLATION magma->translation
#define MAGMA_PARALLEL    magma->parallel_enabled
#define MAGMA_PROFILING   magma->profiling_enabled
#if defined(MORSE_USE_MPI)
#define MAGMA_MPI_RANK    magma->my_mpi_rank
#define MAGMA_MPI_SIZE    magma->mpi_comm_size
#endif

/***************************************************************************//**
 *  Scheduler properties
 **/
#define PRIORITY        16
#define CALLBACK        17
#define REDUX           18

extern char *plasma_lapack_constants[];

/*
 * Warning and error
 */
#define magma_warning( __f, __msg, ... )     fprintf(stderr, "%s (WARNING): "     __msg,  __f, ##__VA_ARGS__); 
#define magma_error( __f, __msg, ... )       fprintf(stderr, "%s (ERROR): "       __msg,  __f, ##__VA_ARGS__); abort();
#define magma_fatal_error( __f, __msg, ... ) fprintf(stderr, "%s (FATAL ERROR): " __msg,  __f, ##__VA_ARGS__); abort();

void MAGMA_error(const char *func_name, const char *err_message);

#endif /* _COMMON_H_ */

