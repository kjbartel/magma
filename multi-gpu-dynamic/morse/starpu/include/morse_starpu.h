/**
 *
 * @file morse_starpu.h
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
#ifndef _MORSE_STARPU_H_
#define _MORSE_STARPU_H_

#include <starpu.h>
#include <starpu_profiling.h>
#if defined(MORSE_USE_CUDA)
#include <starpu_scheduler.h>
#include <starpu_cuda.h>
#endif

#include "common.h"
#include "codelets.h"
#include "profiling.h"
#include "codelet_profile.h"
#include "workspace.h"

/***************************************************************************//**
 * StarPU functions
 */
void starpu_Insert_Task(starpu_codelet *cl, ...);
void starpu_unpack_cl_args(void *_cl_arg, ...);

/* 
 * Access to block pointer and leading dimension
 */
#define BLKADDR( desc, type, m, n ) ( (starpu_data_handle)morse_desc_getaddr( desc, m, n ) )

void splagma_set_reduction_methods(starpu_data_handle handle, PLASMA_enum dtyp);

#endif /* _MORSE_STARPU_H_ */
