/**
 *
 *  @file codelet_zpotrf.c
 *
 *  MAGMA codelets kernel
 *  MAGMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver,
 *  and INRIA Bordeaux Sud-Ouest
 *
 *  @version 2.3.1
 *  @author Mathieu Faverge
 *  @author Cedric Augonnet
 *  @date 2011-06-01
 *  @precisions normal z -> c d s
 *
 **/
#include "morse_starpu.h"

#define PRECISION_z

CL_CALLBACK(zgemm,  starpu_matrix_get_nx(task->buffers[2].handle),     2. *size*size*size);
CL_CALLBACK(zgeqrt, starpu_matrix_get_nx(task->buffers[0].handle), (4./3.)*size*size*size);
CL_CALLBACK(zgessm, starpu_matrix_get_nx(task->buffers[2].handle),     2. *size*size*size);
CL_CALLBACK(zgetrl, starpu_matrix_get_nx(task->buffers[0].handle), (2./3.)*size*size*size);
CL_CALLBACK(zherk,  starpu_matrix_get_nx(task->buffers[0].handle), (   1.+size)*size*size);
CL_CALLBACK(zlacpy, starpu_matrix_get_nx(task->buffers[0].handle),              size*size);
#if defined(PRECISION_z) || defined(PRECISION_c)
CL_CALLBACK(zplghe, starpu_matrix_get_nx(task->buffers[0].handle),              size*size);
#endif
CL_CALLBACK(zplgsy, starpu_matrix_get_nx(task->buffers[0].handle),              size*size);
CL_CALLBACK(zplrnt, starpu_matrix_get_nx(task->buffers[0].handle),              size*size);
CL_CALLBACK(zpotrf, starpu_matrix_get_nx(task->buffers[0].handle), (1./3.)*size*size*size);
CL_CALLBACK(zssssm, starpu_matrix_get_nx(task->buffers[0].handle),     2. *size*size*size);
CL_CALLBACK(ztrsm,  starpu_matrix_get_nx(task->buffers[0].handle),         size*size*size);
CL_CALLBACK(ztsmqr, starpu_matrix_get_nx(task->buffers[0].handle), (4.0*size+starpu_matrix_get_nx(task->buffers[3].handle))*size*size); 
CL_CALLBACK(ztsqrt, starpu_matrix_get_nx(task->buffers[0].handle),     2. *size*size*size);
CL_CALLBACK(ztstrf, starpu_matrix_get_nx(task->buffers[0].handle),         size*size*size);
CL_CALLBACK(zunmqr, starpu_matrix_get_nx(task->buffers[0].handle),     2. *size*size*size);

/* TODO : fix the following macro */
//CL_CALLBACK(zttqrt, starpu_matrix_get_nx(task->buffers[0].handle), 3.0*size*size*size);
//CL_CALLBACK(zttmqr, starpu_matrix_get_nx(task->buffers[0].handle), 4.0*size*size*size + (starpu_matrix_get_nx(task->buffers[3].handle))*size*size); 
//CL_CALLBACK(zcopy_lower_tile, starpu_matrix_get_nx(task->buffers[0].handle), 0);

/* First formula according to equivalent GEMM, second one is the real one */
//CL_CALLBACK(zssssm, starpu_matrix_get_nx(task->buffers[0].handle), size*size*(2.*size+starpu_matrix_get_nx(task->buffers[2].handle))); 
