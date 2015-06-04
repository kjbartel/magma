/**
 *
 * @file workspace.h
 *
 *  MAGMA codelets kernel
 *  MAGMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver,
 *  and INRIA Bordeaux Sud-Ouest
 *
 * @version 2.3.1
 * @author Cedric Augonnet
 * @date 2011-06-01
 *
 **/

#ifndef _MORSE_STARPU_WORKSPACE_H_
#define _MORSE_STARPU_WORKSPACE_H_

/* 
 * Allocate workspace in host memory: CPU for any worker 
 * or allocate workspace in worker's memory: main memory for cpu workers,
 * and embedded memory for CUDA devices. 
 */
#define MAGMA_HOST_MEM    0
#define MAGMA_WORKER_MEM  1

struct morse_starpu_ws_s {
    size_t size;
    int    memory_location;
    int    which_workers;
    void  *workspaces[STARPU_NMAXWORKERS];
};

typedef struct morse_starpu_ws_s morse_starpu_ws_t;

/*
 * This function creates a workspace on each type of worker in "which_workers"
 * (eg. MAGMA_CUDA|MAGMA_CPU for all CPU and GPU workers).  The
 * memory_location argument indicates whether this should be a buffer in host
 * memory or in worker's memory (MAGMA_HOST_MEM or MAGMA_WORKER_MEM). This function
 * returns 0 upon successful completion. 
 */
int   morse_starpu_ws_alloc   ( morse_starpu_ws_t **workspace, size_t size, int which_workers, int memory_location);
int   morse_starpu_ws_free    ( morse_starpu_ws_t  *workspace);
void *morse_starpu_ws_getlocal( morse_starpu_ws_t  *workspace);

#endif /* _MORSE_STARPU_WORKSPACE_H_ */
