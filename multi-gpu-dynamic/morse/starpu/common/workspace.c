/**
 *
 * @file workspace.c
 *
 *  MAGMA codelets kernel
 *  MAGMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver,
 *  and INRIA Bordeaux Sud-Ouest
 *
 * @version 2.3.1
 * @author Cedric Augonnet
 * @author Mathieu Faverge
 * @date 2011-06-01
 *
 **/
#include "morse_starpu.h"

static void allocate_workspace_on_workers(void *arg)
{
    struct morse_starpu_ws_s *workspace = arg;
  
    int id;
    enum starpu_archtype type;
    
    id = starpu_worker_get_id();
    type = starpu_worker_get_type(id);
    
#ifdef MORSE_USE_CUDA
    if (type == STARPU_CUDA_WORKER)
    {
        int memory_location = workspace->memory_location;
        
        if (memory_location == MAGMA_HOST_MEM)
        {
            /* Use pinned memory because the kernel is very likely
             * to transfer these data between the CPU and the GPU.
             * */
            cudaMallocHost(&workspace->workspaces[id], workspace->size);
        }
        else { 
            /* Allocate on the device */
            cudaMalloc(&workspace->workspaces[id], workspace->size);
        }
    }
    else
#endif
    {
        /* This buffer should only be used within the CPU kernel, so
         * there is no point in using pinned memory here. */
        workspace->workspaces[id] = malloc(workspace->size);
    }
        
    assert(workspace->workspaces[id]);
}


static void free_workspace_on_workers(void *arg)
{
    struct morse_starpu_ws_s *workspace = arg;
    
    int id;
    enum starpu_archtype type;
    
    id = starpu_worker_get_id();
    type = starpu_worker_get_type(id);
    
#ifdef MORSE_USE_CUDA
    if (type == STARPU_CUDA_WORKER)
    {
        int memory_location = workspace->memory_location;
        
        if (memory_location == MAGMA_HOST_MEM)
        {
            cudaFreeHost(workspace->workspaces[id]);
        }
        else {
            cudaFree(workspace->workspaces[id]);
        }
    }
    else
#endif
    {
        free(workspace->workspaces[id]);
    }
    
    workspace->workspaces[id] = NULL;
}

/*
 * This function creates a workspace on each type of worker in "which_workers"
 * (eg. MAGMA_CUDA|MAGMA_CPU for all CPU and GPU workers).  The
 * memory_location argument indicates whether this should be a buffer in host
 * memory or in GPU memory (MAGMA_HOST_MEM or MAGMA_GPU_MEM). This function
 * returns 0 upon successful completion.: 
 */
int morse_starpu_ws_alloc(morse_starpu_ws_t **workspace,
                          size_t size, int which_workers, int memory_location)
{
    if (!workspace)
        return -EINVAL;
    
    struct morse_starpu_ws_s *descr = calloc(1, sizeof(struct morse_starpu_ws_s));
    
    *workspace = descr;
    
    if (!descr)
        return -ENOMEM;
    
    descr->size = size;
    descr->memory_location = memory_location;
    descr->which_workers = which_workers;
    
    starpu_execute_on_each_worker(allocate_workspace_on_workers, descr, which_workers);
    
    return 0;
}

int morse_starpu_ws_free(morse_starpu_ws_t *workspace)
{
    if (!workspace)
        return -EINVAL;
    
    starpu_execute_on_each_worker(free_workspace_on_workers, workspace, workspace->which_workers);
    
    free(workspace);
    
    return 0;
}

void *morse_starpu_ws_getlocal(morse_starpu_ws_t *workspace)
{
    struct morse_starpu_ws_s *descr = workspace;
    int id = starpu_worker_get_id();
    return descr->workspaces[id];
}
