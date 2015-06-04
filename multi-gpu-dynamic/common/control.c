/**
 *
 * @file control.c
 *
 *  MAGMA auxiliary routines
 *  MAGMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 1.1.0
 * @author Jakub Kurzak
 * @author Mathieu Faverge
 * @date 2010-11-15
 *
 **/
#include <stdio.h>
#include <stdlib.h>
#include "common.h"

/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  MAGMA_Init - Initialize MAGMA.
 *
 *******************************************************************************
 *
 * @param[in] cores
 *          Number of cores to use (threads to launch).
 *          If cores = 0, cores = MAGMA_NUM_THREADS if it is set, the
 *          system number of core otherwise.
 *
 *******************************************************************************
 *
 * @return
 *          \retval MAGMA_SUCCESS successful exit
 *
 ******************************************************************************/
int MAGMA_Init(int cores, int ncudas)
{
    return MAGMA_InitPar(cores, ncudas, -1);
}

/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  MAGMA_Init_Affinity - Initialize MAGMA.
 *
 *******************************************************************************
 *
 * @param[in] cores
 *          Number of cores to use (threads to launch).
 *          If cores = 0, cores = MAGMA_NUM_THREADS if it is set, the
 *          system number of core otherwise.
 *
 * @param[in] coresbind
 *          Array to specify where to bind each thread.
 *          Each thread i is binded to coresbind[hwloc(i)] if hwloc is
 *          provided, or to coresbind[i] otherwise.
 *          If coresbind = NULL, coresbind = MAGMA_AFF_THREADS if it
 *          is set, the identity function otherwise.
 *
 *******************************************************************************
 *
 * @return
 *          \retval MAGMA_SUCCESS successful exit
 *
 ******************************************************************************/
int MAGMA_InitPar(int nworkers, int ncudas, int nthreads_per_worker)
{
    magma_context_t *magma;

    /* Create context and insert in the context map */
    magma = magma_context_create();
    if (magma == NULL) {
        magma_fatal_error("MAGMA_Init", "magma_context_create() failed");
        return MAGMA_ERR_OUT_OF_RESOURCES;
    }

#if 0    
    /* Init number of cores and topology */
    magma_topology_init();

    /* Set number of nworkers */
    if ( nworkers < 1 ) {
        magma->world_size = magma_get_numthreads();
        if ( magma->world_size == -1 ) {
            magma->world_size = 1;
            magma_warning("MAGMA_Init", "Could not find the number of cores: the thread number is set to 1");
        }
    }
    else
      magma->world_size = nworkers;

    if (magma->world_size <= 0) {
        magma_fatal_error("MAGMA_Init", "failed to get system size");
        return MAGMA_ERR_NOT_FOUND;
    }
    nworkers = magma->world_size;
    
    /* Get the size of each NUMA node */
    magma->group_size = magma_get_numthreads_numa();
    while ( ((magma->world_size)%(magma->group_size)) != 0 ) 
        (magma->group_size)--;
#endif

    morse_init_scheduler( magma, nworkers, ncudas, nthreads_per_worker );

    return MAGMA_SUCCESS;
}

/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  MAGMA_Finalize - Finalize MAGMA.
 *
 *******************************************************************************
 *
 * @return
 *          \retval MAGMA_SUCCESS successful exit
 *
 ******************************************************************************/
int MAGMA_Finalize(void)
{
    magma_context_t *magma = magma_context_self();
    if (magma == NULL) {
        magma_error("MAGMA_Finalize()", "MAGMA not initialized");
        return MAGMA_ERR_NOT_INITIALIZED;
    }
    
    morse_finalize_scheduler( magma );
    magma_context_destroy();

    return MAGMA_SUCCESS;
}

/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  MAGMA_my_mpi_rank - Return the MPI rank of the calling process.
 *
 *******************************************************************************
 *
 *******************************************************************************
 *
 * @return
 *          \retval MPI rank
 *
 ******************************************************************************/
int MAGMA_my_mpi_rank(void)
{
#if defined(MORSE_USE_MPI)
    magma_context_t *magma = magma_context_self();
    if (magma == NULL) {
        magma_error("MAGMA_Finalize()", "MAGMA not initialized");
        return MAGMA_ERR_NOT_INITIALIZED;
    }
    
    return MAGMA_MPI_RANK;
#else
    return 0;
#endif
}
