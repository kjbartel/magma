/**
 *
 * @file workspace.c
 *
 *  MAGMA auxiliary routines
 *  MAGMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.4.5
 * @author Mathieu Faverge
 * @date 2010-11-15
 *
 **/
#include <stdlib.h>
#include "common.h"

/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  MAGMA_Dealloc_Handle - Deallocate workspace descriptor allocated by
 *  any workspace allocation routine.
 *
 *******************************************************************************
 *
 * @param[in] desc
 *          Workspace descriptor
 *
 *******************************************************************************
 *
 * @return
 *          \retval MAGMA_SUCCESS successful exit
 *
 ******************************************************************************/
int MAGMA_Dealloc_Workspace(magma_desc_t **desc)
{
    magma_context_t *magma;

    magma = magma_context_self();
    if (magma == NULL) {
        magma_fatal_error("MAGMA_Dealloc_Workspace", "MAGMA not initialized");
        return MAGMA_ERR_NOT_INITIALIZED;
    }
    if (*desc == NULL) {
        magma_error("MAGMA_Dealloc_Workspace", "attempting to deallocate a NULL descriptor");
        return MAGMA_ERR_UNALLOCATED;
    }

    magma_desc_mat_free( *desc );
    free(*desc);
    *desc = NULL;
    return MAGMA_SUCCESS;
}
