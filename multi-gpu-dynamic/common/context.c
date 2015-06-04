/**
 *
 * @file context.c
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
#include <stdlib.h>
#include "common.h"
#include "context.h"

/***************************************************************************//**
 *  Global data
 **/
/* master threads context lookup table */
static magma_context_t *magma_ctxt = NULL;

/***************************************************************************//**
 *  Create new context
 **/
magma_context_t *magma_context_create()
{
    magma_context_t *magma;

    if ( magma_ctxt != NULL ) {
        magma_error("magma_context_create", "a context is already existing\n");
        return NULL;
    }

    magma = (magma_context_t*)malloc(sizeof(magma_context_t));
    if (magma == NULL) {
        magma_error("magma_context_create", "malloc() failed");
        return NULL;
    }

    magma->scheduler          = MAGMA_SCHED_QUARK;
    magma->nworkers           = 1;
    magma->ncudas             = 0;
    magma->nthreads_per_worker= 1;

    magma->errors_enabled     = MAGMA_FALSE;
    magma->warnings_enabled   = MAGMA_FALSE;
    magma->autotuning_enabled = MAGMA_TRUE;
    magma->parallel_enabled   = MAGMA_FALSE;
    magma->profiling_enabled  = MAGMA_FALSE;

    magma->householder = MAGMA_FLAT_HOUSEHOLDER;
    magma->translation = MAGMA_OUTOFPLACE;

    /* These initializations are just in case the user
       disables autotuning and does not set nb and ib */
    magma->nb = 128;
    magma->ib = 32;
    magma->rhblock = 4;

    /* Initialize scheduler */
    morse_context_create(magma);

    magma_ctxt = magma;
    return magma;
}


/***************************************************************************//**
 *  Return context for a thread
 **/
magma_context_t *magma_context_self()
{
    return magma_ctxt;
}

/***************************************************************************//**
 *  Clean the context
 **/

int magma_context_destroy(){

    morse_context_destroy(magma_ctxt);
    free(magma_ctxt);
    magma_ctxt = NULL;
    
    return MAGMA_SUCCESS;
}

/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  MAGMA_Enable - Enable MAGMA feature.
 *
 *******************************************************************************
 *
 * @param[in] lever
 *          Feature to be enabled:
 *          @arg MAGMA_WARNINGS   printing of warning messages,
 *          @arg MAGMA_ERRORS     printing of error messages,
 *          @arg MAGMA_AUTOTUNING autotuning for tile size and inner block size.
 *          @arg MAGMA_PROFILING_MODE  activate profiling of kernels
 *
 *******************************************************************************
 *
 * @return
 *          \retval MAGMA_SUCCESS successful exit
 *
 ******************************************************************************/
int MAGMA_Enable(MAGMA_enum option)
{
    magma_context_t *magma;

    magma = magma_context_self();
    if (magma == NULL) {
        magma_error("MAGMA_Enable", "MAGMA not initialized");
        return MAGMA_ERR_NOT_INITIALIZED;
    }

    switch (option)
    {
        case MAGMA_WARNINGS:
            magma->warnings_enabled = MAGMA_TRUE;
            break;
        case MAGMA_ERRORS:
            magma->errors_enabled = MAGMA_TRUE;
            break;
        case MAGMA_AUTOTUNING:
            magma->autotuning_enabled = MAGMA_TRUE;
            break;
        case MAGMA_PROFILING_MODE:
            magma->profiling_enabled = MAGMA_TRUE;
            break;
        /* case MAGMA_PARALLEL: */
        /*     magma->parallel_enabled = MAGMA_TRUE; */
        /*     break; */
        default:
            magma_error("MAGMA_Enable", "illegal parameter value");
            return MAGMA_ERR_ILLEGAL_VALUE;
    }
    
    /* Enable at the lower level if required */
    morse_enable( option );

    return MAGMA_SUCCESS;
}

/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  MAGMA_Disable - Disable MAGMA feature.
 *
 *******************************************************************************
 *
 * @param[in] lever
 *          Feature to be disabled:
 *          @arg MAGMA_WARNINGS   printing of warning messages,
 *          @arg MAGMA_ERRORS     printing of error messages,
 *          @arg MAGMA_AUTOTUNING autotuning for tile size and inner block size.
 *          @arg MAGMA_PROFILING_MODE  activate profiling of kernels
 *
 *******************************************************************************
 *
 * @return
 *          \retval MAGMA_SUCCESS successful exit
 *
 ******************************************************************************/
int MAGMA_Disable(MAGMA_enum option)
{
    magma_context_t *magma;

    magma = magma_context_self();
    if (magma == NULL) {
        magma_error("MAGMA_Disable", "MAGMA not initialized");
        return MAGMA_ERR_NOT_INITIALIZED;
    }
    switch ( option )
    {
        case MAGMA_WARNINGS:
            magma->warnings_enabled = MAGMA_FALSE;
            break;
        case MAGMA_ERRORS:
            magma->errors_enabled = MAGMA_FALSE;
            break;
        case MAGMA_AUTOTUNING:
            magma->autotuning_enabled = MAGMA_FALSE;
            break;
        case MAGMA_PROFILING_MODE:
            magma->profiling_enabled = MAGMA_FALSE;
            break;
        case MAGMA_PARALLEL_MODE:
            magma->parallel_enabled = MAGMA_FALSE;
            break;
        default:
            magma_error("MAGMA_Disable", "illegal parameter value");
            return MAGMA_ERR_ILLEGAL_VALUE;
    }
    
    /* Disable at the lower level if required */
    morse_disable( option );

    return MAGMA_SUCCESS;
}

/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  MAGMA_Set - Set MAGMA parameter.
 *
 *******************************************************************************
 *
 * @param[in] param
 *          Feature to be enabled:
 *          @arg MAGMA_TILE_SIZE:        size matrix tile,
 *          @arg MAGMA_INNER_BLOCK_SIZE: size of tile inner block,
 *
 * @param[in] value
 *          Value of the parameter.
 *
 *******************************************************************************
 *
 * @return
 *          \retval MAGMA_SUCCESS successful exit
 *
 ******************************************************************************/
int MAGMA_Set(MAGMA_enum param, int value)
{
    magma_context_t *magma;

    magma = magma_context_self();
    if (magma == NULL) {
        magma_error("MAGMA_Set", "MAGMA not initialized");
        return MAGMA_ERR_NOT_INITIALIZED;
    }
    switch (param) {
        case MAGMA_TILE_SIZE:
            if (value <= 0) {
                magma_error("MAGMA_Set", "negative tile size");
                return MAGMA_ERR_ILLEGAL_VALUE;
            }
            magma->nb = value;
            if ( magma->autotuning_enabled ) {
                magma->autotuning_enabled = MAGMA_FALSE;
                magma_warning("MAGMA_Set", "autotuning has been automatically disable\n");
            }
            /* Limit ib to nb */
            magma->ib = min( magma->nb, magma->ib );
            break;
        case MAGMA_INNER_BLOCK_SIZE:
            if (value <= 0) {
                magma_error("MAGMA_Set", "negative inner block size");
                return MAGMA_ERR_ILLEGAL_VALUE;
            }
            if (value > magma->nb) {
                magma_error("MAGMA_Set", "inner block larger than tile");
                return MAGMA_ERR_ILLEGAL_VALUE;
            }
            /* if (magma->nb % value != 0) { */
            /*     magma_error("MAGMA_Set", "inner block does not divide tile"); */
            /*     return MAGMA_ERR_ILLEGAL_VALUE; */
            /* } */
            magma->ib = value;

            if ( magma->autotuning_enabled ) {
                magma->autotuning_enabled = MAGMA_FALSE;
                magma_warning("MAGMA_Set", "autotuning has been automatically disable\n");
            }
            break;
        case MAGMA_HOUSEHOLDER_MODE:
            if (value != MAGMA_FLAT_HOUSEHOLDER && value != MAGMA_TREE_HOUSEHOLDER) {
                magma_error("MAGMA_Set", "illegal value of MAGMA_HOUSEHOLDER_MODE");
                return MAGMA_ERR_ILLEGAL_VALUE;
            }
            magma->householder = value;
            break;
        case MAGMA_HOUSEHOLDER_SIZE:
            if (value <= 0) {
                magma_error("MAGMA_Set", "negative householder size");
                return MAGMA_ERR_ILLEGAL_VALUE;
            }
            magma->rhblock = value;
            break;
        case MAGMA_TRANSLATION_MODE:
            if (value != MAGMA_INPLACE && value != MAGMA_OUTOFPLACE) {
                magma_error("MAGMA_Set", "illegal value of MAGMA_TRANSLATION_MODE");
                return MAGMA_ERR_ILLEGAL_VALUE;
            }
            magma->translation = value;
            break;
        default:
            magma_error("MAGMA_Set", "unknown parameter");
            return MAGMA_ERR_ILLEGAL_VALUE;
    }

    return MAGMA_SUCCESS;
}

/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  MAGMA_Get - Get value of MAGMA parameter.
 *
 *******************************************************************************
 *
 * @param[in] param
 *          Feature to be enabled:
 *          @arg MAGMA_TILE_SIZE:        size matrix tile,
 *          @arg MAGMA_INNER_BLOCK_SIZE: size of tile inner block,
 *
 * @param[out] value
 *          Value of the parameter.
 *
 *******************************************************************************
 *
 * @return
 *          \retval MAGMA_SUCCESS successful exit
 *
 ******************************************************************************/
int MAGMA_Get(MAGMA_enum param, int *value)
{
    magma_context_t *magma;

    magma = magma_context_self();
    if (magma == NULL) {
        magma_error("MAGMA_Get", "MAGMA not initialized");
        return MAGMA_ERR_NOT_INITIALIZED;
    }
    switch (param) {
        case MAGMA_TILE_SIZE:
            *value = magma->nb;
            return MAGMA_SUCCESS;
        case MAGMA_INNER_BLOCK_SIZE:
            *value = magma->ib;
            return MAGMA_SUCCESS;
        case MAGMA_HOUSEHOLDER_MODE:
            *value = magma->householder;
            return MAGMA_SUCCESS;
        case MAGMA_HOUSEHOLDER_SIZE:
            *value = magma->rhblock;
            return MAGMA_SUCCESS;
        case MAGMA_TRANSLATION_MODE:
            *value = magma->translation;
            return MAGMA_SUCCESS;
        default:
            magma_error("MAGMA_Get", "unknown parameter");
            return MAGMA_ERR_ILLEGAL_VALUE;
    }

    return MAGMA_SUCCESS;
}
