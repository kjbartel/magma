/**
 *
 * @file context.c
 *
 *  MAGMA auxiliary routines
 *  MAGMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 
 * @author Vijay Joshi
 * @date 2011-10-29 
 *
 **/
#include <stdlib.h>
#include "morse_quark.h"

/***************************************************************************//**
 *  Create new context
 **/
void morse_context_create(magma_context_t *magma)
{
    magma->scheduler = MAGMA_SCHED_QUARK;
    /* Will require the static initialization if we want to use it in this code */
    return;
}

/***************************************************************************//**
 *  Clean the context
 **/

void morse_context_destroy(magma_context_t *magma)
{
    return;
}

/***************************************************************************//**
 *
 */
void morse_enable( magma_context_t *magma, MAGMA_enum lever )
{
    switch (lever)
    {
        case MAGMA_PROFILING_MODE:
            fprintf(stderr, "Profiling is not available with Quark\n");
            break;
        default:
            return;
    }
    return;
}

/***************************************************************************//**
 *
 **/
void morse_disable( magma_context_t *magma, MAGMA_enum lever )
{
    switch (lever)
    {
        case MAGMA_PROFILING_MODE:
            fprintf(stderr, "Profiling is not available with Quark\n");
            break;
        default:
            return;
    }
    return;
}
