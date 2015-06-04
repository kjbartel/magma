/**
 *
 * @file context.c
 *
 *  MAGMA auxiliary routines
 *  MAGMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 1.1.0
 * @author Mathieu Faverge
 * @author Cedric Augonnet
 * @date 2010-11-15
 *
 **/
#include <stdlib.h>
#include "morse_starpu.h"

/***************************************************************************//**
 *  Create new context
 **/
void morse_context_create( magma_context_t *magma )
{
    magma->scheduler = MAGMA_SCHED_STARPU;
    magma->schedopt.starpu = (struct starpu_conf*) malloc (sizeof(struct starpu_conf));

    starpu_conf_init( magma->schedopt.starpu );
    magma->schedopt.starpu->nopencl = 0;

    /* By default, use the dmda strategy */
    if (!getenv("STARPU_SCHED"))
        magma->schedopt.starpu->sched_policy_name = "dmda";

    /* By default, enable calibration */
    if (!getenv("STARPU_CALIBRATE"))
        magma->schedopt.starpu->calibrate = 1;

    return;
}

/***************************************************************************//**
 *  Clean the context
 **/

void morse_context_destroy( magma_context_t *magma )
{
    free(magma->schedopt.starpu);
    return;
}

/***************************************************************************//**
 *
 */
void morse_enable( MAGMA_enum lever )
{
    switch (lever)
    {
        case MAGMA_PROFILING_MODE:
            starpu_profiling_status_set(STARPU_PROFILING_ENABLE);
            break;
        default:
            return;
    }
    return;
}

/***************************************************************************//**
 *
 **/
void morse_disable( MAGMA_enum lever )
{
    switch (lever)
    {
        case MAGMA_PROFILING_MODE:
            starpu_profiling_status_set(STARPU_PROFILING_DISABLE);
            break;
        default:
            return;
    }
    return;
}
