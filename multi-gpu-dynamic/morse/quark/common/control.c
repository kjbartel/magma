/**
 *
 * @file control.c
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 
 * @author Vijay Joshi
 * @date 2011-10-29
 *
 **/
#include <stdio.h>
#include <stdlib.h>
#include "morse_quark.h"

/***************************************************************************//**
 *  Busy-waiting barrier
 **/
void morse_barrier( magma_context_t *magma )
{
    QUARK_Barrier(magma->schedopt.quark);
}

int morse_init_scheduler( magma_context_t *magma, int nworkers, int ncudas, int nthreads_per_worker)
{
    hres = 0;
    if ( ncudas > 0 )
        magma_warning( "morse_init_scheduler(quark)", "GPUs are not supported for now");

    if ( nthreads_per_worker > 0 )
        magma_warning( "morse_init_scheduler(quark)", "Multi-threaded kernels are not supported for now");

    magma->schedopt.quark = QUARK_New(nworkers); 

    return hres;
}

/***************************************************************************//**
 *
 */
void morse_finalize_scheduler( magma_context_t *magma )
{
    QUARK_Delete(magma->schedopt.quark);
    return;
}

