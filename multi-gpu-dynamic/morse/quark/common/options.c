/**
 *
 * @file options.c
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

void morse_options_init( MorseOption_t *options, magma_context_t *magma, 
                         magma_sequence_t *sequence, magma_request_t *request )
{
    /* Create the task flag */
    Quark_Task_Flags *task_flags = (Quark_Task_Flags*) malloc(sizeof(Quark_Task_Flags));
    
    /* Initialize task_flags */
    memset( task_flags, 0, sizeof(Quark_Task_Flags) );
    task_flags->task_lock_to_thread = -1;
    task_flags->task_thread_count   =  1;
    task_flags->thread_set_to_manual_scheduling = -1;
    
    /* Initialize options */
    options->sequence   = sequence;
    options->request    = request;
    options->profiling  = MAGMA_PROFILING == MAGMA_TRUE;
    options->parallel   = MAGMA_PARALLEL == MAGMA_TRUE;
    options->priority   = MORSE_PRIORITY_MIN;

    /* quark in options */
    options->quark      = magma->schedopt.quark;
    options->nb         = MAGMA_NB;
    
    options->ws_hsize   = 0;
    options->ws_dsize   = 0;
    options->ws_host    = NULL;
    options->ws_device  = NULL;

    options->task_flags = task_flags;
    
    QUARK_Task_Flag_Set(task_flags, TASK_SEQUENCE, (intptr_t)sequence->schedopt.quark_sequence);

    return;
}

void morse_options_finalize( MorseOption_t *options, magma_context_t *magma )
{
    /* we can free the task_flags without waiting for quark 
       because they should have been copied for every task */
    free( options->task_flags );
    return;
}

int morse_options_ws_alloc( MorseOption_t *options, size_t hsize, size_t dsize )
{
    options->ws_hsize = hsize; 
    options->ws_dsize = dsize; 
    return MAGMA_SUCCESS;
}

int morse_options_ws_free( MorseOption_t *options )
{
    options->ws_hsize = 0; 
    options->ws_dsize = 0; 
    return MAGMA_SUCCESS;
}
