/**
 *
 * @file control.c
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 1.1.0
 * @author Mathieu Faverge
 * @author Cedric Augonnet
 * @date 2010-11-15
 *
 **/
#include <stdio.h>
#include <stdlib.h>
#include "morse_starpu.h"

void morse_options_init( MorseOption_t *option, magma_context_t *magma, 
                         magma_sequence_t *sequence, magma_request_t *request )
{
    option->sequence   = sequence;
    option->request    = request;
    option->profiling  = MAGMA_PROFILING == MAGMA_TRUE;
    option->parallel   = MAGMA_PARALLEL == MAGMA_TRUE;
    option->priority   = MORSE_PRIORITY_MIN;
    option->quark      = NULL;
    option->task_flags = NULL;
    option->nb         = MAGMA_NB;
    option->ws_hsize   = 0; 
    option->ws_dsize   = 0; 
    option->ws_host    = NULL;  
    option->ws_device  = NULL;
    return;
}

void morse_options_finalize( MorseOption_t *option, magma_context_t *magma )
{
    (void)option; (void)magma;
    return;
}

int morse_options_ws_alloc( MorseOption_t *options, size_t hsize, size_t dsize )
{
    int ret = 0;
    if ( hsize > 0 ) { 
        options->ws_hsize = hsize;
        ret = morse_starpu_ws_alloc((morse_starpu_ws_t**)&(options->ws_host), 
                                    hsize, MAGMA_CPU|MAGMA_CUDA, MAGMA_HOST_MEM  );
    }
#if defined(MORSE_USE_CUDA)
    if ( ret == 0 && dsize > 0 ) { 
        options->ws_dsize = dsize;
        ret = morse_starpu_ws_alloc((morse_starpu_ws_t**)&(options->ws_device), 
                                    dsize, MAGMA_CUDA, MAGMA_WORKER_MEM);
    }
#else
    (void)dsize;
#endif
    return ret;
}

int morse_options_ws_free( MorseOption_t *options )
{
    int ret = 0;
    if ( options->ws_host != NULL ) {
        starpu_task_wait_for_all();
        ret = morse_starpu_ws_free( (morse_starpu_ws_t*)(options->ws_host) );
        options->ws_host = NULL;
    }
#if defined(MORSE_USE_CUDA)
    if ( ret == 0 && options->ws_device != NULL ) {
        starpu_task_wait_for_all();
        ret = morse_starpu_ws_free( (morse_starpu_ws_t*)(options->ws_device) );
        options->ws_device = NULL;
    }
#endif
    return ret;
}
