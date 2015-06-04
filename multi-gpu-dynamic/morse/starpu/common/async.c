/**
 *
 * @file async.c
 *
 *  MAGMA auxiliary routines
 *  MAGMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 1.1.0
 * @author Mathieu Faverge
 * @date 2010-11-15
 *
 **/
#include <stdlib.h>
#include "morse_starpu.h"

/***************************************************************************//**
 *  Create a sequence
 **/
int morse_sequence_create( magma_context_t *magma, magma_sequence_t *sequence )
{
    (void)magma; (void)sequence;
    return MAGMA_SUCCESS;
}

/***************************************************************************//**
 *  Destroy a sequence
 **/
int morse_sequence_destroy( magma_context_t *magma, magma_sequence_t *sequence )
{
    (void)magma; (void)sequence;
    return MAGMA_SUCCESS;
}

/***************************************************************************//**
 *  Wait for the completion of a sequence
 **/
int morse_sequence_wait( magma_context_t *magma, magma_sequence_t *sequence )
{
    (void)magma; (void)sequence;
    starpu_task_wait_for_all();
    return MAGMA_SUCCESS;
}

/***************************************************************************//**
 *  Terminate a sequence
 **/
void morse_sequence_flush( void *schedopt, 
                           magma_sequence_t *sequence, 
                           magma_request_t  *request, int status)
{
    (void)schedopt; 
    sequence->request = request;
    sequence->status = status;
    request->status = status;
    starpu_task_wait_for_all();
    return;
}

