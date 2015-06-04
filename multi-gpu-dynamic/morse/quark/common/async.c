/**
 *
 * @file async.c
 *
 *  MAGMA auxiliary routines
 *  MAGMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 
 * @author Jakub Kurzak
 * @author Vijay Joshi
 * @date 2011-10-29
 *
 **/
#include <stdlib.h>
#include "morse_quark.h"

/***************************************************************************//**
 *  Create a sequence
 **/
int morse_sequence_create(magma_context_t *magma, magma_sequence_t *sequence)
{
    if((sequence->schedopt.quark_sequence = QUARK_Sequence_Create(magma->schedopt.quark)) == NULL){
        magma_error("MAGMA_Sequence_Create", "QUARK_Sequence_Create() failed");
        return MAGMA_ERR_OUT_OF_RESOURCES;
    }
    sequence->status = MAGMA_SUCCESS;
    return MAGMA_SUCCESS;
}

/***************************************************************************//**
 *  Destroy a sequence
 **/
int morse_sequence_destroy(magma_context_t *magma, magma_sequence_t *sequence)
{
    QUARK_Sequence_Destroy(magma->schedopt.quark, sequence->schedopt.quark_sequence);
    return MAGMA_SUCCESS;
}

/***************************************************************************//**
 *  Wait for the completion of a sequence
 **/
int morse_sequence_wait( magma_context_t *magma, magma_sequence_t *sequence )
{
    QUARK_Sequence_Wait(magma->schedopt.quark, sequence->schedopt.quark_sequence);
    return MAGMA_SUCCESS;
}

/***************************************************************************//**
 *  Terminate a sequence
 **/
void morse_sequence_flush(void *quark, magma_sequence_t *sequence, magma_request_t *request, int status)
{
    sequence->request = request;
    sequence->status = status;
    request->status = status;
    QUARK_Sequence_Cancel( (Quark *)quark, sequence->schedopt.quark_sequence);
}

