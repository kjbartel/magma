/**
 *
 * @file async.c
 *
 *  MAGMA auxiliary routines
 *  MAGMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 1.1.0
 * @author Jakub Kurzak
 * @author Matheu Faverge
 * @date 2010-11-15
 *
 **/
#include <stdlib.h>
#include "common.h"

/***************************************************************************//**
 *  Register an exception.
 **/
int magma_request_fail(magma_sequence_t *sequence, magma_request_t *request, int status)
{
    sequence->request = request;
    sequence->status = status;
    request->status = status;
    return status;
}

/***************************************************************************//**
 *  Create a sequence
 **/
int magma_sequence_create(magma_context_t *magma, magma_sequence_t **sequence)
{
    if ((*sequence = malloc(sizeof(magma_sequence_t))) == NULL) {
        magma_error("MAGMA_Sequence_Create", "malloc() failed");
        return MAGMA_ERR_OUT_OF_RESOURCES;
    }
    
    morse_sequence_create( magma, *sequence );

    (*sequence)->status = MAGMA_SUCCESS;
    return MAGMA_SUCCESS;
}

/***************************************************************************//**
 *  Destroy a sequence
 **/
int magma_sequence_destroy(magma_context_t *magma, magma_sequence_t *sequence)
{
    morse_sequence_destroy( magma, sequence );
    free(sequence);
    return MAGMA_SUCCESS;
}

/***************************************************************************//**
 *  Wait for the completion of a sequence
 **/
int magma_sequence_wait(magma_context_t *magma, magma_sequence_t *sequence)
{
    morse_sequence_wait( magma, sequence );
    return MAGMA_SUCCESS;
}

/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  MAGMA_Sequence_Create - Create a squence.
 *
 *******************************************************************************
 *
 * @param[out] sequence
 *          Identifies a set of routines sharing common exception handling.
 *
 *******************************************************************************
 *
 * @return
 *          \retval MAGMA_SUCCESS successful exit
 *
 ******************************************************************************/
int MAGMA_Sequence_Create(magma_sequence_t **sequence)
{
    magma_context_t *magma;
    int status;

    magma = magma_context_self();
    if (magma == NULL) {
        magma_fatal_error("MAGMA_Sequence_Create", "MAGMA not initialized");
        return MAGMA_ERR_NOT_INITIALIZED;
    }
    status = magma_sequence_create(magma, sequence);
    return status;
}

/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  MAGMA_Sequence_Destroy - Destroy a sequence.
 *
 *******************************************************************************
 *
 * @param[in] sequence
 *          Identifies a set of routines sharing common exception handling.
 *
 *******************************************************************************
 *
 * @return
 *          \retval MAGMA_SUCCESS successful exit
 *
 ******************************************************************************/
int MAGMA_Sequence_Destroy(magma_sequence_t *sequence)
{
    magma_context_t *magma;
    int status;

    magma = magma_context_self();
    if (magma == NULL) {
        magma_fatal_error("MAGMA_Sequence_Destroy", "MAGMA not initialized");
        return MAGMA_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        magma_fatal_error("MAGMA_Sequence_Destroy", "NULL sequence");
        return MAGMA_ERR_UNALLOCATED;
    }
    status = magma_sequence_destroy(magma, sequence);
    return status;
}

/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  MAGMA_Sequence_Wait - Wait for the completion of a sequence.
 *
 *******************************************************************************
 *
 * @param[in] sequence
 *          Identifies a set of routines sharing common exception handling.
 *
 *******************************************************************************
 *
 * @return
 *          \retval MAGMA_SUCCESS successful exit
 *
 ******************************************************************************/
int MAGMA_Sequence_Wait(magma_sequence_t *sequence)
{
    magma_context_t *magma;
    int status;

    magma = magma_context_self();
    if (magma == NULL) {
        magma_fatal_error("MAGMA_Sequence_Wait", "MAGMA not initialized");
        return MAGMA_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        magma_fatal_error("MAGMA_Sequence_Wait", "NULL sequence");
        return MAGMA_ERR_UNALLOCATED;
    }
    status = magma_sequence_wait(magma, sequence);
    return status;
}

/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  MAGMA_Sequence_Flush - Terminate a sequence.
 *
 *******************************************************************************
 *
 * @param[in] sequence
 *          Identifies a set of routines sharing common exception handling.
 *
 * @param[in] request
 *          The flush request.
 *
 *******************************************************************************
 *
 * @return
 *          \retval MAGMA_SUCCESS successful exit
 *
 ******************************************************************************/
int MAGMA_Sequence_Flush(magma_sequence_t *sequence, magma_request_t *request)
{
    magma_context_t *magma;

    magma = magma_context_self();
    if (magma == NULL) {
        magma_fatal_error("MAGMA_Sequence_Flush", "MAGMA not initialized");
        return MAGMA_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        magma_fatal_error("MAGMA_Sequence_Flush", "NULL sequence");
        return MAGMA_ERR_UNALLOCATED;
    }
    morse_sequence_flush( (void *)magma->schedopt.quark, sequence, request, MAGMA_ERR_SEQUENCE_FLUSHED);
    return MAGMA_SUCCESS;
}
