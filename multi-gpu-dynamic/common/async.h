/***
 * 
 *
 * @file async.h
 *
 *  MAGMA auxiliary routines
 *  MAGMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.3.1
 * @author Jakub Kurzak
 * @date 2010-11-15
 *
 **/
#ifndef _MAGMA_ASYNC_H_
#define _MAGMA_ASYNC_H_

#ifdef __cplusplus
extern "C" {
#endif

/***************************************************************************//**
 *  Internal routines
 **/
int  magma_request_fail(magma_sequence_t *sequence, magma_request_t *request, int error);
int  magma_sequence_create(magma_context_t *magma, magma_sequence_t **sequence);
int  magma_sequence_destroy(magma_context_t *magma, magma_sequence_t *sequence);
int  magma_sequence_wait(magma_context_t *magma, magma_sequence_t *sequence);
void magma_sequence_flush(void *quark, magma_sequence_t *sequence, magma_request_t *request, int status);

/***************************************************************************//**
 *  User routines
 **/
int MAGMA_Sequence_Create(magma_sequence_t **sequence);
int MAGMA_Sequence_Destroy(magma_sequence_t *sequence);
int MAGMA_Sequence_Wait(magma_sequence_t *sequence);
int MAGMA_Sequence_Flush(magma_sequence_t *sequence, magma_request_t *request);

#ifdef __cplusplus
}
#endif

#endif
