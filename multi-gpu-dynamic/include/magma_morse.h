/**
 *
 * @file magma_starpu.h
 *
 *  MAGMA codelets kernel
 *  MAGMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver,
 *  and INRIA Bordeaux Sud-Ouest
 *
 * @version 2.3.1
 * @author Cedric Augonnet
 * @author Mathieu Faverge
 * @date 2011-06-01
 *
 **/

/***************************************************************************//**
 *  MAGMA facilities of interest to both src and magmablas directories
 **/
#ifndef _MAGMA_MORSE_H_
#define _MAGMA_MORSE_H_

#if defined(MORSE_USE_CUDA)
#include <magma.h>
#endif
#include <plasma.h>

typedef int MAGMA_enum;
typedef int MAGMA_bool;

#if !defined(__STARPU_DATA_H__)
struct starpu_data_state_t;
typedef struct starpu_data_state_t * starpu_data_handle;
#endif

typedef struct magma_desc_s {
    PLASMA_desc desc;
    int occurences;
    union {
        starpu_data_handle *starpu_handles;
    } schedopt;
} magma_desc_t;

struct magma_context_s;
typedef struct magma_context_s magma_context_t;

/** ****************************************************************************
 *  MAGMA request uniquely identifies each asynchronous function call.
 **/
typedef struct magma_request_s {
    MAGMA_enum status; // MAGMA_SUCCESS or appropriate error code
} magma_request_t;

#define MAGMA_REQUEST_INITIALIZER {MAGMA_SUCCESS}

/** ****************************************************************************
 *  MAGMA sequence uniquely identifies a set of asynchronous function calls
 *  sharing common exception handling.
 **/

typedef struct magma_sequence_s {
    MAGMA_bool       status;    /* MAGMA_SUCCESS or appropriate error code */
    magma_request_t *request;   /* failed request                          */
    union {
        Quark_Sequence *quark_sequence; /* QUARK sequence associated with MAGMA sequence */
    } schedopt;
} magma_sequence_t;

#include "morse.h"

/** ****************************************************************************
 *  MAGMA constants - boolean
 **/
#define MAGMA_FALSE 0
#define MAGMA_TRUE  1

#define MAGMA_CPU        ((1ULL)<<1)
#define MAGMA_CUDA        ((1ULL)<<3)

/** ****************************************************************************
 *  State machine switches
 **/
#define MAGMA_WARNINGS        1
#define MAGMA_ERRORS          2
#define MAGMA_AUTOTUNING      3
#define MAGMA_PROFILING_MODE  4
#define MAGMA_PARALLEL_MODE   5

/** ****************************************************************************
 *  MAGMA constants - configuration parameters
 **/
#define MAGMA_CONCURRENCY      1
#define MAGMA_TILE_SIZE        2
#define MAGMA_INNER_BLOCK_SIZE 3
#define MAGMA_HOUSEHOLDER_MODE 5
#define MAGMA_HOUSEHOLDER_SIZE 6
#define MAGMA_TRANSLATION_MODE 7

#define MAGMA_FLAT_HOUSEHOLDER 1
#define MAGMA_TREE_HOUSEHOLDER 2

#define MAGMA_INPLACE    1
#define MAGMA_OUTOFPLACE 2

/** ****************************************************************************
 *  MAGMA constants - success & error codes
 **/
#define MAGMA_SUCCESS                 0
#define MAGMA_ERR_NOT_INITIALIZED  -101
#define MAGMA_ERR_REINITIALIZED    -102
#define MAGMA_ERR_NOT_SUPPORTED    -103
#define MAGMA_ERR_ILLEGAL_VALUE    -104
#define MAGMA_ERR_NOT_FOUND        -105
#define MAGMA_ERR_OUT_OF_RESOURCES -106
#define MAGMA_ERR_INTERNAL_LIMIT   -107
#define MAGMA_ERR_UNALLOCATED      -108
#define MAGMA_ERR_FILESYSTEM       -109
#define MAGMA_ERR_UNEXPECTED       -110
#define MAGMA_ERR_SEQUENCE_FLUSHED -111

int MAGMA_Init(int nworkers, int ncudas);
int MAGMA_InitPar(int nworkers, int ncudas, int nthreads_per_worker);
int MAGMA_Finalize(void);

/* Descriptor */
int MAGMA_Desc_Create(magma_desc_t **desc, void *mat, MAGMA_enum dtyp, 
                      int mb, int nb, int bsiz, int lm, int ln, 
                      int i, int j, int m, int n);
int MAGMA_Desc_Destroy(magma_desc_t **desc);
int MAGMA_Desc_acquire(magma_desc_t  *desc);
int MAGMA_Desc_release(magma_desc_t  *desc);

/*int MAGMA_Dealloc_Handle_Tile(magma_desc_t **desc);*/

int MAGMA_Enable (MAGMA_enum option);
int MAGMA_Disable(MAGMA_enum option);
int MAGMA_Set(MAGMA_enum param, int  value);
int MAGMA_Get(MAGMA_enum param, int *value);

/*void MAGMA_profile_display(void);*/

#include "magma_morse_z.h"
#include "magma_morse_c.h"
#include "magma_morse_d.h"
#include "magma_morse_s.h"

#endif /* __MAGMA_MORSE_H__ */
