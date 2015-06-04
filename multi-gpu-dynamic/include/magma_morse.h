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
typedef struct starpu_data_state_t * starpu_data_handle_t;
#endif

struct magma_desc_s;
typedef struct magma_desc_s  magma_desc_t;

struct magma_desc_s {
    void *(*get_blkaddr)( const magma_desc_t*, int, int );
    int   (*get_blkldd )( const magma_desc_t*, int );
    void *mat;          // pointer to the beginning of the matrix
    size_t A21;        // pointer to the beginning of the matrix A21
    size_t A12;        // pointer to the beginning of the matrix A12
    size_t A22;        // pointer to the beginning of the matrix A22
    PLASMA_enum styp;   // storage layout of the matrix
    PLASMA_enum dtyp;   // precision of the matrix
    int mb;             // number of rows in a tile
    int nb;             // number of columns in a tile
    int bsiz;           // size in elements including padding
    int lm;             // number of rows of the entire matrix
    int ln;             // number of columns of the entire matrix
    int lm1;            // number of tile rows of the A11 matrix - derived parameter
    int ln1;            // number of tile columns of the A11 matrix - derived parameter
    int lmt;            // number of tile rows of the entire matrix - derived parameter
    int lnt;            // number of tile columns of the entire matrix - derived parameter
    int i;              // row index to the beginning of the submatrix
    int j;              // column index to the beginning of the submatrix
    int m;              // number of rows of the submatrix
    int n;              // number of columns of the submatrix
    int mt;             // number of tile rows of the submatrix - derived parameter
    int nt;             // number of tile columns of the submatrix - derived parameter
    int occurences;
    int myrank;
    union {
        starpu_data_handle_t *starpu_handles;
    } schedopt;
};

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
#ifndef _MAGMA_
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
#endif

int MAGMA_Init(int nworkers, int ncudas);
int MAGMA_InitPar(int nworkers, int ncudas, int nthreads_per_worker);
int MAGMA_Finalize(void);
int MAGMA_my_mpi_rank(void);

/* Descriptor */
int MAGMA_Desc_Create(magma_desc_t **desc, void *mat, MAGMA_enum dtyp, 
                      int mb, int nb, int bsiz, int lm, int ln, 
                      int i, int j, int m, int n);
int MAGMA_Desc_Destroy(magma_desc_t **desc);
int MAGMA_Desc_acquire(magma_desc_t  *desc);
int MAGMA_Desc_release(magma_desc_t  *desc);

/* Workspaces */
int MAGMA_Dealloc_Workspace(magma_desc_t **desc);

/* Options */
int MAGMA_Enable (MAGMA_enum option);
int MAGMA_Disable(MAGMA_enum option);
int MAGMA_Set(MAGMA_enum param, int  value);
int MAGMA_Get(MAGMA_enum param, int *value);

/* Sequences */
int MAGMA_Sequence_Create (magma_sequence_t **sequence);
int MAGMA_Sequence_Destroy(magma_sequence_t *sequence);
int MAGMA_Sequence_Wait   (magma_sequence_t *sequence);
int MAGMA_Sequence_Flush  (magma_sequence_t *sequence, magma_request_t *request);

/*void MAGMA_profile_display(void);*/

#include "magma_morse_z.h"
#include "magma_morse_c.h"
#include "magma_morse_d.h"
#include "magma_morse_s.h"

#endif /* __MAGMA_MORSE_H__ */
