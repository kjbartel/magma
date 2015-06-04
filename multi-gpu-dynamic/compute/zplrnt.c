/**
 *
 *  @file zplrnt.c
 *
 *  MAGMA compute
 *  MAGMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver,
 *  and INRIA Bordeaux Sud-Ouest
 *
 *  @version 2.3.1
 *  @author Mathieu Faverge
 *  @date 2011-06-01
 *  @precisions normal z -> c d s
 *
 **/
#include "common.h"

/***************************************************************************//**
 *
 * @ingroup MAGMA_Complex64_t
 *
 *  MAGMA_zplrnt - Generate a random matrix by tiles.
 *
 *******************************************************************************
 *
 * @param[in] M
 *          The number of rows of A.
 *
 * @param[in] N
 *          The order of the matrix A. N >= 0.
 *
 * @param[out] A
 *          On exit, The random matrix A generated.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,M).
 *
 * @param[in] seed
 *          The seed used in the random generation.
 *
 *******************************************************************************
 *
 * @return
 *          \retval MAGMA_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *
 *******************************************************************************
 *
 * @sa MAGMA_zplrnt_Tile
 * @sa MAGMA_zplrnt_Tile_Async
 * @sa MAGMA_cplrnt
 * @sa MAGMA_dplrnt
 * @sa MAGMA_splrnt
 * @sa MAGMA_zplghe
 * @sa MAGMA_zplgsy
 *
 ******************************************************************************/
int MAGMA_zplrnt( int M, int N,
                   PLASMA_Complex64_t *A, int LDA,
                   unsigned long long int seed )
{
    int NB;
    int status;
    magma_context_t *magma;
    magma_sequence_t *sequence = NULL;
    magma_request_t request = MAGMA_REQUEST_INITIALIZER;
    magma_desc_t descA;

    magma = magma_context_self();
    if (magma == NULL) {
        magma_fatal_error("MAGMA_zplrnt", "MAGMA not initialized");
        return MAGMA_ERR_NOT_INITIALIZED;
    }
    /* Check input arguments */
    if (M < 0) {
        magma_error("MAGMA_zplrnt", "illegal value of M");
        return -1;
    }
    if (N < 0) {
        magma_error("MAGMA_zplrnt", "illegal value of N");
        return -2;
    }
    if (LDA < max(1, M)) {
        magma_error("MAGMA_zplrnt", "illegal value of LDA");
        return -4;
    }
    /* Quick return */
    if (min(M, N) == 0)
        return MAGMA_SUCCESS;

    /* Tune NB depending on M, N && NRHS; Set NBNB */
    /* status = magma_tune(MAGMA_FUNC_ZGEMM, M, N, 0); */
    /* if (status != MAGMA_SUCCESS) { */
    /*     magma_error("MAGMA_zplrnt", "magma_tune() failed"); */
    /*     return status; */
    /* } */
    
    /* Set NT */
    NB = MAGMA_NB;
    magma_sequence_create(magma, &sequence);

    magma_zdesc_alloc( descA, NB, NB, LDA, N, 0, 0, N, N, magma_desc_mat_free(&(descA)) );

    /* Call the tile interface */
    MAGMA_zplrnt_Tile_Async( &descA, seed, sequence, &request );

    morse_barrier( magma );
    
    magma_zooptile2lap( descA, A, NB, NB, LDA, N );
    morse_barrier( magma );

    status = sequence->status;
    magma_sequence_destroy(magma, sequence);

    return status;
}

/***************************************************************************//**
 *
 * @ingroup MAGMA_Complex64_t_Tile
 *
 *  MAGMA_zplrnt_Tile - Generate a random matrix by tiles.
 *  Tile equivalent of MAGMA_zplrnt().
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] A
 *          On exit, The random matrix A generated.
 *
 * @param[in] seed
 *          The seed used in the random generation.
 *
 *******************************************************************************
 *
 * @return
 *          \retval MAGMA_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa MAGMA_zplrnt
 * @sa MAGMA_zplrnt_Tile_Async
 * @sa MAGMA_cplrnt_Tile
 * @sa MAGMA_dplrnt_Tile
 * @sa MAGMA_splrnt_Tile
 * @sa MAGMA_zplghe_Tile
 * @sa MAGMA_zplgsy_Tile
 *
 ******************************************************************************/
int MAGMA_zplrnt_Tile( magma_desc_t *A,
                       unsigned long long int seed )
{
    magma_context_t *magma;
    magma_sequence_t *sequence = NULL;
    magma_request_t request = MAGMA_REQUEST_INITIALIZER;
    int status;

    magma = magma_context_self();
    if (magma == NULL) {
        magma_fatal_error("MAGMA_zplrnt_Tile", "MAGMA not initialized");
        return MAGMA_ERR_NOT_INITIALIZED;
    }
    magma_sequence_create(magma, &sequence);
    MAGMA_zplrnt_Tile_Async( A, seed, sequence, &request );
    morse_barrier( magma );
    morse_desc_getoncpu( A );
    status = sequence->status;
    magma_sequence_destroy(magma, sequence);
    return status;
}

/***************************************************************************//**
 *
 * @ingroup MAGMA_Complex64_t_Tile_Async
 *
 *  MAGMA_zplrnt_Tile_Async - Generate a random matrix by tiles.
 *  Non-blocking equivalent of MAGMA_zplrnt_Tile().
 *  May return before the computation is finished.
 *  Allows for pipelining of operations ar runtime.
 *
 *******************************************************************************
 *
 * @param[in] sequence
 *          Identifies the sequence of function calls that this call belongs to
 *          (for completion checks and exception handling purposes).
 *
 * @param[out] request
 *          Identifies this function call (for exception handling purposes).
 *
 *******************************************************************************
 *
 * @sa MAGMA_zplrnt
 * @sa MAGMA_zplrnt_Tile
 * @sa MAGMA_cplrnt_Tile_Async
 * @sa MAGMA_dplrnt_Tile_Async
 * @sa MAGMA_splrnt_Tile_Async
 * @sa MAGMA_zplghe_Tile_Async
 * @sa MAGMA_zplgsy_Tile_Async
 *
 ******************************************************************************/
int MAGMA_zplrnt_Tile_Async( magma_desc_t     *A,
                             unsigned long long int seed,
                             magma_sequence_t *sequence, 
                             magma_request_t  *request)
{
    magma_context_t *magma;

    magma = magma_context_self();
    if (magma == NULL) {
        magma_fatal_error("MAGMA_zplrnt_Tile", "MAGMA not initialized");
        return MAGMA_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        magma_fatal_error("MAGMA_zplrnt_Tile", "NULL sequence");
        return MAGMA_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        magma_fatal_error("MAGMA_zplrnt_Tile", "NULL request");
        return MAGMA_ERR_UNALLOCATED;
    }
    /* Check sequence status */
    if (sequence->status == MAGMA_SUCCESS)
        request->status = MAGMA_SUCCESS;
    else
        return magma_request_fail(sequence, request, MAGMA_ERR_SEQUENCE_FLUSHED);

    /* Check descriptors for correctness */
    if (magma_desc_check( A ) != MAGMA_SUCCESS) {
        magma_error("MAGMA_zplrnt_Tile", "invalid descriptor");
        return magma_request_fail(sequence, request, MAGMA_ERR_ILLEGAL_VALUE);
    }
    /* Check input arguments */
    if (A->nb != A->mb) {
        magma_error("MAGMA_zplrnt_Tile", "only square tiles supported");
        return magma_request_fail(sequence, request, MAGMA_ERR_ILLEGAL_VALUE);
    }

    /* Quick return */
    if (min( A->m, A->n ) == 0)
        return MAGMA_SUCCESS;

    magma_pzplrnt( A, seed, sequence, request);

    return PLASMA_SUCCESS;
}
