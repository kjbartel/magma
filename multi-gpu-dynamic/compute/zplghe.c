/**
 *
 *  @file zplghe.c
 *
 *  MAGMA compute
 *  MAGMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver,
 *  and INRIA Bordeaux Sud-Ouest
 *
 *  @version 2.3.1
 *  @author Mathieu Faverge
 *  @date 2011-06-01
 *  @precisions normal z -> c
 *
 **/
#include "common.h"

/***************************************************************************//**
 *
 * @ingroup MAGMA_Complex64_t
 *
 *  MAGMA_zplghe - Generate a random hermitian matrix by tiles.
 *
 *******************************************************************************
 *
 * @param[in] bump
 *          The value to add to the diagonal to be sure 
 *          to have a positive definite matrix.
 *
 * @param[in] N
 *          The order of the matrix A. N >= 0.
 *
 * @param[out] A
 *          On exit, The random hermitian matrix A generated.
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
 * @sa MAGMA_zplghe_Tile
 * @sa MAGMA_zplghe_Tile_Async
 * @sa MAGMA_cplghe
 * @sa MAGMA_dplghe
 * @sa MAGMA_splghe
 * @sa MAGMA_zplrnt
 * @sa MAGMA_zplgsy
 *
 ******************************************************************************/
int MAGMA_zplghe( double bump, int N,
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
        magma_fatal_error("MAGMA_zplghe", "MAGMA not initialized");
        return MAGMA_ERR_NOT_INITIALIZED;
    }
    /* Check input arguments */
    if (N < 0) {
        magma_error("MAGMA_zplghe", "illegal value of N");
        return -2;
    }
    if (LDA < max(1, N)) {
        magma_error("MAGMA_zplghe", "illegal value of LDA");
        return -4;
    }
    /* Quick return */
    if (max(0, N) == 0)
        return MAGMA_SUCCESS;

    /* Tune NB depending on M, N && NRHS; Set NBNB */
    /* status = magma_tune(MAGMA_FUNC_ZGEMM, N, N, 0); */
    /* if (status != MAGMA_SUCCESS) { */
    /*     magma_error("MAGMA_zplghe", "magma_tune() failed"); */
    /*     return status; */
    /* } */
    
    /* Set NT */
    NB = MAGMA_NB;
    magma_sequence_create(magma, &sequence);
    
    magma_zdesc_alloc( descA, NB, NB, LDA, N, 0, 0, N, N, magma_desc_mat_free(&(descA)) );

    /* Call the tile interface */
    MAGMA_zplghe_Tile_Async( bump, &descA, seed, sequence, &request );

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
 *  MAGMA_zplghe_Tile - Generate a random hermitian matrix by tiles.
 *  Tile equivalent of MAGMA_zplghe().
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] bump
 *          The value to add to the diagonal to be sure 
 *          to have a positive definite matrix.
 *
 * @param[in] A
 *          On exit, The random hermitian matrix A generated.
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
 * @sa MAGMA_zplghe
 * @sa MAGMA_zplghe_Tile_Async
 * @sa MAGMA_cplghe_Tile
 * @sa MAGMA_dplghe_Tile
 * @sa MAGMA_splghe_Tile
 * @sa MAGMA_zplrnt_Tile
 * @sa MAGMA_zplgsy_Tile
 *
 ******************************************************************************/
int MAGMA_zplghe_Tile( double bump, magma_desc_t *A,
                        unsigned long long int seed )
{
    magma_context_t *magma;
    magma_sequence_t *sequence = NULL;
    magma_request_t request = MAGMA_REQUEST_INITIALIZER;
    int status;

    magma = magma_context_self();
    if (magma == NULL) {
        magma_fatal_error("MAGMA_zplghe_Tile", "MAGMA not initialized");
        return MAGMA_ERR_NOT_INITIALIZED;
    }
    magma_sequence_create(magma, &sequence);
    MAGMA_zplghe_Tile_Async( bump, A, seed, sequence, &request );
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
 *  MAGMA_zplghe_Tile_Async - Generate a random hermitian matrix by tiles.
 *  Non-blocking equivalent of MAGMA_zplghe_Tile().
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
 * @sa MAGMA_zplghe
 * @sa MAGMA_zplghe_Tile
 * @sa MAGMA_cplghe_Tile_Async
 * @sa MAGMA_dplghe_Tile_Async
 * @sa MAGMA_splghe_Tile_Async
 * @sa MAGMA_zplghe_Tile_Async
 * @sa MAGMA_zplgsy_Tile_Async
 *
 ******************************************************************************/
int MAGMA_zplghe_Tile_Async( double          bump,
                             magma_desc_t     *A,
                             unsigned long long int seed,
                             magma_sequence_t *sequence, 
                             magma_request_t  *request)
{
    PLASMA_desc descA = A->desc;
    magma_context_t *magma;

    magma = magma_context_self();
    if (magma == NULL) {
        magma_fatal_error("MAGMA_zplghe_Tile", "MAGMA not initialized");
        return MAGMA_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        magma_fatal_error("MAGMA_zplghe_Tile", "NULL sequence");
        return MAGMA_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        magma_fatal_error("MAGMA_zplghe_Tile", "NULL request");
        return MAGMA_ERR_UNALLOCATED;
    }
    /* Check sequence status */
    if (sequence->status == MAGMA_SUCCESS)
        request->status = MAGMA_SUCCESS;
    else
        return magma_request_fail(sequence, request, MAGMA_ERR_SEQUENCE_FLUSHED);

    /* Check descriptors for correctness */
    if (magma_desc_check( A ) != MAGMA_SUCCESS) {
        magma_error("MAGMA_zplghe_Tile", "invalid descriptor");
        return magma_request_fail(sequence, request, MAGMA_ERR_ILLEGAL_VALUE);
    }
    /* Check input arguments */
    if (descA.nb != descA.mb) {
        magma_error("MAGMA_zplghe_Tile", "only square tiles supported");
        return magma_request_fail(sequence, request, MAGMA_ERR_ILLEGAL_VALUE);
    }

    /* Quick return */
    if (min( descA.m, descA.n ) == 0)
        return MAGMA_SUCCESS;

    magma_pzplghe( bump, A, seed, sequence, request);

    return MAGMA_SUCCESS;
}
