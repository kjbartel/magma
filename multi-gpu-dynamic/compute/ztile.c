/**
 *
 * @file ztile.c
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.3.1
 * @author Jakub Kurzak
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @precisions normal z -> s d c
 *
 **/
#include "common.h"

/***************************************************************************//**
 *
 * @ingroup MAGMA_Complex64_t
 *
 *  MAGMA_zLapack_to_Tile - Conversion from LAPACK layout to tile layout.
 *
 *******************************************************************************
 *
 * @param[in] Af77
 *          LAPACK matrix.
 *
 * @param[in] LDA
 *          The leading dimension of the matrix Af77.
 *
 * @param[in,out] A
 *          Descriptor of the MAGMA matrix in tile layout.
 *          If MAGMA_TRANSLATION_MODE is set to MAGMA_INPLACE,
 *          A->mat is not used and set to Af77 when returns, else if
 *          MAGMA_TRANSLATION_MODE is set to MAGMA_OUTOFPLACE,
 *          A->mat has to be allocated before.
 *
 *******************************************************************************
 *
 * @return
 *          \retval MAGMA_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa MAGMA_zLapack_to_Tile_Async
 * @sa MAGMA_zTile_to_Lapack
 * @sa MAGMA_cLapack_to_Tile
 * @sa MAGMA_dLapack_to_Tile
 * @sa MAGMA_sLapack_to_Tile
 *
 ******************************************************************************/
int MAGMA_zLapack_to_Tile(PLASMA_Complex64_t *Af77, int LDA, magma_desc_t *A)
{
    magma_context_t *magma;
    magma_sequence_t *sequence = NULL;
    magma_request_t request;
    int status;

    magma = magma_context_self();
    if (magma == NULL) {
        magma_fatal_error("MAGMA_zLapack_to_Tile", "MAGMA not initialized");
        return MAGMA_ERR_NOT_INITIALIZED;
    }
    /* Check descriptor for correctness */
    if (magma_desc_check( A ) != MAGMA_SUCCESS) {
        magma_error("MAGMA_zLapack_to_Tile", "invalid descriptor");
        return MAGMA_ERR_ILLEGAL_VALUE;
    }
    magma_sequence_create(magma, &sequence);

    magma_pzlapack_to_tile( Af77, LDA, A, sequence, &request);

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
 *  MAGMA_zLapack_to_Tile_Async - Conversion from LAPACK layout to tile layout.
 *  Non-blocking equivalent of MAGMA_zLapack_to_Tile().
 *  May return before the computation is finished.
 *  Allows for pipelining of operations ar runtime.
 *
 *
 *******************************************************************************
 *
 * @param[in] Af77
 *          LAPACK matrix.
 *
 * @param[in] LDA
 *          The leading dimension of the matrix Af77.
 *
 * @param[in,out] A
 *          Descriptor of the MAGMA matrix in tile layout.
 *          If MAGMA_TRANSLATION_MODE is set to MAGMA_INPLACE,
 *          A->mat is not used and set to Af77 when returns, else if
 *          MAGMA_TRANSLATION_MODE is set to MAGMA_OUTOFPLACE,
 *          A->mat has to be allocated before.
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
 * @sa MAGMA_zTile_to_Lapack_Async
 * @sa MAGMA_zLapack_to_Tile
 * @sa MAGMA_cLapack_to_Tile_Async
 * @sa MAGMA_dLapack_to_Tile_Async
 * @sa MAGMA_sLapack_to_Tile_Async
 *
 ******************************************************************************/
int MAGMA_zLapack_to_Tile_Async(PLASMA_Complex64_t *Af77, int LDA, magma_desc_t *A,
                                  magma_sequence_t *sequence, magma_request_t *request)
{
    magma_context_t *magma;

    magma = magma_context_self();
    if (magma == NULL) {
        magma_fatal_error("MAGMA_zLapack_to_Tile", "MAGMA not initialized");
        return MAGMA_ERR_NOT_INITIALIZED;
    }
    /* Check descriptor for correctness */
    if (magma_desc_check( A ) != MAGMA_SUCCESS) {
        magma_error("MAGMA_zLapack_to_Tile", "invalid descriptor");
        return MAGMA_ERR_ILLEGAL_VALUE;
    }

    magma_pzlapack_to_tile( Af77, LDA, A, sequence, request);

    return MAGMA_SUCCESS;
}

/***************************************************************************//**
 *
 * @ingroup MAGMA_Complex64_t
 *
 *  MAGMA_Tile_to_Lapack - Conversion from tile layout to LAPACK layout.
 *
 *******************************************************************************
 *
 * @param[in] A
 *          Descriptor of the MAGMA matrix in tile layout.
 *
 * @param[in,out] Af77
 *          LAPACK matrix.
 *          If MAGMA_TRANSLATION_MODE is set to MAGMA_INPLACE,
 *          Af77 has to be A->mat, else if
 *          MAGMA_TRANSLATION_MODE is set to MAGMA_OUTOFPLACE,
 *          Af77 has to be allocated before.
 *
 * @param[in] LDA
 *          The leading dimension of the matrix Af77.
 *
 *******************************************************************************
 *
 * @return
 *          \retval MAGMA_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa MAGMA_zTile_to_Lapack_Async
 * @sa MAGMA_zLapack_to_Tile
 * @sa MAGMA_cTile_to_Lapack
 * @sa MAGMA_dTile_to_Lapack
 * @sa MAGMA_sTile_to_Lapack
 *
******************************************************************************/
int MAGMA_zTile_to_Lapack(magma_desc_t *A, PLASMA_Complex64_t *Af77, int LDA)
{
    magma_context_t *magma;
    magma_sequence_t *sequence = NULL;
    magma_request_t request;
    int status;

    magma = magma_context_self();
    if (magma == NULL) {
        magma_fatal_error("MAGMA_zTile_to_Lapack", "MAGMA not initialized");
        return MAGMA_ERR_NOT_INITIALIZED;
    }
    /* Check descriptor for correctness */
    if (magma_desc_check( A ) != MAGMA_SUCCESS) {
        magma_error("MAGMA_zTile_to_Lapack", "invalid descriptor");
        return MAGMA_ERR_ILLEGAL_VALUE;
    }
    magma_sequence_create(magma, &sequence);

    magma_pztile_to_lapack( A, Af77, LDA, sequence, &request);
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
 *  MAGMA_zTile_to_Lapack_Async - Conversion from LAPACK layout to tile layout.
 *  Non-blocking equivalent of MAGMA_zTile_to_Lapack().
 *  May return before the computation is finished.
 *  Allows for pipelining of operations ar runtime.
 *
 *
 *******************************************************************************
 *
 * @param[in] A
 *          Descriptor of the MAGMA matrix in tile layout.
 *
 * @param[in,out] Af77
 *          LAPACK matrix.
 *          If MAGMA_TRANSLATION_MODE is set to MAGMA_INPLACE,
 *          Af77 has to be A->mat, else if
 *          MAGMA_TRANSLATION_MODE is set to MAGMA_OUTOFPLACE,
 *          Af77 has to be allocated before.
 *
 * @param[in] LDA
 *          The leading dimension of the matrix Af77.
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
 * @sa MAGMA_zLapack_to_Tile_Async
 * @sa MAGMA_zTile_to_Lapack
 * @sa MAGMA_cTile_to_Lapack_Async
 * @sa MAGMA_dTile_to_Lapack_Async
 * @sa MAGMA_sTile_to_Lapack_Async
 *
 ******************************************************************************/
int MAGMA_zTile_to_Lapack_Async(magma_desc_t *A, PLASMA_Complex64_t *Af77, int LDA,
                                magma_sequence_t *sequence, magma_request_t *request)
{
    magma_context_t *magma;

    magma = magma_context_self();
    if (magma == NULL) {
        magma_fatal_error("MAGMA_zTile_to_Lapack", "MAGMA not initialized");
        return MAGMA_ERR_NOT_INITIALIZED;
    }
    /* Check descriptor for correctness */
    if (magma_desc_check( A ) != MAGMA_SUCCESS) {
        magma_error("MAGMA_zTile_to_Lapack", "invalid descriptor");
        return MAGMA_ERR_ILLEGAL_VALUE;
    }

    magma_pztile_to_lapack( A, Af77, LDA, sequence, request );

    return MAGMA_SUCCESS;
}
