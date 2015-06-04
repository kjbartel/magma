/**
 *
 * @file zlauum.c
 *
 *  PLASMA computational routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.4.1
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
 *  MAGMA_zlauum - Computes the product U * U' or L' * L, where the triangular
 *  factor U or L is stored in the upper or lower triangular part of
 *  the array A.
 *
 *  If UPLO = 'U' or 'u' then the upper triangle of the result is stored,
 *  overwriting the factor U in A.
 *  If UPLO = 'L' or 'l' then the lower triangle of the result is stored,
 *  overwriting the factor L in A.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          = SplagmaUpper: Upper triangle of A is stored;
 *          = SplagmaLower: Lower triangle of A is stored.
 *
 * @param[in] N
 *          The order of the triangular factor U or L.  N >= 0.
 *
 * @param[in,out] A
 *          On entry, the triangular factor U or L.
 *          On exit, if UPLO = 'U', the upper triangle of A is
 *          overwritten with the upper triangle of the product U * U';
 *          if UPLO = 'L', the lower triangle of A is overwritten with
 *          the lower triangle of the product L' * L.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,N).
 *
 *******************************************************************************
 *
 * @return
 *          \retval MAGMA_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *
 *******************************************************************************
 *
 * @sa MAGMA_zlauum_Tile
 * @sa MAGMA_zlauum_Tile_Async
 * @sa MAGMA_clauum
 * @sa MAGMA_dlauum
 * @sa MAGMA_slauum
 * @sa MAGMA_zpotri
 *                                              
 ******************************************************************************/

int MAGMA_zlauum(PLASMA_enum uplo, int N,
                   PLASMA_Complex64_t *A, int LDA)
{
    int NB;
    int status;
    magma_context_t *magma;
    magma_sequence_t *sequence = NULL;
    magma_request_t request = MAGMA_REQUEST_INITIALIZER;
    magma_desc_t descA;

    magma = magma_context_self();
    if (magma == NULL) {
        magma_fatal_error("MAGMA_zlauum", "MAGMA not initialized");
        return MAGMA_ERR_NOT_INITIALIZED;
    }
    /* Check input arguments */
    if (uplo != PlasmaUpper && uplo != PlasmaLower) {
        magma_error("MAGMA_zlauum", "illegal value of uplo");
        return -1;
    }
    if (N < 0) {
         magma_error("MAGMA_zlauum", "illegal value of N");
        return -2;
    }
    if (LDA < max(1, N)) {
        magma_error("MAGMA_zlauum", "illegal value of LDA");
        return -4;
    }
    /* Quick return */
    if (max(N, 0) == 0)
        return MAGMA_SUCCESS;


    /* Set NT */
    NB   = MAGMA_NB;

    magma_sequence_create(magma, &sequence);

    /* if ( MAGMA_TRANSLATION == MAGMA_OUTOFPLACE ) { */
        magma_zooplap2tile( descA, A, NB, NB, LDA, N, 0, 0, N, N, magma_desc_mat_free(&(descA)) );
    /* } else { */
    /*     magma_ziplap2tile(  descA, A, NB, NB, LDA, N, 0, 0, N, N); */
    /* } */

    /* Call the tile interface */
    MAGMA_zlauum_Tile_Async(uplo, &descA, sequence, &request);

    /* if ( MAGMA_TRANSLATION == MAGMA_OUTOFPLACE ) { */
        magma_zooptile2lap( descA, A, NB, NB, LDA, N );
        morse_barrier( magma );
        morse_desc_getoncpu( &descA );
        magma_desc_mat_free(&descA);
    /* } else { */
    /*     magma_ziptile2lap( descA, A, NB, NB, LDA, N ); */
    /*     morse_barrier( magma ); */
    /* } */
        
    status = sequence->status;
    magma_sequence_destroy(magma, sequence);

    return status;
}

/***************************************************************************//**
 * @ingroup MAGMA_Complex64_t_Tile
 *  
 *  MAGMA_zlauum_Tile - Computes the product U * U' or L' * L, where
 *  the triangular factor U or L is stored in the upper or lower
 *  triangular part of the array A.
 *  Tile equivalent of MAGMA_zlauum().
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *          
 *******************************************************************************
 *            
 * @param[in] uplo
 *          = SplagmaUpper: Upper triangle of A is stored;
 *          = SplagmaLower: Lower triangle of A is stored.
 *                
 * @param[in] A
 *          On entry, the triangular factor U or L.
 *          On exit, if UPLO = 'U', the upper triangle of A is
 *          overwritten with the upper triangle of the product U * U';
 *          if UPLO = 'L', the lower triangle of A is overwritten with
 *          the lower triangle of the product L' * L.
 *                       
 *******************************************************************************
 *                         
 * @return
 *          \retval MAGMA_SUCCESS successful exit
 *                            
 *******************************************************************************
 *                              
 * @sa MAGMA_zlauum
 * @sa MAGMA_zlauum_Tile_Async
 * @sa MAGMA_clauum_Tile
 * @sa MAGMA_dlauum_Tile
 * @sa MAGMA_slauum_Tile
 * @sa MAGMA_zpotri_Tile
 *                                     
 ******************************************************************************/
int MAGMA_zlauum_Tile(PLASMA_enum uplo, magma_desc_t *A)
{
    magma_context_t *magma;
    magma_sequence_t *sequence = NULL;
    magma_request_t request = MAGMA_REQUEST_INITIALIZER;
    int status;

    magma = magma_context_self();
    if (magma == NULL) {
        magma_fatal_error("MAGMA_zlauum_Tile", "MAGMA not initialized");
        return MAGMA_ERR_NOT_INITIALIZED;
    }
    magma_sequence_create(magma, &sequence);
    MAGMA_zlauum_Tile_Async(uplo, A, sequence, &request);
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
 *  MAGMA_zlauum_Tile_Async - Computes the product U * U' or L' * L, where the
 *  triangular factor U or L is stored in the upper or lower triangular part of
 *  the array A.
 *  Non-blocking equivalent of MAGMA_zlauum_Tile().
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
 * @sa MAGMA_zlauum
 * @sa MAGMA_zlauum_Tile
 * @sa MAGMA_clauum_Tile_Async
 * @sa MAGMA_dlauum_Tile_Async
 * @sa MAGMA_slauum_Tile_Async
 * @sa MAGMA_zpotri_Tile_Async
 *
 ******************************************************************************/
int MAGMA_zlauum_Tile_Async(PLASMA_enum uplo, magma_desc_t *A,
                              magma_sequence_t *sequence, magma_request_t *request)
{
    magma_context_t *magma;

    magma = magma_context_self();
    if (magma == NULL) {
        magma_fatal_error("MAGMA_zlauum_Tile", "MAGMA not initialized");
        return MAGMA_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        magma_fatal_error("MAGMA_zlauum_Tile", "NULL sequence");
        return MAGMA_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        magma_fatal_error("MAGMA_zlauum_Tile", "NULL request");
        return MAGMA_ERR_UNALLOCATED;
    }
    /* Check sequence status */
    if (sequence->status == MAGMA_SUCCESS)
        request->status = MAGMA_SUCCESS;
    else
        return magma_request_fail(sequence, request, MAGMA_ERR_SEQUENCE_FLUSHED);

    /* Check input arguments */
    if (A->nb != A->mb) {
        magma_error("MAGMA_zlauum_Tile", "only square tiles supported");
        return magma_request_fail(sequence, request, MAGMA_ERR_ILLEGAL_VALUE);
    }
    if (uplo != PlasmaUpper && uplo != PlasmaLower) {
        magma_error("MAGMA_zlauum_Tile", "illegal value of uplo");
        return magma_request_fail(sequence, request, -1);
    }
    magma_pzlauum(uplo, A, sequence, request);

    return MAGMA_SUCCESS;
}
