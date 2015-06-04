/**
 *
 * @file zpotri.c
 *
 *  MAGMA computational routines
 *  MAGMA is a software package provided by Univ. of Tennessee,
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
 *  MAGMA_zpotri - Computes the inverse of a complex Hermitian positive definite
 *  matrix A using the Cholesky factorization A = U**H*U or A = L*L**H
 *  computed by MAGMA_zpotrf.  
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          = SplagmaUpper: Upper triangle of A is stored;
 *          = SplagmaLower: Lower triangle of A is stored.
 *
 * @param[in] N
 *          The order of the matrix A. N >= 0.
 *
 * @param[in,out] A
 *          On entry, the triangular factor U or L from the Cholesky
 *          factorization A = U**H*U or A = L*L**H, as computed by
 *          MAGMA_zpotrf.
 *          On exit, the upper or lower triangle of the (Hermitian)
 *          inverse of A, overwriting the input factor U or L.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,N).
 *
 *******************************************************************************
 *
 * @return
 *          \retval MAGMA_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *          \retval >0 if i, the (i,i) element of the factor U or L is
 *                zero, and the inverse could not be computed.
 *
 *******************************************************************************
 *
 * @sa MAGMA_zpotri_Tile
 * @sa MAGMA_zpotri_Tile_Async
 * @sa MAGMA_cpotri
 * @sa MAGMA_dpotri
 * @sa MAGMA_spotri
 * @sa MAGMA_zpotrf
 *                           
 ******************************************************************************/

int MAGMA_zpotri(PLASMA_enum uplo, int N,
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
        magma_fatal_error("MAGMA_zpotri", "MAGMA not initialized");
        return MAGMA_ERR_NOT_INITIALIZED;
    }
    /* Check input arguments */
    if (uplo != PlasmaUpper && uplo != PlasmaLower) {
        magma_error("MAGMA_zpotri", "illegal value of uplo");
        return -1;
    }
    if (N < 0) {
         magma_error("MAGMA_zpotri", "illegal value of N");
        return -2;
    }
    if (LDA < max(1, N)) {
        magma_error("MAGMA_zpotri", "illegal value of LDA");
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
    MAGMA_zpotri_Tile_Async(uplo, &descA, sequence, &request);

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
 * 
 * @ingroup MAGMA_Complex64_t_Tile
 * 
 *  MAGMA_zpotri_Tile - Computes the inverse of a complex Hermitian
 *  positive definite matrix A using the Cholesky factorization
 *  A = U**H*U or A = L*L**H computed by MAGMA_zpotrf.
 *  Tile equivalent of MAGMA_zpotri().
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
 *          On entry, the triangular factor U or L from the Cholesky
 *          factorization A = U**H*U or A = L*L**H, as computed by
 *          MAGMA_zpotrf.
 *          On exit, the upper or lower triangle of the (Hermitian)
 *          inverse of A, overwriting the input factor U or L.
 *
 *******************************************************************************
 *
 * @return
 *          \retval MAGMA_SUCCESS successful exit
 *          \retval >0 if i, the leading minor of order i of A is not
 *               positive definite, so the factorization could not be
 *               completed, and the solution has not been computed.
 *
 *******************************************************************************
 *
 * @sa MAGMA_zpotri
 * @sa MAGMA_zpotri_Tile_Async
 * @sa MAGMA_cpotri_Tile
 * @sa MAGMA_dpotri_Tile
 * @sa MAGMA_spotri_Tile
 * @sa MAGMA_zpotrf_Tile
 *
 ******************************************************************************/

int MAGMA_zpotri_Tile(PLASMA_enum uplo, magma_desc_t *A)
{
    magma_context_t *magma;
    magma_sequence_t *sequence = NULL;
    magma_request_t request = MAGMA_REQUEST_INITIALIZER;
    int status;

    magma = magma_context_self();
    if (magma == NULL) {
        magma_fatal_error("MAGMA_zpotri_Tile", "MAGMA not initialized");
        return MAGMA_ERR_NOT_INITIALIZED;
    }
    magma_sequence_create(magma, &sequence);
    MAGMA_zpotri_Tile_Async(uplo, A, sequence, &request);
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
 *  MAGMA_zpotri_Tile_Async - Computes the inverse of a complex Hermitian
 *  positive definite matrix A using the Cholesky factorization A = U**H*U
 *  or A = L*L**H computed by MAGMA_zpotrf.
 *  Non-blocking equivalent of MAGMA_zpotri_Tile().
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
 * @sa MAGMA_zpotri
 * @sa MAGMA_zpotri_Tile
 * @sa MAGMA_cpotri_Tile_Async
 * @sa MAGMA_dpotri_Tile_Async
 * @sa MAGMA_spotri_Tile_Async
 * @sa MAGMA_zpotrf_Tile_Async
 *
 ******************************************************************************/

int MAGMA_zpotri_Tile_Async(PLASMA_enum uplo, magma_desc_t *A,
                            magma_sequence_t *sequence, magma_request_t *request)
{
    magma_context_t *magma;

    magma = magma_context_self();
    if (magma == NULL) {
        magma_fatal_error("MAGMA_zpotri_Tile", "MAGMA not initialized");
        return MAGMA_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        magma_fatal_error("MAGMA_zpotri_Tile", "NULL sequence");
        return MAGMA_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        magma_fatal_error("MAGMA_zpotri_Tile", "NULL request");
        return MAGMA_ERR_UNALLOCATED;
    }
    /* Check sequence status */
    if (sequence->status == MAGMA_SUCCESS)
        request->status = MAGMA_SUCCESS;
    else
        return magma_request_fail(sequence, request, MAGMA_ERR_SEQUENCE_FLUSHED);

    /* Check input arguments */
    if (A->nb != A->mb) {
        magma_error("MAGMA_zpotri_Tile", "only square tiles supported");
        return magma_request_fail(sequence, request, MAGMA_ERR_ILLEGAL_VALUE);
    }
    if (uplo != PlasmaUpper && uplo != PlasmaLower) {
        magma_error("MAGMA_zpotri_Tile", "illegal value of uplo");
        return magma_request_fail(sequence, request, -1);
    }
    
    magma_pztrtri(uplo, PlasmaNonUnit, A, sequence, request);
    magma_pzlauum(uplo, A, sequence, request);

    return MAGMA_SUCCESS;
}
