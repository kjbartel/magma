/**	
 *
 * @file ztrtri.c
 *
 *  MAGMA computational routines
 *  MAGMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.3.1
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
 *  MAGMA_ztrtri - Computes the inverse of a complex upper or lower
 *  triangular matrix A.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          = MagmaUpper: Upper triangle of A is stored;
 *          = MagmaLower: Lower triangle of A is stored.
 *
 * @param[in] diag
 *          = MagmaNonUnit: A is non-unit triangular;
 *          = MagmaUnit:    A is unit triangular.
 *
 * @param[in] N
 *          The order of the matrix A. N >= 0.
 *
 * @param[in,out] A
 *          On entry, the triangular matrix A.  If UPLO = 'U', the
 *          leading N-by-N upper triangular part of the array A
 *          contains the upper triangular matrix, and the strictly
 *          lower triangular part of A is not referenced.  If UPLO =
 *          'L', the leading N-by-N lower triangular part of the array
 *          A contains the lower triangular matrix, and the strictly
 *          upper triangular part of A is not referenced.  If DIAG =
 *          'U', the diagonal elements of A are also not referenced and
 *          are assumed to be 1.  On exit, the (triangular) inverse of
 *          the original matrix, in the same storage format.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,N).
 *
 *******************************************************************************
 *
 * @return
 *          \retval MAGMA_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *          \retval >0 if i, A(i,i) is exactly zero.  The triangular
 *               matrix is singular and its inverse can not be computed.
 *
 *******************************************************************************
 *
 * @sa MAGMA_ztrtri_Tile
 * @sa MAGMA_ztrtri_Tile_Async
 * @sa MAGMA_ctrtri
 * @sa MAGMA_dtrtri
 * @sa MAGMA_strtri
 * @sa MAGMA_zpotri
 *
 ******************************************************************************/


int MAGMA_ztrtri(PLASMA_enum uplo, PLASMA_enum diag, int N,
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
        magma_fatal_error("MAGMA_ztrtri", "MAGMA not initialized");
        return MAGMA_ERR_NOT_INITIALIZED;
    }
    /* Check input arguments */
    if (uplo != PlasmaUpper && uplo != PlasmaLower) {
        magma_error("MAGMA_ztrtri", "illegal value of uplo");
        return -1;
    }
    
    if (diag != PlasmaUnit && diag != PlasmaNonUnit) {
        magma_error("MAGMA_ztrtri", "illegal value of diag");
        return -2;
    }

    if (N < 0) {
         magma_error("MAGMA_ztrtri", "illegal value of N");
        return -3;
    }
    if (LDA < max(1, N)) {
        magma_error("MAGMA_ztrtri", "illegal value of LDA");
        return -5;
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
    MAGMA_ztrtri_Tile_Async(uplo, diag, &descA, sequence, &request);

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
 *  MAGMA_ztrtri_Tile - Computes the inverse of a complex upper or
 *  lower triangular matrix A.
 *  Tile equivalent of MAGMA_ztrtri().
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          = MagmaUpper: Upper triangle of A is stored;
 *          = MagmaLower: Lower triangle of A is stored.
 *
 * @param[in] diag
 *          = MagmaNonUnit: A is non-unit triangular;
 *          = MagmaUnit:    A us unit triangular.
 *
 * @param[in] A
 *          On entry, the triangular matrix A.  If UPLO = 'U', the
 *          leading N-by-N upper triangular part of the array A
 *          contains the upper triangular matrix, and the strictly
 *          lower triangular part of A is not referenced.  If UPLO =
 *          'L', the leading N-by-N lower triangular part of the array
 *          A contains the lower triangular matrix, and the strictly
 *          upper triangular part of A is not referenced.  If DIAG =
 *          'U', the diagonal elements of A are also not referenced and
 *          are assumed to be 1.  On exit, the (triangular) inverse of
 *          the original matrix, in the same storage format.
 * 
 *******************************************************************************
 *
 * @return
 *          \retval MAGMA_SUCCESS successful exit
 *          \retval >0 if i, A(i,i) is exactly zero.  The triangular
 *               matrix is singular and its inverse can not be computed.
 *
 *******************************************************************************
 *
 * @sa MAGMA_ztrtri
 * @sa MAGMA_ztrtri_Tile_Async
 * @sa MAGMA_ctrtri_Tile
 * @sa MAGMA_dtrtri_Tile
 * @sa MAGMA_strtri_Tile
 * @sa MAGMA_zpotri_Tile
 *
 ******************************************************************************/

int MAGMA_ztrtri_Tile(PLASMA_enum uplo, PLASMA_enum diag, magma_desc_t *A)
{
    magma_context_t *magma;
    magma_sequence_t *sequence = NULL;
    magma_request_t request = MAGMA_REQUEST_INITIALIZER;
    int status;

    magma = magma_context_self();
    if (magma == NULL) {
        magma_fatal_error("MAGMA_ztrtri_Tile", "MAGMA not initialized");
        return MAGMA_ERR_NOT_INITIALIZED;
    }
    magma_sequence_create(magma, &sequence);
    MAGMA_ztrtri_Tile_Async(uplo, diag, A, sequence, &request);
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
 *  MAGMA_ztrtri_Tile_Async - Computes the inverse of a complex upper or lower
 *  triangular matrix A.
 *  Non-blocking equivalent of MAGMA_ztrtri_Tile().
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
 * @sa MAGMA_ztrtri
 * @sa MAGMA_ztrtri_Tile
 * @sa MAGMA_ctrtri_Tile_Async
 * @sa MAGMA_dtrtri_Tile_Async
 * @sa MAGMA_strtri_Tile_Async
 * @sa MAGMA_zpotri_Tile_Async
 *
 ******************************************************************************/

int MAGMA_ztrtri_Tile_Async(PLASMA_enum uplo, PLASMA_enum diag, magma_desc_t *A,
                              magma_sequence_t *sequence, magma_request_t *request)
{
    magma_context_t *magma;

    magma = magma_context_self();
    if (magma == NULL) {
        magma_fatal_error("MAGMA_ztrtri_Tile", "MAGMA not initialized");
        return MAGMA_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        magma_fatal_error("MAGMA_ztrtri_Tile", "NULL sequence");
        return MAGMA_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        magma_fatal_error("MAGMA_ztrtri_Tile", "NULL request");
        return MAGMA_ERR_UNALLOCATED;
    }
    /* Check sequence status */
    if (sequence->status == MAGMA_SUCCESS)
        request->status = MAGMA_SUCCESS;
    else
        return magma_request_fail(sequence, request, MAGMA_ERR_SEQUENCE_FLUSHED);

    /* Check input arguments */
    if (A->nb != A->mb) {
        magma_error("MAGMA_ztrtri_Tile", "only square tiles supported");
        return magma_request_fail(sequence, request, MAGMA_ERR_ILLEGAL_VALUE);
    }
    if (uplo != PlasmaUpper && uplo != PlasmaLower) {
        magma_error("MAGMA_ztrtri_Tile", "illegal value of uplo");
        return magma_request_fail(sequence, request, -1);
    }

    if (diag != PlasmaUnit && diag != PlasmaNonUnit) {
        magma_error("MAGMA_ztrtri_Tile", "illegal value of diag");
        return magma_request_fail(sequence, request, -2);
    }

    magma_pztrtri(uplo, diag, A, sequence, request);

    return MAGMA_SUCCESS;
}
