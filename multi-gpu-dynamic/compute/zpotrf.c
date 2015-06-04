/**
 *
 * @file zpotrf.c
 *
 *  PLASMA computational routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
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
 *  MAGMA_zpotrf - Computes the Cholesky factorization of a symmetric positive definite
 *  (or Hermitian positive definite in the complex case) matrix A.
 *  The factorization has the form
 *
 *    \f[ A = \{_{L\times L^H, if uplo = PlasmaLower}^{U^H\times U, if uplo = PlasmaUpper} \f]
 *
 *  where U is an upper triangular matrix and L is a lower triangular matrix.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          = PlasmaUpper: Upper triangle of A is stored;
 *          = PlasmaLower: Lower triangle of A is stored.
 *
 * @param[in] N
 *          The order of the matrix A. N >= 0.
 *
 * @param[in,out] A
 *          On entry, the symmetric positive definite (or Hermitian) matrix A.
 *          If uplo = PlasmaUpper, the leading N-by-N upper triangular part of A
 *          contains the upper triangular part of the matrix A, and the strictly lower triangular
 *          part of A is not referenced.
 *          If UPLO = 'L', the leading N-by-N lower triangular part of A contains the lower
 *          triangular part of the matrix A, and the strictly upper triangular part of A is not
 *          referenced.
 *          On exit, if return value = 0, the factor U or L from the Cholesky factorization
 *          A = U**H*U or A = L*L**H.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,N).
 *
 *******************************************************************************
 *
 * @return
 *          \retval MAGMA_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *          \retval >0 if i, the leading minor of order i of A is not positive definite, so the
 *               factorization could not be completed, and the solution has not been computed.
 *
 *******************************************************************************
 *
 * @sa MAGMA_zpotrf_Tile
 * @sa MAGMA_zpotrf_Tile_Async
 * @sa MAGMA_cpotrf
 * @sa MAGMA_dpotrf
 * @sa MAGMA_spotrf
 * @sa MAGMA_zpotrs
 *
 ******************************************************************************/
int MAGMA_zpotrf(PLASMA_enum uplo, int N,
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
        magma_fatal_error("MAGMA_zpotrf", "MAGMA not initialized");
        return MAGMA_ERR_NOT_INITIALIZED;
    }
    /* Check input arguments */
    if (uplo != PlasmaUpper && uplo != PlasmaLower) {
        magma_error("MAGMA_zpotrf", "illegal value of uplo");
        return -1;
    }
    if (N < 0) {
         magma_error("MAGMA_zpotrf", "illegal value of N");
        return -2;
    }
    if (LDA < max(1, N)) {
        magma_error("MAGMA_zpotrf", "illegal value of LDA");
        return -4;
    }
    /* Quick return */
    if (max(N, 0) == 0)
        return MAGMA_SUCCESS;

    /* Tune NB depending on M, N & NRHS; Set NBNB */
    /* status = magma_tune(MAGMA_FUNC_ZPOSV, N, N, 0); */
    /* if (status != MAGMA_SUCCESS) { */
    /*     magma_error("MAGMA_zpotrf", "magma_tune() failed"); */
    /*     return status; */
    /* } */

    /* Set NT */
    NB   = MAGMA_NB;

    magma_sequence_create(magma, &sequence);

    /* if ( MAGMA_TRANSLATION == MAGMA_OUTOFPLACE ) { */
        magma_zooplap2tile( descA, A, NB, NB, LDA, N, 0, 0, N, N, magma_desc_mat_free(&(descA)) );
    /* } else { */
    /*     magma_ziplap2tile(  descA, A, NB, NB, LDA, N, 0, 0, N, N); */
    /* } */

    /* Call the tile interface */
    MAGMA_zpotrf_Tile_Async(uplo, &descA, sequence, &request);

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
 *  MAGMA_zpotrf_Tile - Computes the Cholesky factorization of a symmetric positive definite
 *  or Hermitian positive definite matrix.
 *  Tile equivalent of MAGMA_zpotrf().
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
 * @param[in] A
 *          On entry, the symmetric positive definite (or Hermitian) matrix A.
 *          If uplo = MagmaUpper, the leading N-by-N upper triangular part of A
 *          contains the upper triangular part of the matrix A, and the strictly lower triangular
 *          part of A is not referenced.
 *          If UPLO = 'L', the leading N-by-N lower triangular part of A contains the lower
 *          triangular part of the matrix A, and the strictly upper triangular part of A is not
 *          referenced.
 *          On exit, if return value = 0, the factor U or L from the Cholesky factorization
 *          A = U**H*U or A = L*L**H.
 *
 *******************************************************************************
 *
 * @return
 *          \retval MAGMA_SUCCESS successful exit
 *          \retval >0 if i, the leading minor of order i of A is not positive definite, so the
 *               factorization could not be completed, and the solution has not been computed.
 *
 *******************************************************************************
 *
 * @sa MAGMA_zpotrf
 * @sa MAGMA_zpotrf_Tile_Async
 * @sa MAGMA_cpotrf_Tile
 * @sa MAGMA_dpotrf_Tile
 * @sa MAGMA_spotrf_Tile
 * @sa MAGMA_zpotrs_Tile
 *
 ******************************************************************************/
int MAGMA_zpotrf_Tile(PLASMA_enum uplo, magma_desc_t *A)
{
    magma_context_t *magma;
    magma_sequence_t *sequence = NULL;
    magma_request_t request = MAGMA_REQUEST_INITIALIZER;
    int status;

    magma = magma_context_self();
    if (magma == NULL) {
        magma_fatal_error("MAGMA_zpotrf_Tile", "MAGMA not initialized");
        return MAGMA_ERR_NOT_INITIALIZED;
    }
    magma_sequence_create(magma, &sequence);
    MAGMA_zpotrf_Tile_Async(uplo, A, sequence, &request);
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
 *  MAGMA_zpotrf_Tile_Async - Computes the Cholesky factorization of a symmetric
 *  positive definite or Hermitian positive definite matrix.
 *  Non-blocking equivalent of MAGMA_zpotrf_Tile().
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
 * @sa MAGMA_zpotrf
 * @sa MAGMA_zpotrf_Tile
 * @sa MAGMA_cpotrf_Tile_Async
 * @sa MAGMA_dpotrf_Tile_Async
 * @sa MAGMA_spotrf_Tile_Async
 * @sa MAGMA_zpotrs_Tile_Async
 *
 ******************************************************************************/
int MAGMA_zpotrf_Tile_Async(PLASMA_enum uplo, magma_desc_t *A,
                              magma_sequence_t *sequence, magma_request_t *request)
{
    magma_context_t *magma;

    magma = magma_context_self();
    if (magma == NULL) {
        magma_fatal_error("MAGMA_zpotrf_Tile", "MAGMA not initialized");
        return MAGMA_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        magma_fatal_error("MAGMA_zpotrf_Tile", "NULL sequence");
        return MAGMA_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        magma_fatal_error("MAGMA_zpotrf_Tile", "NULL request");
        return MAGMA_ERR_UNALLOCATED;
    }
    /* Check sequence status */
    if (sequence->status == MAGMA_SUCCESS)
        request->status = MAGMA_SUCCESS;
    else
        return magma_request_fail(sequence, request, MAGMA_ERR_SEQUENCE_FLUSHED);

    /* Check descriptors for correctness */
    /* if (magma_desc_check( A ) != MAGMA_SUCCESS) { */
    /*     magma_error("MAGMA_zpotrf_Tile", "invalid descriptor"); */
    /*     return magma_request_fail(sequence, request, MAGMA_ERR_ILLEGAL_VALUE); */
    /* } */
    /* Check input arguments */
    if (A->nb != A->mb) {
        magma_error("MAGMA_zpotrf_Tile", "only square tiles supported");
        return magma_request_fail(sequence, request, MAGMA_ERR_ILLEGAL_VALUE);
    }
    if (uplo != PlasmaUpper && uplo != PlasmaLower) {
        magma_error("MAGMA_zpotrf_Tile", "illegal value of uplo");
        return magma_request_fail(sequence, request, -1);
    }
    /* Quick return */
/*
    if (max(N, 0) == 0)
        return MAGMA_SUCCESS;
*/
    magma_pzpotrf(uplo, A, sequence, request);

    return MAGMA_SUCCESS;
}
