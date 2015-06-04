/**
 *
 * @file ztrsm.c
 *
 *  PLASMA computational routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.3.1
 * @author Jakub Kurzak
 * @date 2010-11-15
 * @precisions normal z -> s d c
 *
 **/
#include "common.h"

/***************************************************************************//**
 *
 * @ingroup MAGMA_Complex64_t
 *
 *  MAGMA_ztrsm - Computes triangular solve A*X = B or X*A = B.
 *
 *******************************************************************************
 *
 * @param[in] side
 *          Specifies whether A appears on the left or on the right of X:
 *          = PlasmaLeft:  A*X = B
 *          = PlasmaRight: X*A = B
 *
 * @param[in] uplo
 *          Specifies whether the matrix A is upper triangular or lower triangular:
 *          = PlasmaUpper: Upper triangle of A is stored;
 *          = PlasmaLower: Lower triangle of A is stored.
 *
 * @param[in] transA
 *          Specifies whether the matrix A is transposed, not transposed or conjugate transposed:
 *          = PlasmaNoTrans:   A is transposed;
 *          = PlasmaTrans:     A is not transposed;
 *          = PlasmaConjTrans: A is conjugate transposed.
 *
 * @param[in] diag
 *          Specifies whether or not A is unit triangular:
 *          = PlasmaNonUnit: A is non unit;
 *          = PlasmaUnit:    A us unit.
 *
 * @param[in] N
 *          The order of the matrix A. N >= 0.
 *
 * @param[in] NRHS
 *          The number of right hand sides, i.e., the number of columns of the matrix B. NRHS >= 0.
 *
 * @param[in] alpha
 *          alpha specifies the scalar alpha.
 *
 * @param[in] A
 *          The triangular matrix A. If uplo = PlasmaUpper, the leading N-by-N upper triangular
 *          part of the array A contains the upper triangular matrix, and the strictly lower
 *          triangular part of A is not referenced. If uplo = PlasmaLower, the leading N-by-N
 *          lower triangular part of the array A contains the lower triangular matrix, and the
 *          strictly upper triangular part of A is not referenced. If diag = PlasmaUnit, the
 *          diagonal elements of A are also not referenced and are assumed to be 1.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,N).
 *
 * @param[in,out] B
 *          On entry, the N-by-NRHS right hand side matrix B.
 *          On exit, if return value = 0, the N-by-NRHS solution matrix X.
 *
 * @param[in] LDB
 *          The leading dimension of the array B. LDB >= max(1,N).
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *
 *******************************************************************************
 *
 * @sa MAGMA_ztrsm_Tile
 * @sa MAGMA_ztrsm_Tile_Async
 * @sa MAGMA_ctrsm
 * @sa MAGMA_dtrsm
 * @sa MAGMA_strsm
 *
 ******************************************************************************/
int MAGMA_ztrsm(PLASMA_enum side, PLASMA_enum uplo,
                 PLASMA_enum transA, PLASMA_enum diag,
                 int N, int NRHS, PLASMA_Complex64_t alpha,
                 PLASMA_Complex64_t *A, int LDA,
                 PLASMA_Complex64_t *B, int LDB)
{
    int NB, NA;
    int status;
    magma_context_t *magma;
    magma_sequence_t *sequence = NULL;
    magma_request_t request = MAGMA_REQUEST_INITIALIZER;
    magma_desc_t descA, descB;

    magma = magma_context_self();
    if (magma == NULL) {
        magma_fatal_error("MAGMA_ztrsm", "MAGMA not initialized");
        return MAGMA_ERR_NOT_INITIALIZED;
    }
    /* Check input arguments */
    if (side != PlasmaLeft && side != PlasmaRight) {
        magma_error("MAGMA_ztrsm", "illegal value of side");
        return -1;
    }
    if (uplo != PlasmaUpper && uplo != PlasmaLower) {
        magma_error("MAGMA_ztrsm", "illegal value of uplo");
        return -2;
    }
    if (transA != PlasmaConjTrans && transA != PlasmaNoTrans && transA != PlasmaTrans ) {
        magma_error("MAGMA_ztrsm", "illegal value of transA");
        return -3;
    }
    if (diag != PlasmaUnit && diag != PlasmaNonUnit) {
        magma_error("MAGMA_ztrsm", "illegal value of diag");
        return -4;
    }
    if (N < 0) {
        magma_error("MAGMA_ztrsm", "illegal value of N");
        return -5;
    }
    if (NRHS < 0) {
        magma_error("MAGMA_ztrsm", "illegal value of NRHS");
        return -6;
    }
    if (LDA < max(1, N)) {
        magma_error("MAGMA_ztrsm", "illegal value of LDA");
        return -8;
    }
    if (LDB < max(1, N)) {
        magma_error("MAGMA_ztrsm", "illegal value of LDB");
        return -10;
    }
    /* Quick return */
    if (min(N, NRHS) == 0)
        return MAGMA_SUCCESS;

    /* Tune NB depending on M, N && NRHS; Set NBNB */
    /* status = magma_tune(MAGMA_FUNC_ZPOSV, N, N, NRHS); */
    /* if (status != MAGMA_SUCCESS) { */
    /*     magma_error("MAGMA_ztrsm", "magma_tune() failed"); */
    /*     return status; */
    /* } */

    /* Set NT && NTRHS */
    NB = MAGMA_NB;
    if (side == PlasmaLeft) {
      NA = N;
    } else {
      NA = NRHS;
    }

    magma_sequence_create(magma, &sequence);

    /* if ( MAGMA_TRANSLATION == MAGMA_OUTOFPLACE ) { */
        magma_zooplap2tile( descA, A, NB, NB, LDA, NA,   0, 0, NA, NA,   magma_desc_mat_free(&(descA)) );
        magma_zooplap2tile( descB, B, NB, NB, LDB, NRHS, 0, 0, N,  NRHS, magma_desc_mat_free(&(descA)); magma_desc_mat_free(&(descB)));
    /* } else { */
    /*     magma_ziplap2tile( descA, A, NB, NB, LDA, NA,   0, 0, NA, NA  ); */
    /*     magma_ziplap2tile( descB, B, NB, NB, LDB, NRHS, 0, 0, N,  NRHS); */
    /* } */

    /* Call the tile interface */
    MAGMA_ztrsm_Tile_Async(
        side, uplo, transA, diag, alpha, &descA, &descB, sequence, &request);

    /* if ( MAGMA_TRANSLATION == MAGMA_OUTOFPLACE ) { */
        magma_zooptile2lap( descB, B, NB, NB, LDB, NRHS );
        morse_barrier( magma );
        magma_desc_mat_free(&descA);
        magma_desc_mat_free(&descB);
    /* } else { */
    /*     magma_ziptile2lap( descA, A, NB, NB, LDA, NA   ); */
    /*     magma_ziptile2lap( descB, B, NB, NB, LDB, NRHS ); */
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
 *  MAGMA_ztrsm_Tile - Computes triangular solve.
 *  Tile equivalent of MAGMA_ztrsm().
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] side
 *          Specifies whether A appears on the left or on the right of X:
 *          = PlasmaLeft:  A*X = B
 *          = PlasmaRight: X*A = B
 *
 * @param[in] uplo
 *          Specifies whether the matrix A is upper triangular or lower triangular:
 *          = PlasmaUpper: Upper triangle of A is stored;
 *          = PlasmaLower: Lower triangle of A is stored.
 *
 * @param[in] transA
 *          Specifies whether the matrix A is transposed, not transposed or conjugate transposed:
 *          = PlasmaNoTrans:   A is transposed;
 *          = PlasmaTrans:     A is not transposed;
 *          = PlasmaConjTrans: A is conjugate transposed.
 *
 * @param[in] diag
 *          Specifies whether or not A is unit triangular:
 *          = PlasmaNonUnit: A is non unit;
 *          = PlasmaUnit:    A us unit.
 *
 * @param[in] alpha
 *          alpha specifies the scalar alpha.
 *
 * @param[in] A
 *          The triangular matrix A. If uplo = PlasmaUpper, the leading N-by-N upper triangular
 *          part of the array A contains the upper triangular matrix, and the strictly lower
 *          triangular part of A is not referenced. If uplo = PlasmaLower, the leading N-by-N
 *          lower triangular part of the array A contains the lower triangular matrix, and the
 *          strictly upper triangular part of A is not referenced. If diag = PlasmaUnit, the
 *          diagonal elements of A are also not referenced and are assumed to be 1.
 *
 * @param[in,out] B
 *          On entry, the N-by-NRHS right hand side matrix B.
 *          On exit, if return value = 0, the N-by-NRHS solution matrix X.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa MAGMA_ztrsm
 * @sa MAGMA_ztrsm_Tile_Async
 * @sa MAGMA_ctrsm_Tile
 * @sa MAGMA_dtrsm_Tile
 * @sa MAGMA_strsm_Tile
 *
 ******************************************************************************/
int MAGMA_ztrsm_Tile(PLASMA_enum side, PLASMA_enum uplo,
                      PLASMA_enum transA, PLASMA_enum diag,
                      PLASMA_Complex64_t alpha, magma_desc_t *A, magma_desc_t *B)
{
    magma_context_t *magma;
    magma_sequence_t *sequence = NULL;
    magma_request_t request = MAGMA_REQUEST_INITIALIZER;
    int status;

    magma = magma_context_self();
    if (magma == NULL) {
        magma_fatal_error("MAGMA_ztrsm_Tile", "MAGMA not initialized");
        return MAGMA_ERR_NOT_INITIALIZED;
    }
    magma_sequence_create(magma, &sequence);
    MAGMA_ztrsm_Tile_Async(side, uplo, transA, diag, alpha, A, B, sequence, &request);
    morse_desc_getoncpu( A );
    morse_desc_getoncpu( B );
    morse_barrier( magma );
    status = sequence->status;
    magma_sequence_destroy(magma, sequence);
    return status;
}

/***************************************************************************//**
 *
 * @ingroup MAGMA_Complex64_t_Tile_Async
 *
 *  MAGMA_ztrsm_Tile_Async - Computes triangular solve.
 *  Non-blocking equivalent of MAGMA_ztrsm_Tile().
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
 * @sa MAGMA_ztrsm
 * @sa MAGMA_ztrsm_Tile
 * @sa MAGMA_ctrsm_Tile_Async
 * @sa MAGMA_dtrsm_Tile_Async
 * @sa MAGMA_strsm_Tile_Async
 *
 ******************************************************************************/
int MAGMA_ztrsm_Tile_Async(PLASMA_enum side, PLASMA_enum uplo,
                            PLASMA_enum transA, PLASMA_enum diag,
                            PLASMA_Complex64_t alpha, magma_desc_t *A, magma_desc_t *B,
                            magma_sequence_t *sequence, magma_request_t *request)
{
    magma_context_t *magma;

    magma = magma_context_self();
    if (magma == NULL) {
        magma_fatal_error("MAGMA_ztrsm_Tile", "MAGMA not initialized");
        return MAGMA_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        magma_fatal_error("MAGMA_ztrsm_Tile", "NULL sequence");
        return MAGMA_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        magma_fatal_error("MAGMA_ztrsm_Tile", "NULL request");
        return MAGMA_ERR_UNALLOCATED;
    }
    /* Check sequence status */
    if (sequence->status == MAGMA_SUCCESS)
        request->status = MAGMA_SUCCESS;
    else
        return magma_request_fail(sequence, request, MAGMA_ERR_SEQUENCE_FLUSHED);

    /* Check descriptors for correctness */
    if (magma_desc_check( A ) != MAGMA_SUCCESS) {
        magma_error("MAGMA_ztrsm_Tile", "invalid first descriptor");
        return magma_request_fail(sequence, request, MAGMA_ERR_ILLEGAL_VALUE);
    }
    if (magma_desc_check( B ) != MAGMA_SUCCESS) {
        magma_error("MAGMA_ztrsm_Tile", "invalid second descriptor");
        return magma_request_fail(sequence, request, MAGMA_ERR_ILLEGAL_VALUE);
    }
    /* Check input arguments */
    if (A->nb != A->mb || B->nb != B->mb) {
        magma_error("MAGMA_ztrsm_Tile", "only square tiles supported");
        return magma_request_fail(sequence, request, MAGMA_ERR_ILLEGAL_VALUE);
    }
    if (side != PlasmaLeft && side != PlasmaRight) {
        magma_error("MAGMA_ztrsm_Tile", "illegal value of side");
        return magma_request_fail(sequence, request, -1);
    }
    if (uplo != PlasmaUpper && uplo != PlasmaLower) {
        magma_error("MAGMA_ztrsm_Tile", "illegal value of uplo");
        return magma_request_fail(sequence, request, -2);
    }
    if (transA != PlasmaConjTrans && transA != PlasmaNoTrans && transA != PlasmaTrans) {
        magma_error("MAGMA_ztrsm_Tile", "illegal value of transA");
        return magma_request_fail(sequence, request, -3);
    }
    if (diag != PlasmaUnit && diag != PlasmaNonUnit) {
        magma_error("MAGMA_ztrsm_Tile", "illegal value of diag");
        return magma_request_fail(sequence, request, -4);
    }

    /* Quick return */
    magma_pztrsm( side, uplo, transA, diag,
                  alpha, A, B,
                  sequence, request);

    return PLASMA_SUCCESS;
}
