/**
 *
 * @file zpotrs.c
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
 *  MAGMA_zpotrs - Solves a system of linear equations A * X = B with a symmetric positive
 *  definite (or Hermitian positive definite in the complex case) matrix A using the Cholesky
 *  factorization A = U**H*U or A = L*L**H computed by MAGMA_zpotrf.
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
 * @param[in] NRHS
 *          The number of right hand sides, i.e., the number of columns of the matrix B. NRHS >= 0.
 *
 * @param[in] A
 *          The triangular factor U or L from the Cholesky factorization A = U**H*U or A = L*L**H,
 *          computed by PLASMA_zpotrf.
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
 *          \retval MAGMA_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *
 *******************************************************************************
 *
 * @sa MAGMA_zpotrs_Tile
 * @sa MAGMA_zpotrs_Tile_Async
 * @sa MAGMA_cpotrs
 * @sa MAGMA_dpotrs
 * @sa MAGMA_spotrs
 * @sa MAGMA_zpotrf
 *
 ******************************************************************************/
int MAGMA_zpotrs(PLASMA_enum uplo, int N, int NRHS,
                  PLASMA_Complex64_t *A, int LDA,
                  PLASMA_Complex64_t *B, int LDB)
{
    int NB;
    int status;
    magma_context_t *magma;
    magma_sequence_t *sequence = NULL;
    magma_request_t request = MAGMA_REQUEST_INITIALIZER;
    magma_desc_t descA, descB;

    magma = magma_context_self();
    if (magma == NULL) {
        magma_fatal_error("MAGMA_zpotrs", "MAGMA not initialized");
        return MAGMA_ERR_NOT_INITIALIZED;
    }
    /* Check input arguments */
    if (uplo != PlasmaUpper && uplo != PlasmaLower) {
        magma_error("MAGMA_zpotrs", "illegal value of uplo");
        return -1;
    }
    if (N < 0) {
        magma_error("MAGMA_zpotrs", "illegal value of N");
        return -2;
    }
    if (NRHS < 0) {
        magma_error("MAGMA_zpotrs", "illegal value of NRHS");
        return -3;
    }
    if (LDA < max(1, N)) {
        magma_error("MAGMA_zpotrs", "illegal value of LDA");
        return -5;
    }
    if (LDB < max(1, N)) {
        magma_error("MAGMA_zpotrs", "illegal value of LDB");
        return -7;
    }
    /* Quick return */
    if (min(N, NRHS) == 0)
        return MAGMA_SUCCESS;

    /* Tune NB depending on M, N && NRHS; Set NBNB */
    /* status = magma_tune(MAGMA_FUNC_ZPOSV, N, N, NRHS); */
    /* if (status != MAGMA_SUCCESS) { */
    /*     magma_error("MAGMA_zpotrs", "magma_tune() failed"); */
    /*     return status; */
    /* } */

    /* Set NT && NTRHS */
    NB    = MAGMA_NB;

    magma_sequence_create(magma, &sequence);

    /* if ( MAGMA_TRANSLATION == MAGMA_OUTOFPLACE ) { */
        magma_zooplap2tile( descA, A, NB, NB, LDA, N,    0, 0, N, N   , magma_desc_mat_free(&(descA)) );
        magma_zooplap2tile( descB, B, NB, NB, LDB, NRHS, 0, 0, N, NRHS, magma_desc_mat_free(&(descA)); magma_desc_mat_free(&(descB)));
    /* } else { */
    /*     magma_ziplap2tile( descA, A, NB, NB, LDA, N,    0, 0, N, N   ); */
    /*     magma_ziplap2tile( descB, B, NB, NB, LDB, NRHS, 0, 0, N, NRHS); */
    /* } */

    /* Call the tile interface */
    MAGMA_zpotrs_Tile_Async(uplo, &descA, &descB, sequence, &request);

    /* if ( MAGMA_TRANSLATION == MAGMA_OUTOFPLACE ) { */
        magma_zooptile2lap( descB, B, NB, NB, LDB, NRHS );
        morse_barrier( magma );
        morse_desc_getoncpu( &descB );
        magma_desc_mat_free(&descA);
        magma_desc_mat_free(&descB);
    /* } else { */
    /*     magma_ziptile2lap( descA, A, NB, NB, LDA, N    ); */
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
 *  MAGMA_zpotrs_Tile - Solves a system of linear equations using previously
 *  computed Cholesky factorization.
 *  Tile equivalent of MAGMA_zpotrs().
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          = PlasmaUpper: Upper triangle of A is stored;
 *          = PlasmaLower: Lower triangle of A is stored.
 *
 * @param[in] A
 *          The triangular factor U or L from the Cholesky factorization A = U**H*U or A = L*L**H,
 *          computed by MAGMA_zpotrf.
 *
 * @param[in,out] B
 *          On entry, the N-by-NRHS right hand side matrix B.
 *          On exit, if return value = 0, the N-by-NRHS solution matrix X.
 *
 *******************************************************************************
 *
 * @return
 *          \retval MAGMA_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa MAGMA_zpotrs
 * @sa MAGMA_zpotrs_Tile_Async
 * @sa MAGMA_cpotrs_Tile
 * @sa MAGMA_dpotrs_Tile
 * @sa MAGMA_spotrs_Tile
 * @sa MAGMA_zpotrf_Tile
 *
 ******************************************************************************/
int MAGMA_zpotrs_Tile(PLASMA_enum uplo, magma_desc_t *A, magma_desc_t *B)
{
    magma_context_t *magma;
    magma_sequence_t *sequence = NULL;
    magma_request_t request = MAGMA_REQUEST_INITIALIZER;
    int status;

    magma = magma_context_self();
    if (magma == NULL) {
        magma_fatal_error("MAGMA_zpotrs_Tile", "MAGMA not initialized");
        return MAGMA_ERR_NOT_INITIALIZED;
    }
    magma_sequence_create(magma, &sequence);
    MAGMA_zpotrs_Tile_Async(uplo, A, B, sequence, &request);
    morse_barrier( magma );
    morse_desc_getoncpu( A );
    morse_desc_getoncpu( B );
    status = sequence->status;
    magma_sequence_destroy(magma, sequence);
    return status;
}

/***************************************************************************//**
 *
 * @ingroup MAGMA_Complex64_t_Tile_Async
 *
 *  MAGMA_zpotrs_Tile_Async - Solves a system of linear equations using previously
 *  computed Cholesky factorization.
 *  Non-blocking equivalent of MAGMA_zpotrs_Tile().
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
 * @sa MAGMA_zpotrs
 * @sa MAGMA_zpotrs_Tile
 * @sa MAGMA_cpotrs_Tile_Async
 * @sa MAGMA_dpotrs_Tile_Async
 * @sa MAGMA_spotrs_Tile_Async
 * @sa MAGMA_zpotrf_Tile_Async
 *
 ******************************************************************************/
int MAGMA_zpotrs_Tile_Async(PLASMA_enum uplo, magma_desc_t *A, magma_desc_t *B,
                             magma_sequence_t *sequence, magma_request_t *request)
{
    PLASMA_desc descA = A->desc;
    PLASMA_desc descB = B->desc;
    magma_context_t *magma;

    magma = magma_context_self();
    if (magma == NULL) {
        magma_fatal_error("MAGMA_zpotrs_Tile", "MAGMA not initialized");
        return MAGMA_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        magma_fatal_error("MAGMA_zpotrs_Tile", "NULL sequence");
        return MAGMA_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        magma_fatal_error("MAGMA_zpotrs_Tile", "NULL request");
        return MAGMA_ERR_UNALLOCATED;
    }
    /* Check sequence status */
    if (sequence->status == MAGMA_SUCCESS)
        request->status = MAGMA_SUCCESS;
    else
        return magma_request_fail(sequence, request, MAGMA_ERR_SEQUENCE_FLUSHED);

    /* Check descriptors for correctness */
    if (magma_desc_check( A ) != MAGMA_SUCCESS) {
        magma_error("MAGMA_zpotrs_Tile", "invalid first descriptor");
        return magma_request_fail(sequence, request, MAGMA_ERR_ILLEGAL_VALUE);
    }
    if (magma_desc_check( B ) != MAGMA_SUCCESS) {
        magma_error("MAGMA_zpotrs_Tile", "invalid second descriptor");
        return magma_request_fail(sequence, request, MAGMA_ERR_ILLEGAL_VALUE);
    }
    /* Check input arguments */
    if (descA.nb != descA.mb || descB.nb != descB.mb) {
        magma_error("MAGMA_zpotrs_Tile", "only square tiles supported");
        return magma_request_fail(sequence, request, MAGMA_ERR_ILLEGAL_VALUE);
    }
    if (uplo != PlasmaUpper && uplo != PlasmaLower) {
        magma_error("MAGMA_zpotrs_Tile", "illegal value of uplo");
        return magma_request_fail(sequence, request, -1);
    }
    /* Quick return */
/*
    if (min(N, NRHS) == 0)
        return MAGMA_SUCCESS;
*/

    magma_pztrsm( PlasmaLeft, uplo,
                  uplo == PlasmaUpper ? PlasmaConjTrans : PlasmaNoTrans,
                  PlasmaNonUnit,
                  1.0, A, B,
                  sequence, request);

    magma_pztrsm( PlasmaLeft, uplo,
                  uplo == PlasmaUpper ? PlasmaNoTrans : PlasmaConjTrans,
                  PlasmaNonUnit,
                  1.0, A, B,
                  sequence, request);

    return PLASMA_SUCCESS;
}
