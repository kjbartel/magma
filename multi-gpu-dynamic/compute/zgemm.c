/**
 *
 * @file zgemm.c
 *
 *  MAGMA computational routines
 *  MAGMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.3.1
 * @author Emmanuel Agullo
 * @date 2010-11-15
 * @precisions normal z -> s d c
 *
 **/
#include "common.h"

/***************************************************************************//**
 *
 * @ingroup PLASMA_Complex64_t
 *
 *  MAGMA_zgemm - Performs one of the matrix-matrix operations
 *
 *    \f[ C = \alpha [op( A )\times op( B )] + \beta C \f],
 *
 *  where op( X ) is one of
 *
 *    op( X ) = X  or op( X ) = X' or op( X ) = conjg( X' )
 *
 *  alpha and beta are scalars, and A, B and C  are matrices, with op( A )
 *  an m by k matrix, op( B ) a k by n matrix and C an m by n matrix.
 *
 *******************************************************************************
 *
 * @param[in] transA
 *          Specifies whether the matrix A is transposed, not transposed or conjugate transposed:
 *          = PlasmaNoTrans:   A is not transposed;
 *          = PlasmaTrans:     A is transposed;
 *          = PlasmaConjTrans: A is conjugate transposed.
 *
 * @param[in] transB
 *          Specifies whether the matrix B is transposed, not transposed or conjugate transposed:
 *          = PlasmaNoTrans:   B is not transposed;
 *          = PlasmaTrans:     B is transposed;
 *          = PlasmaConjTrans: B is conjugate transposed.
 *
 * @param[in] M
 *          M specifies the number of rows of the matrix op( A ) and of the matrix C. M >= 0.
 *
 * @param[in] N
 *          N specifies the number of columns of the matrix op( B ) and of the matrix C. N >= 0.
 *
 * @param[in] K
 *          K specifies the number of columns of the matrix op( A ) and the number of rows of
 *          the matrix op( B ). K >= 0.
 *
 * @param[in] alpha
 *          alpha specifies the scalar alpha
 *
 * @param[in] A
 *          A is a LDA-by-ka matrix, where ka is K when  transA = PlasmaNoTrans,
 *          and is  M  otherwise.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,M).
 *
 * @param[in] B
 *          B is a LDB-by-kb matrix, where kb is N when  transB = PlasmaNoTrans,
 *          and is  K  otherwise.
 *
 * @param[in] LDB
 *          The leading dimension of the array B. LDB >= max(1,N).
 *
 * @param[in] beta
 *          beta specifies the scalar beta
 *
 * @param[in,out] C
 *          C is a LDC-by-N matrix.
 *          On exit, the array is overwritten by the M by N matrix ( alpha*op( A )*op( B ) + beta*C )
 *
 * @param[in] LDC
 *          The leading dimension of the array C. LDC >= max(1,M).
 *
 *******************************************************************************
 *
 * @return
 *          \retval MAGMA_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa MAGMA_zgemm_Tile
 * @sa MAGMA_cgemm
 * @sa MAGMA_dgemm
 * @sa MAGMA_sgemm
 *
 ******************************************************************************/
int MAGMA_zgemm(PLASMA_enum transA, PLASMA_enum transB, int M, int N, int K,
                 PLASMA_Complex64_t alpha, PLASMA_Complex64_t *A, int LDA,
                                           PLASMA_Complex64_t *B, int LDB,
                 PLASMA_Complex64_t beta,  PLASMA_Complex64_t *C, int LDC)
{
    int NB;
    int Am, An, Bm, Bn;
    int status;
    magma_desc_t descA, descB, descC;
    magma_context_t *magma;
    magma_sequence_t *sequence = NULL;
    magma_request_t request = MAGMA_REQUEST_INITIALIZER;

    magma = magma_context_self();
    if (magma == NULL) {
        magma_fatal_error("MAGMA_zgemm", "MAGMA not initialized");
        return MAGMA_ERR_NOT_INITIALIZED;
    }

    /* Check input arguments */
    if ((transA != PlasmaNoTrans) && (transA != PlasmaTrans) && (transA != PlasmaConjTrans)) {
        magma_error("MAGMA_zgemm", "illegal value of transA");
        return -1;
    }
    if ((transB != PlasmaNoTrans) && (transB != PlasmaTrans) && (transB != PlasmaConjTrans)) {
        magma_error("MAGMA_zgemm", "illegal value of transB");
        return -2;
    }
    if ( transA == PlasmaNoTrans ) { 
        Am = M; An = K;
    } else {
        Am = K; An = M;
    }
    if ( transB == PlasmaNoTrans ) { 
        Bm = K; Bn = N;
    } else {
        Bm = N; Bn = K;
    }
    if (M < 0) {
        magma_error("MAGMA_zgemm", "illegal value of M");
        return -3;
    }
    if (N < 0) {
        magma_error("MAGMA_zgemm", "illegal value of N");
        return -4;
    }
    if (K < 0) {
        magma_error("MAGMA_zgemm", "illegal value of N");
        return -5;
    }
    if (LDA < max(1, Am)) {
        magma_error("MAGMA_zgemm", "illegal value of LDA");
        return -8;
    }
    if (LDB < max(1, Bm)) {
        magma_error("MAGMA_zgemm", "illegal value of LDB");
        return -10;
    }
    if (LDC < max(1, M)) {
        magma_error("MAGMA_zgemm", "illegal value of LDC");
        return -13;
    }

    /* Quick return */
    if (M == 0 || N == 0 ||
        ((alpha == (PLASMA_Complex64_t)0.0 || K == 0) && beta == (PLASMA_Complex64_t)1.0))
        return MAGMA_SUCCESS;

    /* Tune NB depending on M, N && NRHS; Set NBNBSIZE */
    /* status = magma_tune(MAGMA_FUNC_ZGEMM, M, N, 0); */
    /* if (status != MAGMA_SUCCESS) { */
    /*     magma_error("MAGMA_zgemm", "magma_tune() failed"); */
    /*     return status; */
    /* } */

    /* Set MT && NT && KT */
    NB = MAGMA_NB;

    magma_sequence_create(magma, &sequence);

    /* if ( MAGMA_TRANSLATION == MAGMA_OUTOFPLACE ) { */
        magma_zooplap2tile( descA, A, NB, NB, LDA, An, 0, 0, Am, An, magma_desc_mat_free(&(descA)) );
        magma_zooplap2tile( descB, B, NB, NB, LDB, Bn, 0, 0, Bm, Bn, magma_desc_mat_free(&(descA)); magma_desc_mat_free(&(descB)));
        magma_zooplap2tile( descC, C, NB, NB, LDC, N,  0, 0, M,  N,  magma_desc_mat_free(&(descA)); magma_desc_mat_free(&(descB)); magma_desc_mat_free(&(descC)));
    /* } else { */
    /*     magma_ziplap2tile( descA, A, NB, NB, LDA, An, 0, 0, Am, An ); */
    /*     magma_ziplap2tile( descB, B, NB, NB, LDB, Bn, 0, 0, Bm, Bn ); */
    /*     magma_ziplap2tile( descC, C, NB, NB, LDC, N,  0, 0, M,  N  ); */
    /* } */

    /* Call the tile interface */
    MAGMA_zgemm_Tile_Async(
        transA, transB, alpha, &descA, &descB, beta, &descC, sequence, &request);

    /* if ( MAGMA_TRANSLATION == MAGMA_OUTOFPLACE ) { */
        magma_zooptile2lap( descC, C, NB, NB, LDC, N );
        morse_barrier( magma );
        morse_desc_getoncpu( &descC );
        magma_desc_mat_free(&descA);
        magma_desc_mat_free(&descB);
        magma_desc_mat_free(&descC);
    /* } else { */
    /*     magma_ziptile2lap( descA, A, NB, NB, LDA, An ); */
    /*     magma_ziptile2lap( descB, B, NB, NB, LDB, Bn ); */
    /*     magma_ziptile2lap( descC, C, NB, NB, LDC, N  ); */
    /*     morse_barrier( magma ); */
    /* } */

    status = sequence->status;
    magma_sequence_destroy(magma, sequence);
    return status;
}

/***************************************************************************//**
 *
 * @ingroup PLASMA_Complex64_t_Tile
 *
 *  MAGMA_zgemm_Tile - Performs matrix multiplication.
 *  Tile equivalent of MAGMA_zgemm().
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] transA
 *          Specifies whether the matrix A is transposed, not transposed or conjugate transposed:
 *          = PlasmaNoTrans:   A is not transposed;
 *          = PlasmaTrans:     A is transposed;
 *          = PlasmaConjTrans: A is conjugate transposed.
 *
 * @param[in] transB
 *          Specifies whether the matrix B is transposed, not transposed or conjugate transposed:
 *          = PlasmaNoTrans:   B is not transposed;
 *          = PlasmaTrans:     B is transposed;
 *          = PlasmaConjTrans: B is conjugate transposed.
 *
 * @param[in] alpha
 *          alpha specifies the scalar alpha
 *
 * @param[in] A
 *          A is a LDA-by-ka matrix, where ka is K when  transA = PlasmaNoTrans,
 *          and is  M  otherwise.
 *
 * @param[in] B
 *          B is a LDB-by-kb matrix, where kb is N when  transB = PlasmaNoTrans,
 *          and is  K  otherwise.
 *
 * @param[in] beta
 *          beta specifies the scalar beta
 *
 * @param[in,out] C
 *          C is a LDC-by-N matrix.
 *          On exit, the array is overwritten by the M by N matrix ( alpha*op( A )*op( B ) + beta*C )
 *
 *******************************************************************************
 *
 * @return
 *          \retval MAGMA_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa MAGMA_zgemm
 * @sa MAGMA_zgemm_Tile_Async
 * @sa MAGMA_cgemm_Tile
 * @sa MAGMA_dgemm_Tile
 * @sa MAGMA_sgemm_Tile
 *
 ******************************************************************************/
int MAGMA_zgemm_Tile(PLASMA_enum transA, PLASMA_enum transB,
                       PLASMA_Complex64_t alpha, magma_desc_t *A, magma_desc_t *B,
                       PLASMA_Complex64_t beta,  magma_desc_t *C)
{
    magma_context_t *magma;
    magma_sequence_t *sequence = NULL;
    magma_request_t request = MAGMA_REQUEST_INITIALIZER;
    int status;

    magma = magma_context_self();
    if (magma == NULL) {
        magma_fatal_error("MAGMA_zgemm_Tile", "MAGMA not initialized");
        return MAGMA_ERR_NOT_INITIALIZED;
    }
    magma_sequence_create(magma, &sequence);
    MAGMA_zgemm_Tile_Async(transA, transB, alpha, A, B, beta, C, sequence, &request);
    morse_barrier( magma );
    morse_desc_getoncpu( A );
    morse_desc_getoncpu( B );
    morse_desc_getoncpu( C );
    status = sequence->status;
    magma_sequence_destroy(magma, sequence);
    return status;
}

/***************************************************************************//**
 *
 * @ingroup PLASMA_Complex64_t_Tile_Async
 *
 *  MAGMA_zgemm_Tile_Async - Performs matrix multiplication.
 *  Non-blocking equivalent of MAGMA_zgemm_Tile().
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
 * @sa MAGMA_zgemm
 * @sa MAGMA_zgemm_Tile
 * @sa MAGMA_cgemm_Tile_Async
 * @sa MAGMA_dgemm_Tile_Async
 * @sa MAGMA_sgemm_Tile_Async
 *
 ******************************************************************************/
int MAGMA_zgemm_Tile_Async(PLASMA_enum transA, PLASMA_enum transB,
                             PLASMA_Complex64_t alpha, magma_desc_t *A, magma_desc_t *B,
                             PLASMA_Complex64_t beta,  magma_desc_t *C,
                             magma_sequence_t *sequence, magma_request_t *request)
{
    magma_context_t *magma;
    int M, N, K;
    int Am, An, Ai, Aj, Amb, Anb;
    int Bm, Bn, Bi, Bj, Bmb, Bnb;

    magma = magma_context_self();
    if (magma == NULL) {
        magma_fatal_error("MAGMA_zgemm_Tile_Async", "MAGMA not initialized");
        return MAGMA_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        magma_fatal_error("MAGMA_zgemm_Tile_Async", "NULL sequence");
        return MAGMA_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        magma_fatal_error("MAGMA_zgemm_Tile_Async", "NULL request");
        return MAGMA_ERR_UNALLOCATED;
    }
    /* Check sequence status */
    if (sequence->status == MAGMA_SUCCESS)
        request->status = MAGMA_SUCCESS;
    else
        return magma_request_fail(sequence, request, MAGMA_ERR_SEQUENCE_FLUSHED);

    /* Check descriptors for correctness */
    if (magma_desc_check( A ) != MAGMA_SUCCESS) {
        magma_error("MAGMA_zgemm_Tile_Async", "invalid first descriptor");
        return magma_request_fail(sequence, request, MAGMA_ERR_ILLEGAL_VALUE);
    }
    if (magma_desc_check( B ) != MAGMA_SUCCESS) {
        magma_error("MAGMA_zgemm_Tile_Async", "invalid second descriptor");
        return magma_request_fail(sequence, request, MAGMA_ERR_ILLEGAL_VALUE);
    }
    if (magma_desc_check( C ) != MAGMA_SUCCESS) {
        magma_error("MAGMA_zgemm_Tile_Async", "invalid third descriptor");
        return magma_request_fail(sequence, request, MAGMA_ERR_ILLEGAL_VALUE);
    }
    /* Check input arguments */
    if ((transA != PlasmaNoTrans) && (transA != PlasmaTrans) && (transA != PlasmaConjTrans)) {
        magma_error("MAGMA_zgemm_Tile_Async", "illegal value of transA");
        return magma_request_fail(sequence, request, -1);
    }
    if ((transB != PlasmaNoTrans) && (transB != PlasmaTrans) && (transB != PlasmaConjTrans)) {
        magma_error("MAGMA_zgemm_Tile_Async", "illegal value of transB");
        return magma_request_fail(sequence, request, -2);
    }

    if ( transA == PlasmaNoTrans ) {
        Am  = A->m;
        An  = A->n;
        Amb = A->mb;
        Anb = A->nb;
        Ai  = A->i;
        Aj  = A->j;
    } else {
        Am  = A->n;
        An  = A->m;
        Amb = A->nb;
        Anb = A->mb;
        Ai  = A->j;
        Aj  = A->i;
    }

    if ( transB == PlasmaNoTrans ) {
        Bm  = B->m;
        Bn  = B->n;
        Bmb = B->mb;
        Bnb = B->nb;
        Bi  = B->i;
        Bj  = B->j;
    } else {
        Bm  = B->n;
        Bn  = B->m;
        Bmb = B->nb;
        Bnb = B->mb;
        Bi  = B->j;
        Bj  = B->i;
    }

    if ( (Amb != C->mb) || (Anb != Bmb) || (Bnb != C->nb) ) {
        magma_error("MAGMA_zgemm_Tile_Async", "tile sizes have to match");
        return magma_request_fail(sequence, request, MAGMA_ERR_ILLEGAL_VALUE);
    }
    if ( (Am != C->m) || (An != Bm) || (Bn != C->n) ) {
        magma_error("MAGMA_zgemm_Tile_Async", "sizes of matrices have to match");
        return magma_request_fail(sequence, request, MAGMA_ERR_ILLEGAL_VALUE);
    }
    if ( (Ai != C->i) || (Aj != Bi) || (Bj != C->j) ) {
        magma_error("MAGMA_zgemm_Tile_Async", "start indexes have to match");
        return magma_request_fail(sequence, request, MAGMA_ERR_ILLEGAL_VALUE);
    }

    M = C->m;
    N = C->n;
    K = An;

    /* Quick return */
    if (M == 0 || N == 0 ||
        ((alpha == (PLASMA_Complex64_t)0.0 || K == 0) && beta == (PLASMA_Complex64_t)1.0))
        return MAGMA_SUCCESS;

    magma_pzgemm( transA, transB,
                  alpha, A, B,
                  beta, C,
                  sequence, request);

    return MAGMA_SUCCESS;
}
