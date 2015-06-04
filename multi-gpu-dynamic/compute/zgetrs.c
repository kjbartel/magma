/**
 *
 * @file zgetrs.c
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
 *  MAGMA_zgetrs - Solves a system of linear equations A * X = B, with a general N-by-N matrix A
 *  using the tile LU factorization computed by MAGMA_zgetrf.
 *
 *******************************************************************************
 *
 * @param[in] trans
 *          Intended to specify the the form of the system of equations:
 *          = PlasmaNoTrans:   A * X = B     (No transpose)
 *          = PlasmaTrans:     A**T * X = B  (Transpose)
 *          = PlasmaConjTrans: A**H * X = B  (Conjugate transpose)
 *          Currently only PlasmaNoTrans is supported.
 *
 * @param[in] N
 *          The order of the matrix A.  N >= 0.
 *
 * @param[in] NRHS
 *          The number of right hand sides, i.e., the number of columns of the matrix B.
 *          NRHS >= 0.
 *
 * @param[in] A
 *          The tile factors L and U from the factorization, computed by MAGMA_zgetrf.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,N).
 *
 * @param[in] L
 *          Auxiliary factorization data, related to the tile L factor, computed by MAGMA_zgetrf.
 *
 * @param[in] IPIV
 *          The pivot indices from MAGMA_zgetrf (not equivalent to LAPACK).
 *
 * @param[in,out] B
 *          On entry, the N-by-NRHS matrix of right hand side matrix B.
 *          On exit, the solution matrix X.
 *
 * @param[in] LDB
 *          The leading dimension of the array B. LDB >= max(1,N).
 *
 *******************************************************************************
 *
 * @return
 *          \retval MAGMA_SUCCESS successful exit
 *          \return <0 if -i, the i-th argument had an illegal value
 *
 *******************************************************************************
 *
 * @sa MAGMA_zgetrs_Tile
 * @sa MAGMA_zgetrs_Tile_Async
 * @sa MAGMA_cgetrs
 * @sa MAGMA_dgetrs
 * @sa MAGMA_sgetrs
 * @sa MAGMA_zgetrf
 *
 ******************************************************************************/
int MAGMA_zgetrs(PLASMA_enum trans, int N, int NRHS,
                  PLASMA_Complex64_t *A, int LDA,
                  magma_desc_t *L, int *IPIV,
                  PLASMA_Complex64_t *B, int LDB)
{
    int NB, IB, IBNB, NT;
    int status;
    magma_context_t *magma;
    magma_sequence_t *sequence = NULL;
    magma_request_t request = MAGMA_REQUEST_INITIALIZER;
    magma_desc_t descA, descB;

    magma = magma_context_self();
    if (magma == NULL) {
        magma_fatal_error("MAGMA_zgetrs", "MAGMA not initialized");
        return MAGMA_ERR_NOT_INITIALIZED;
    }
    /* Check input arguments */
    if (trans != PlasmaNoTrans) {
        magma_error("MAGMA_zgetrs", "only PlasmaNoTrans supported");
        return MAGMA_ERR_NOT_SUPPORTED;
    }
    if (N < 0) {
        magma_error("MAGMA_zgetrs", "illegal value of N");
        return -2;
    }
    if (NRHS < 0) {
        magma_error("MAGMA_zgetrs", "illegal value of NRHS");
        return -3;
    }
    if (LDA < max(1, N)) {
        magma_error("MAGMA_zgetrs", "illegal value of LDA");
        return -5;
    }
    if (LDB < max(1, N)) {
        magma_error("MAGMA_zgetrs", "illegal value of LDB");
        return -9;
    }
    /* Quick return */
    if (min(N, NRHS) == 0)
        return MAGMA_SUCCESS;

    /* Tune NB && IB depending on N && NRHS; Set NBNBSIZE */
    /* status = magma_tune(PLASMA_FUNC_ZGESV, N, N, NRHS); */
    /* if (status != PLASMA_SUCCESS) { */
    /*     plasma_error("PLASMA_zgetrs", "plasma_tune() failed"); */
    /*     return status; */
    /* } */

    /* Set NT && NTRHS */
    NB    = MAGMA_NB;
    IB    = MAGMA_IB;
    IBNB  = IB*NB;
    NT    = (N%NB==0) ? (N/NB) : (N/NB+1);

    magma_sequence_create(magma, &sequence);

    /* if ( MAGMA_TRANSLATION == MAGMA_OUTOFPLACE ) { */
        magma_zooplap2tile( descA, A, NB, NB, LDA, N,    0, 0, N, N   , magma_desc_mat_free(&(descA)) );
        magma_zooplap2tile( descB, B, NB, NB, LDB, NRHS, 0, 0, N, NRHS, magma_desc_mat_free(&(descA)); magma_desc_mat_free(&(descB)));
    /* } else { */
    /*     magma_ziplap2tile( descA, A, NB, NB, LDA, N,    0, 0, N, N   ); */
    /*     magma_ziplap2tile( descB, B, NB, NB, LDB, NRHS, 0, 0, N, NRHS); */
    /* } */

    /* Call the tile interface */
    MAGMA_zgetrs_Tile_Async(&descA, L, IPIV, &descB, sequence, &request);

    /* if ( MAGMA_TRANSLATION == MAGMA_OUTOFPLACE ) { */
        magma_zooptile2lap( descB, B, NB, NB, LDB, NRHS );
        morse_barrier( magma );
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
 *  MAGMA_zgetrs_Tile - Solves a system of linear equations using previously
 *  computed LU factorization.
 *  Tile equivalent of MAGMA_zgetrs().
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] A
 *          The tile factors L and U from the factorization, computed by MAGMA_zgetrf.
 *
 * @param[in] L
 *          Auxiliary factorization data, related to the tile L factor, computed by MAGMA_zgetrf.
 *
 * @param[in] IPIV
 *          The pivot indices from MAGMA_zgetrf (not equivalent to LAPACK).
 *
 * @param[in,out] B
 *          On entry, the N-by-NRHS matrix of right hand side matrix B.
 *          On exit, the solution matrix X.
 *
 *******************************************************************************
 *
 * @return
 *          \retval MAGMA_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa MAGMA_zgetrs
 * @sa MAGMA_zgetrs_Tile_Async
 * @sa MAGMA_cgetrs_Tile
 * @sa MAGMA_dgetrs_Tile
 * @sa MAGMA_sgetrs_Tile
 * @sa MAGMA_zgetrf_Tile
 *
 ******************************************************************************/
int MAGMA_zgetrs_Tile(magma_desc_t *A, magma_desc_t *L, int *IPIV, magma_desc_t *B)
{
    magma_context_t *magma;
    magma_sequence_t *sequence = NULL;
    magma_request_t request = MAGMA_REQUEST_INITIALIZER;
    int status;

    magma = magma_context_self();
    if (magma == NULL) {
        magma_fatal_error("MAGMA_zgetrs_Tile", "MAGMA not initialized");
        return MAGMA_ERR_NOT_INITIALIZED;
    }
    magma_sequence_create(magma, &sequence);
    MAGMA_zgetrs_Tile_Async(A, L, IPIV, B, sequence, &request);
    morse_barrier( magma );
    morse_desc_getoncpu( A );
    morse_desc_getoncpu( L );
    morse_desc_getoncpu( B );
    status = sequence->status;
    magma_sequence_destroy(magma, sequence);
    return status;
}

/***************************************************************************//**
 *
 * @ingroup MAGMA_Complex64_t_Tile_Async
 *
 *  MAGMA_zgetrs_Tile_Async - Solves a system of linear equations using previously
 *  computed LU factorization.
 *  Non-blocking equivalent of MAGMA_zgetrs_Tile().
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
 * @sa MAGMA_zgetrs
 * @sa MAGMA_zgetrs_Tile
 * @sa MAGMA_cgetrs_Tile_Async
 * @sa MAGMA_dgetrs_Tile_Async
 * @sa MAGMA_sgetrs_Tile_Async
 * @sa MAGMA_zgetrf_Tile_Async
 *
 ******************************************************************************/
int MAGMA_zgetrs_Tile_Async(magma_desc_t *A, magma_desc_t *L, int *IPIV, magma_desc_t *B,
                             magma_sequence_t *sequence, magma_request_t *request)
{
    PLASMA_desc descA = A->desc;
    PLASMA_desc descB = B->desc;
    magma_context_t *magma;

    magma = magma_context_self();
    if (magma == NULL) {
        magma_fatal_error("MAGMA_zgetrs_Tile", "MAGMA not initialized");
        return MAGMA_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        magma_fatal_error("MAGMA_zgetrs_Tile", "NULL sequence");
        return MAGMA_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        magma_fatal_error("MAGMA_zgetrs_Tile", "NULL request");
        return MAGMA_ERR_UNALLOCATED;
    }
    /* Check sequence status */
    if (sequence->status == MAGMA_SUCCESS)
        request->status = MAGMA_SUCCESS;
    else
        return magma_request_fail(sequence, request, MAGMA_ERR_SEQUENCE_FLUSHED);

    /* Check descriptors for correctness */
    if (magma_desc_check( A ) != MAGMA_SUCCESS) {
        magma_error("MAGMA_zgetrs_Tile", "invalid first descriptor");
        return magma_request_fail(sequence, request, MAGMA_ERR_ILLEGAL_VALUE);
    }
    if (magma_desc_check( L ) != MAGMA_SUCCESS) {
        magma_error("MAGMA_zgetrs_Tile", "invalid second descriptor");
        return magma_request_fail(sequence, request, MAGMA_ERR_ILLEGAL_VALUE);
    }
    if (magma_desc_check( B ) != MAGMA_SUCCESS) {
        magma_error("MAGMA_zgetrs_Tile", "invalid third descriptor");
        return magma_request_fail(sequence, request, MAGMA_ERR_ILLEGAL_VALUE);
    }
    /* Check input arguments */
    if (descA.nb != descA.mb || descB.nb != descB.mb) {
        magma_error("MAGMA_zgetrs_Tile", "only square tiles supported");
        return magma_request_fail(sequence, request, MAGMA_ERR_ILLEGAL_VALUE);
    }
    /* Quick return */
/*
    if (min(N, NRHS) == 0)
        return PLASMA_SUCCESS;
*/
    magma_pztrsmpl(A, B, L, IPIV, sequence, request);
    magma_pztrsm( PlasmaLeft, PlasmaUpper, PlasmaNoTrans, PlasmaNonUnit,
                  1.0, A, B,
                  sequence, request);

    return PLASMA_SUCCESS;
}
