/**
 *
 * @file zgetrf_incpiv.c
 *
 *  PLASMA computational routines
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

void plasma_memzero(void *memptr, PLASMA_size size, int type);

/***************************************************************************//**
 *
 * @ingroup MAGMA_Complex64_t
 *
 *  MAGMA_zgetrf_incpiv - Computes an LU factorization of a general M-by-N matrix A
 *  using the tile LU algorithm with partial tile pivoting with row interchanges.
 *
 *******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the matrix A. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix A. N >= 0.
 *
 * @param[in,out] A
 *          On entry, the M-by-N matrix to be factored.
 *          On exit, the tile factors L and U from the factorization.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,M).
 *
 * @param[out] L
 *          On exit, auxiliary factorization data, related to the tile L factor,
 *          required by MAGMA_zgetrs to solve the system of equations.
 *
 * @param[out] IPIV
 *          The pivot indices that define the permutations (not equivalent to LAPACK).
 *
 *******************************************************************************
 *
 * @return
 *          \retval MAGMA_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *          \retval >0 if i, U(i,i) is exactly zero. The factorization has been completed,
 *               but the factor U is exactly singular, and division by zero will occur
 *               if it is used to solve a system of equations.
 *
 *******************************************************************************
 *
 * @sa MAGMA_zgetrf_incpiv_Tile
 * @sa MAGMA_zgetrf_incpiv_Tile_Async
 * @sa MAGMA_cgetrf_incpiv
 * @sa MAGMA_dgetrf_incpiv
 * @sa MAGMA_sgetrf_incpiv
 * @sa MAGMA_zgetrs_incpiv
 *
 ******************************************************************************/
int MAGMA_zgetrf_incpiv(int M, int N,
                  PLASMA_Complex64_t *A, int LDA,
                  magma_desc_t *L, int *IPIV)
{
    int NB, IB, IBNB, MT, NT;
    int status;
    magma_context_t *magma;
    magma_sequence_t *sequence = NULL;
    magma_request_t request = MAGMA_REQUEST_INITIALIZER;
    magma_desc_t descA;

    magma = magma_context_self();
    if (magma == NULL) {
        magma_fatal_error("MAGMA_zgetrf_incpiv", "MAGMA not initialized");
        return MAGMA_ERR_NOT_INITIALIZED;
    }
    /* Check input arguments */
    if (M < 0) {
        magma_error("MAGMA_zgetrf_incpiv", "illegal value of M");
        return -1;
    }
    if (N < 0) {
        magma_error("MAGMA_zgetrf_incpiv", "illegal value of N");
        return -2;
    }
    if (LDA < max(1, M)) {
        magma_error("MAGMA_zgetrf_incpiv", "illegal value of LDA");
        return -4;
    }
    /* Quick return */
    if (min(M, N) == 0)
        return MAGMA_SUCCESS;

    /* Tune NB && IB depending on M, N && NRHS; Set NBNBSIZE */
    /* status = plasma_tune(PLASMA_FUNC_ZGESV, M, N, 0); */
    /* if (status != PLASMA_SUCCESS) { */
    /*     plasma_error("PLASMA_zgetrf_incpiv", "plasma_tune() failed"); */
    /*     return status; */
    /* } */

    /* Set NT && NTRHS */
    NB   = MAGMA_NB;
    IB   = MAGMA_IB;
    IBNB = IB*NB;
    MT   = (M%NB==0) ? (M/NB) : (M/NB+1);
    NT   = (N%NB==0) ? (N/NB) : (N/NB+1);

    magma_sequence_create(magma, &sequence);

    /* if ( MAGMA_TRANSLATION == MAGMA_OUTOFPLACE ) { */
        magma_zooplap2tile( descA, A, NB, NB, LDA, N, 0, 0, M, N, magma_desc_mat_free(&(descA)) );
    /* } else { */
    /*     magma_ziplap2tile( descA, A, NB, NB, LDA, N, 0, 0, M, N); */
    /* } */

    /* Call the tile interface */
    MAGMA_zgetrf_incpiv_Tile_Async(&descA, L, IPIV, sequence, &request);

    /* if ( MAGMA_TRANSLATION == MAGMA_OUTOFPLACE ) { */
        magma_zooptile2lap( descA, A, NB, NB, LDA, N );
        morse_barrier( magma );
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
 *  MAGMA_zgetrf_incpiv_Tile - Computes the tile LU factorization of a matrix.
 *  Tile equivalent of MAGMA_zgetrf_incpiv().
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in,out] A
 *          On entry, the M-by-N matrix to be factored.
 *          On exit, the tile factors L and U from the factorization.
 *
 * @param[out] L
 *          On exit, auxiliary factorization data, related to the tile L factor,
 *          required by MAGMA_zgetrs to solve the system of equations.
 *
 * @param[out] IPIV
 *          The pivot indices that define the permutations (not equivalent to LAPACK).
 *
 *******************************************************************************
 *
 * @return
 *          \retval MAGMA_SUCCESS successful exit
 *          \retval >0 if i, U(i,i) is exactly zero. The factorization has been completed,
 *               but the factor U is exactly singular, and division by zero will occur
 *               if it is used to solve a system of equations.
 *
 *******************************************************************************
 *
 * @sa MAGMA_zgetrf_incpiv
 * @sa MAGMA_zgetrf_incpiv_Tile_Async
 * @sa MAGMA_cgetrf_incpiv_Tile
 * @sa MAGMA_dgetrf_incpiv_Tile
 * @sa MAGMA_sgetrf_incpiv_Tile
 * @sa MAGMA_zgetrs_incpiv_Tile
 *
 ******************************************************************************/
int MAGMA_zgetrf_incpiv_Tile(magma_desc_t *A, magma_desc_t *L, int *IPIV)
{
    magma_context_t *magma;
    magma_sequence_t *sequence = NULL;
    magma_request_t request = MAGMA_REQUEST_INITIALIZER;
    int status;

    magma = magma_context_self();
    if (magma == NULL) {
        magma_fatal_error("MAGMA_zgetrf_incpiv_Tile", "MAGMA not initialized");
        return MAGMA_ERR_NOT_INITIALIZED;
    }
    magma_sequence_create(magma, &sequence);
    MAGMA_zgetrf_incpiv_Tile_Async(A, L, IPIV, sequence, &request);
    morse_barrier( magma );
    morse_desc_getoncpu( A );
    morse_desc_getoncpu( L );
    status = sequence->status;
    magma_sequence_destroy(magma, sequence);
    return status;
}

/***************************************************************************//**
 *
 * @ingroup MAGMA_Complex64_t_Tile_Async
 *
 *  MAGMA_zgetrf_incpiv_Tile_Async - Computes the tile LU factorization of a matrix.
 *  Non-blocking equivalent of MAGMA_zgetrf_incpiv_Tile().
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
 * @sa MAGMA_zgetrf_incpiv
 * @sa MAGMA_zgetrf_incpiv_Tile
 * @sa MAGMA_cgetrf_incpiv_Tile_Async
 * @sa MAGMA_dgetrf_incpiv_Tile_Async
 * @sa MAGMA_sgetrf_incpiv_Tile_Async
 * @sa MAGMA_zgetrs_incpiv_Tile_Async
 *
 ******************************************************************************/
int MAGMA_zgetrf_incpiv_Tile_Async(magma_desc_t *A, magma_desc_t *L, int *IPIV,
                              magma_sequence_t *sequence, magma_request_t *request)
{
    magma_context_t *magma;

    magma = magma_context_self();
    if (magma == NULL) {
        magma_fatal_error("MAGMA_zgetrf_incpiv_Tile", "MAGMA not initialized");
        return MAGMA_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        magma_fatal_error("MAGMA_zgetrf_incpiv_Tile", "NULL sequence");
        return MAGMA_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        magma_fatal_error("MAGMA_zgetrf_incpiv_Tile", "NULL request");
        return MAGMA_ERR_UNALLOCATED;
    }
    /* Check sequence status */
    if (sequence->status == MAGMA_SUCCESS)
        request->status = MAGMA_SUCCESS;
    else
        return magma_request_fail(sequence, request, MAGMA_ERR_SEQUENCE_FLUSHED);

    /* Check descriptors for correctness */
    if (magma_desc_check( A ) != MAGMA_SUCCESS) {
        magma_error("MAGMA_zgetrf_incpiv_Tile", "invalid first descriptor");
        return magma_request_fail(sequence, request, MAGMA_ERR_ILLEGAL_VALUE);
    }
    if (magma_desc_check( L ) != MAGMA_SUCCESS) {
        magma_error("MAGMA_zgetrf_incpiv_Tile", "invalid second descriptor");
        return magma_request_fail(sequence, request, MAGMA_ERR_ILLEGAL_VALUE);
    }
    /* Check input arguments */
    if (A->nb != A->mb) {
        magma_error("MAGMA_zgetrf_incpiv_Tile", "only square tiles supported");
        return magma_request_fail(sequence, request, MAGMA_ERR_ILLEGAL_VALUE);
    }
    /* Quick return */
/*
    if (min(M, N) == 0)
        return MAGMA_SUCCESS;
*/
    /* Clear IPIV and Lbdl */
    plasma_memzero(IPIV,      A->mt*A->nt*A->nb, PlasmaInteger);
    plasma_memzero(L->mat, L->lm*L->ln,          PlasmaComplexDouble);

    magma_pzgetrf_incpiv( A, L, IPIV, sequence, request);

    return PLASMA_SUCCESS;
}
