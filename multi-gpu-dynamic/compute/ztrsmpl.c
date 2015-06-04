/**
 *
 * @file ztrsmpl.c
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
 *  MAGMA_ztrsmpl - Performs the forward substitution step of solving a system of linear equations
 *  after the tile LU factorization of the matrix.
 *
 *******************************************************************************
 *
 * @param[in] N
 *          The order of the matrix A. N >= 0.
 *
 * @param[in] NRHS
 *          The number of right hand sides, i.e., the number of columns of the matrix B. NRHS >= 0.
 *
 * @param[in] A
 *          The tile factor L from the factorization, computed by MAGMA_zgetrf.
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
 * @sa MAGMA_ztrsmpl_Tile
 * @sa MAGMA_ztrsmpl_Tile_Async
 * @sa MAGMA_ctrsmpl
 * @sa MAGMA_dtrsmpl
 * @sa MAGMA_strsmpl
 * @sa MAGMA_zgetrf
 *
 ******************************************************************************/
int MAGMA_ztrsmpl(int N, int NRHS,
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
        magma_fatal_error("MAGMA_ztrsmpl", "MAGMA not initialized");
        return MAGMA_ERR_NOT_INITIALIZED;
    }
    /* Check input arguments */
    if (N < 0) {
        magma_error("MAGMA_ztrsmpl", "illegal value of N");
        return -1;
    }
    if (NRHS < 0) {
        magma_error("MAGMA_ztrsmpl", "illegal value of NRHS");
        return -2;
    }
    if (LDA < max(1, N)) {
        magma_error("MAGMA_ztrsmpl", "illegal value of LDA");
        return -4;
    }
    if (LDB < max(1, N)) {
        magma_error("MAGMA_ztrsmpl", "illegal value of LDB");
        return -8;
    }
    /* Quick return */
    if (min(N, NRHS) == 0)
        return MAGMA_SUCCESS;

    /* Tune NB && IB depending on N && NRHS; Set NBNB */
    /* status = magma_tune(MAGMA_FUNC_ZGESV, N, N, NRHS); */
    /* if (status != MAGMA_SUCCESS) { */
    /*     magma_error("MAGMA_ztrsmpl", "magma_tune() failed"); */
    /*     return status; */
    /* } */

    /* Set Mt, NT && NTRHS */
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
    MAGMA_ztrsmpl_Tile_Async(&descA, L, IPIV, &descB, sequence, &request);

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
 * MAGMA_ztrsmpl_Tile - Performs the forward substitution step of solving a system of linear equations
 * after the tile LU factorization of the matrix.
 * All matrices are passed through descriptors. All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] A
 *          The tile factor L from the factorization, computed by MAGMA_zgetrf.
 *
 * @param[in] L
 *          Auxiliary factorization data, related to the tile L factor, computed by MAGMA_zgetrf.
 *
 * @param[in] IPIV
 *          The pivot indices from MAGMA_zgetrf (not equivalent to LAPACK).
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
 * @sa MAGMA_ztrsmpl
 * @sa MAGMA_ztrsmpl_Tile_Async
 * @sa MAGMA_ctrsmpl_Tile
 * @sa MAGMA_dtrsmpl_Tile
 * @sa MAGMA_strsmpl_Tile
 * @sa MAGMA_zgetrf_Tile
 *
 ******************************************************************************/
int MAGMA_ztrsmpl_Tile(magma_desc_t *A, magma_desc_t *L, int *IPIV, magma_desc_t *B)
{
    magma_context_t *magma;
    magma_sequence_t *sequence = NULL;
    magma_request_t request = MAGMA_REQUEST_INITIALIZER;
    int status;

    magma = magma_context_self();
    if (magma == NULL) {
        magma_fatal_error("MAGMA_ztrsmpl_Tile", "MAGMA not initialized");
        return MAGMA_ERR_NOT_INITIALIZED;
    }
    magma_sequence_create(magma, &sequence);
    MAGMA_ztrsmpl_Tile_Async(A, L, IPIV, B, sequence, &request);
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
 *  MAGMA_ztrsmpl_Tile - Performs the forward substitution step of solving
 *  a system of linear equations after the tile LU factorization of the matrix.
 *  Non-blocking equivalent of MAGMA_ztrsmpl_Tile().
 *  Returns control to the user thread before worker threads finish the computation
 *  to allow for pipelined execution of diferent routines.
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
 * @sa MAGMA_ztrsmpl
 * @sa MAGMA_ztrsmpl_Tile
 * @sa MAGMA_ctrsmpl_Tile_Async
 * @sa MAGMA_dtrsmpl_Tile_Async
 * @sa MAGMA_strsmpl_Tile_Async
 * @sa MAGMA_zgetrf_Tile_Async
 *
 ******************************************************************************/
int MAGMA_ztrsmpl_Tile_Async(magma_desc_t *A, magma_desc_t *L, int *IPIV, magma_desc_t *B,
                              magma_sequence_t *sequence, magma_request_t *request)
{
    PLASMA_desc descA = A->desc;
    PLASMA_desc descB = B->desc;
    magma_context_t *magma;

    magma = magma_context_self();
    if (magma == NULL) {
        magma_fatal_error("MAGMA_ztrsmpl_Tile", "MAGMA not initialized");
        return MAGMA_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        magma_fatal_error("MAGMA_ztrsmpl_Tile", "NULL sequence");
        return MAGMA_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        magma_fatal_error("MAGMA_ztrsmpl_Tile", "NULL request");
        return MAGMA_ERR_UNALLOCATED;
    }
    /* Check sequence status */
    if (sequence->status == MAGMA_SUCCESS)
        request->status = MAGMA_SUCCESS;
    else
        return magma_request_fail(sequence, request, MAGMA_ERR_SEQUENCE_FLUSHED);

    /* Check descriptors for correctness */
    if (magma_desc_check( A ) != MAGMA_SUCCESS) {
        magma_error("MAGMA_ztrsmpl_Tile", "invalid first descriptor");
        return magma_request_fail(sequence, request, MAGMA_ERR_ILLEGAL_VALUE);
    }
    if (magma_desc_check( L ) != MAGMA_SUCCESS) {
        magma_error("MAGMA_ztrsmpl_Tile", "invalid second descriptor");
        return magma_request_fail(sequence, request, MAGMA_ERR_ILLEGAL_VALUE);
    }
    if (magma_desc_check( B ) != MAGMA_SUCCESS) {
        magma_error("MAGMA_ztrsmpl_Tile", "invalid third descriptor");
        return magma_request_fail(sequence, request, MAGMA_ERR_ILLEGAL_VALUE);
    }
    /* Check input arguments */
    if (descA.nb != descA.mb || descB.nb != descB.mb) {
        magma_error("MAGMA_ztrsmpl_Tile", "only square tiles supported");
        return magma_request_fail(sequence, request, MAGMA_ERR_ILLEGAL_VALUE);
    }
    /* Quick return */
/*
    if (min(N, NRHS) == 0)
        return PLASMA_SUCCESS;
*/
    magma_pztrsmpl( A, B, L, IPIV, sequence, request);

    return PLASMA_SUCCESS;
}
