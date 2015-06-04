/**
 * \file zgesv.c
 *
 *  PLASMA computational routines
 *  Release Date: November, 15th 2009
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

void plasma_memzero(void *memptr, PLASMA_size size, int type);

/***************************************************************************//**
 *
 * @ingroup MAGMA_Complex64_t
 *
 *  MAGMA_zgesv - Computes the solution to a system of linear equations A * X = B,
 *  where A is an N-by-N matrix and X and B are N-by-NRHS matrices.
 *  The tile LU decomposition with partial tile pivoting and row interchanges is used to factor A.
 *  The factored form of A is then used to solve the system of equations A * X = B.
 *
 *******************************************************************************
 *
 * @param[in] N
 *          The number of linear equations, i.e., the order of the matrix A. N >= 0.
 *
 * @param[in] NRHS
 *          The number of right hand sides, i.e., the number of columns of the matrix B.
 *          NRHS >= 0.
 *
 * @param[in,out] A
 *          On entry, the N-by-N coefficient matrix A.
 *          On exit, the tile L and U factors from the factorization (not equivalent to LAPACK).
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,N).
 *
 * @param[out] L
 *          On exit, auxiliary factorization data, related to the tile L factor,
 *          necessary to solve the system of equations.
 *
 * @param[out] IPIV
 *          On exit, the pivot indices that define the permutations (not equivalent to LAPACK).
 *
 * @param[in,out] B
 *          On entry, the N-by-NRHS matrix of right hand side matrix B.
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
 *          \retval >0 if i, U(i,i) is exactly zero. The factorization has been completed,
 *               but the factor U is exactly singular, so the solution could not be computed.
 *
 *******************************************************************************
 *
 * @sa MAGMA_zgesv_Tile
 * @sa MAGMA_zgesv_Tile_Async
 * @sa MAGMA_cgesv
 * @sa MAGMA_dgesv
 * @sa MAGMA_sgesv
 * @sa MAGMA_zcgesv
 *
 ******************************************************************************/
int MAGMA_zgesv(int N, int NRHS,
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
        magma_error("MAGMA_zgesv", "MAGMA not initialized");
        return MAGMA_ERR_NOT_INITIALIZED;
    }
    /* Check input arguments */
    if (N < 0) {
        magma_error("MAGMA_zgesv", "illegal value of N");
        return -1;
    }
    if (NRHS < 0) {
        magma_error("MAGMA_zgesv", "illegal value of NRHS");
        return -2;
    }
    if (LDA < max(1, N)) {
        magma_error("MAGMA_zgesv", "illegal value of LDA");
        return -4;
    }
    if (LDB < max(1, N)) {
        magma_error("MAGMA_zgesv", "illegal value of LDB");
        return -8;
    }
    /* Quick return */
    if (min(N, NRHS) == 0)
        return MAGMA_SUCCESS;

    /* Tune NB && IB depending on M, N && NRHS; Set NBNB */
    /* status = plasma_tune(PLASMA_FUNC_ZGESV, N, N, NRHS); */
    /* if (status != PLASMA_SUCCESS) { */
    /*     plasma_error("PLASMA_zgesv", "plasma_tune() failed"); */
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
    MAGMA_zgesv_Tile_Async(&descA, L, IPIV, &descB, sequence, &request);

    /* if ( MAGMA_TRANSLATION == MAGMA_OUTOFPLACE ) { */
        magma_zooptile2lap( descA, A, NB, NB, LDA, N    );
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
 *  MAGMA_zgesv_Tile - Solves a system of linear equations using the tile LU factorization.
 *  Tile equivalent of MAGMA_zgetrf().
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in,out] A
 *          On entry, the N-by-N coefficient matrix A.
 *          On exit, the tile L and U factors from the factorization (not equivalent to LAPACK).
 *
 * @param[in,out] L
 *          On exit, auxiliary factorization data, related to the tile L factor,
 *          necessary to solve the system of equations.
 *
 * @param[out] IPIV
 *          On exit, the pivot indices that define the permutations (not equivalent to LAPACK).
 *
 * @param[in,out] B
 *          On entry, the N-by-NRHS matrix of right hand side matrix B.
 *          On exit, if return value = 0, the N-by-NRHS solution matrix X.
 *
 *******************************************************************************
 *
 * @return
 *          \retval MAGMA_SUCCESS successful exit
 *          \retval >0 if i, U(i,i) is exactly zero. The factorization has been completed,
 *               but the factor U is exactly singular, so the solution could not be computed.
 *
 *******************************************************************************
 *
 * @sa MAGMA_zgesv
 * @sa MAGMA_zgesv_Tile_Async
 * @sa MAGMA_cgesv_Tile
 * @sa MAGMA_dgesv_Tile
 * @sa MAGMA_sgesv_Tile
 * @sa MAGMA_zcgesv_Tile
 *
 ******************************************************************************/
int MAGMA_zgesv_Tile(magma_desc_t *A, magma_desc_t *L, int *IPIV, magma_desc_t *B)
{
    magma_context_t *magma;
    magma_sequence_t *sequence = NULL;
    magma_request_t request = MAGMA_REQUEST_INITIALIZER;
    int status;

    magma = magma_context_self();
    if (magma == NULL) {
        magma_fatal_error("MAGMA_zgesv_Tile", "MAGMA not initialized");
        return MAGMA_ERR_NOT_INITIALIZED;
    }
    magma_sequence_create(magma, &sequence);
    MAGMA_zgesv_Tile_Async(A, L, IPIV, B, sequence, &request);
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
 *  MAGMA_zgesv_Tile_Async - Solves a system of linear equations using the tile
 *  LU factorization.
 *  Non-blocking equivalent of MAGMA_zgesv_Tile().
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
 * @sa MAGMA_zgesv
 * @sa MAGMA_zgesv_Tile
 * @sa MAGMA_cgesv_Tile_Async
 * @sa MAGMA_dgesv_Tile_Async
 * @sa MAGMA_sgesv_Tile_Async
 * @sa MAGMA_zcgesv_Tile_Async
 *
 ******************************************************************************/
int MAGMA_zgesv_Tile_Async(magma_desc_t *A, magma_desc_t *L, int *IPIV, magma_desc_t *B,
                            magma_sequence_t *sequence, magma_request_t *request)
{
    PLASMA_desc descA = A->desc;
    PLASMA_desc descL = L->desc;
    PLASMA_desc descB = B->desc;
    magma_context_t *magma;

    magma = magma_context_self();
    if (magma == NULL) {
        magma_fatal_error("MAGMA_zgesv_Tile", "MAGMA not initialized");
        return MAGMA_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        magma_fatal_error("MAGMA_zgesv_Tile", "NULL sequence");
        return MAGMA_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        magma_fatal_error("MAGMA_zgesv_Tile", "NULL request");
        return MAGMA_ERR_UNALLOCATED;
    }
    /* Check sequence status */
    if (sequence->status == MAGMA_SUCCESS)
        request->status = MAGMA_SUCCESS;
    else
        return magma_request_fail(sequence, request, MAGMA_ERR_SEQUENCE_FLUSHED);

    /* Check descriptors for correctness */
    if (magma_desc_check( A ) != MAGMA_SUCCESS) {
        magma_error("MAGMA_zgesv_Tile", "invalid first descriptor");
        return magma_request_fail(sequence, request, MAGMA_ERR_ILLEGAL_VALUE);
    }
    if (magma_desc_check( L ) != MAGMA_SUCCESS) {
        magma_error("MAGMA_zgesv_Tile", "invalid second descriptor");
        return magma_request_fail(sequence, request, MAGMA_ERR_ILLEGAL_VALUE);
    }
    if (magma_desc_check( B ) != MAGMA_SUCCESS) {
        magma_error("MAGMA_zgesv_Tile", "invalid third descriptor");
        return magma_request_fail(sequence, request, MAGMA_ERR_ILLEGAL_VALUE);
    }
    /* Check input arguments */
    if (descA.nb != descA.mb || descB.nb != descB.mb) {
        magma_error("MAGMA_zgesv_Tile", "only square tiles supported");
        return magma_request_fail(sequence, request, MAGMA_ERR_ILLEGAL_VALUE);
    }
    /* Quick return */
/*
    if (min(N, NRHS) == 0)
        return PLASMA_SUCCESS;
*/
    /* Clear IPIV and Lbdl */
    plasma_memzero(IPIV,      descA.mt*descA.nt*descA.nb,          PlasmaInteger);
    plasma_memzero(descL.mat, descL.mt*descL.nt*descL.mb*descL.nb, PlasmaComplexDouble);

    magma_pzgetrf( A, L, IPIV, sequence, request);
    magma_pztrsmpl(A, B, L, IPIV, sequence, request);
    magma_pztrsm( PlasmaLeft, PlasmaUpper, PlasmaNoTrans, PlasmaNonUnit,
                  1.0, A, B,
                  sequence, request);

    return PLASMA_SUCCESS;
}
