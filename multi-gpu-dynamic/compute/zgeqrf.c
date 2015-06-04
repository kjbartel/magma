/**
 *
 * @file zgeqrf.c
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
 * @ingroup PLASMA_Complex64_t
 *
 *  MAGMA_zgeqrf - Computes the tile QR factorization of a complex M-by-N matrix A: A = Q * R.
 *
 *******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the matrix A. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix A.  N >= 0.
 *
 * @param[in,out] A
 *          On entry, the M-by-N matrix A.
 *          On exit, the elements on and above the diagonal of the array contain the min(M,N)-by-N
 *          upper trapezoidal matrix R (R is upper triangular if M >= N); the elements below the
 *          diagonal represent the unitary matrix Q as a product of elementary reflectors stored
 *          by tiles.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,M).
 *
 * @param[out] T
 *          On exit, auxiliary factorization data, required by MAGMA_zgeqrs to solve the system
 *          of equations.
 *
 *******************************************************************************
 *
 * @return
 *          \retval MAGMA_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *
 *******************************************************************************
 *
 * @sa MAGMA_zgeqrf_Tile
 * @sa MAGMA_zgeqrf_Tile_Async
 * @sa MAGMA_cgeqrf
 * @sa MAGMA_dgeqrf
 * @sa MAGMA_sgeqrf
 * @sa MAGMA_zgeqrs
 *
 ******************************************************************************/
int MAGMA_zgeqrf(int M, int N,
                  PLASMA_Complex64_t *A, int LDA,
                  magma_desc_t *T)
{
    int NB, IB, IBNB, MT, NT;
    int status;
    magma_context_t *magma;
    magma_sequence_t *sequence = NULL;
    magma_request_t request = MAGMA_REQUEST_INITIALIZER;
    magma_desc_t descA;

    magma = magma_context_self();
    if (magma == NULL) {
        magma_fatal_error("MAGMA_zgeqrf", "MAGMA not initialized");
        return MAGMA_ERR_NOT_INITIALIZED;
    }
    /* Check input arguments */
    if (M < 0) {
        magma_error("MAGMA_zgeqrf", "illegal value of M");
        return -1;
    }
    if (N < 0) {
        magma_error("MAGMA_zgeqrf", "illegal value of N");
        return -2;
    }
    if (LDA < max(1, M)) {
        magma_error("MAGMA_zgeqrf", "illegal value of LDA");
        return -4;
    }
    /* Quick return */
    if (min(M, N) == 0)
        return MAGMA_SUCCESS;

    /* Tune NB && IB depending on M, N && NRHS; Set NBNBSIZE */
    /* status = magma_tune(MAGMA_FUNC_ZGELS, M, N, 0); */
    /* if (status != MAGMA_SUCCESS) { */
    /*     magma_error("MAGMA_zgeqrf", "magma_tune() failed"); */
    /*     return status; */
    /* } */

    /* Set MT && NT */
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
    MAGMA_zgeqrf_Tile_Async(&descA, T, sequence, &request);

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
 * @ingroup PLASMA_Complex64_t_Tile
 *
 *  MAGMA_zgeqrf_Tile - Computes the tile QR factorization of a matrix.
 *  Tile equivalent of MAGMA_zgeqrf().
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in,out] A
 *          On entry, the M-by-N matrix A.
 *          On exit, the elements on and above the diagonal of the array contain the min(M,N)-by-N
 *          upper trapezoidal matrix R (R is upper triangular if M >= N); the elements below the
 *          diagonal represent the unitary matrix Q as a product of elementary reflectors stored
 *          by tiles.
 *
 * @param[out] T
 *          On exit, auxiliary factorization data, required by MAGMA_zgeqrs to solve the system
 *          of equations.
 *
 *******************************************************************************
 *
 * @return
 *          \retval MAGMA_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa MAGMA_zgeqrf
 * @sa MAGMA_zgeqrf_Tile_Async
 * @sa MAGMA_cgeqrf_Tile
 * @sa MAGMA_dgeqrf_Tile
 * @sa MAGMA_sgeqrf_Tile
 * @sa MAGMA_zgeqrs_Tile
 *
 ******************************************************************************/
int MAGMA_zgeqrf_Tile(magma_desc_t *A, magma_desc_t *T)
{
    magma_context_t *magma;
    magma_sequence_t *sequence = NULL;
    magma_request_t request = MAGMA_REQUEST_INITIALIZER;
    int status;

    magma = magma_context_self();
    if (magma == NULL) {
        magma_fatal_error("MAGMA_zgeqrf_Tile", "MAGMA not initialized");
        return MAGMA_ERR_NOT_INITIALIZED;
    }
    magma_sequence_create(magma, &sequence);
    MAGMA_zgeqrf_Tile_Async(A, T, sequence, &request);
    morse_barrier( magma );
    morse_desc_getoncpu( A );
    morse_desc_getoncpu( T );
    status = sequence->status;
    magma_sequence_destroy(magma, sequence);
    return status;
}

/***************************************************************************//**
 *
 * @ingroup PLASMA_Complex64_t_Tile_Async
 *
 *  MAGMA_zgeqrf_Tile_Async - Computes the tile QR factorization of a matrix.
 *  Non-blocking equivalent of MAGMA_zgeqrf_Tile().
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
 * @sa MAGMA_zgeqrf
 * @sa MAGMA_zgeqrf_Tile
 * @sa MAGMA_cgeqrf_Tile_Async
 * @sa MAGMA_dgeqrf_Tile_Async
 * @sa MAGMA_sgeqrf_Tile_Async
 * @sa MAGMA_zgeqrs_Tile_Async
 *
 ******************************************************************************/
int MAGMA_zgeqrf_Tile_Async(magma_desc_t *A, magma_desc_t *T,
                             magma_sequence_t *sequence, magma_request_t *request)
{
    PLASMA_desc descA = A->desc;
    magma_context_t *magma;

    magma = magma_context_self();
    if (magma == NULL) {
        magma_error("MAGMA_zgeqrf_Tile", "MAGMA not initialized");
        return MAGMA_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        magma_fatal_error("MAGMA_zgeqrf_Tile", "NULL sequence");
        return MAGMA_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        magma_fatal_error("MAGMA_zgeqrf_Tile", "NULL request");
        return MAGMA_ERR_UNALLOCATED;
    }
    /* Check sequence status */
    if (sequence->status == MAGMA_SUCCESS)
        request->status = MAGMA_SUCCESS;
    else
        return magma_request_fail(sequence, request, MAGMA_ERR_SEQUENCE_FLUSHED);

    /* Check descriptors for correctness */
    if (magma_desc_check( A ) != MAGMA_SUCCESS) {
        magma_error("MAGMA_zgeqrf_Tile", "invalid first descriptor");
        return magma_request_fail(sequence, request, MAGMA_ERR_ILLEGAL_VALUE);
    }
    if (magma_desc_check( T ) != MAGMA_SUCCESS) {
        magma_error("MAGMA_zgeqrf_Tile", "invalid second descriptor");
        return magma_request_fail(sequence, request, MAGMA_ERR_ILLEGAL_VALUE);
    }
    /* Check input arguments */
    if (descA.nb != descA.mb) {
        magma_error("MAGMA_zgeqrf_Tile", "only square tiles supported");
        return magma_request_fail(sequence, request, MAGMA_ERR_ILLEGAL_VALUE);
    }
    /* Quick return */
/*
    if (min(M, N) == 0)
        return MAGMA_SUCCESS;
*/
    /* if (magma->householder == MAGMA_FLAT_HOUSEHOLDER) { */

    magma_pzgeqrf( A, T, sequence, request);

    /* } */
    /* else { */
    /*     magma_dynamic_call_5(magma_pzgeqrfrh, */
    /*         magma_desc_t, descA, */
    /*         magma_desc_t, descT, */
    /*         MAGMA_enum, MAGMA_RHBLK, */
    /*         magma_sequence_t*, sequence, */
    /*         magma_request_t*, request); */
    /* } */

    return MAGMA_SUCCESS;
}
