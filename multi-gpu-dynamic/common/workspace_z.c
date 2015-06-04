/**
 *
 * @file workspace_z.c
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.4.5
 * @author Jakub Kurzak
 * @author Hatem Ltaief
 * @author Azzam Haidar
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 **/
#include "common.h"
#include "compute_z.h"

/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  MAGMA_Alloc_Workspace_zgeqrf - Allocates workspace for MAGMA_zgeqrf.
 *
 *******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the matrix A. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix A.  N >= 0.
 *
 * @param[out] T
 *          On exit, workspace handle for storage of the extra T factors required by the tile QR
 *          factorization.
 *
 *******************************************************************************
 *
 * @return
 *          \retval MAGMA_SUCCESS successful exit
 *
 ******************************************************************************/
int MAGMA_Alloc_Workspace_zgeqrf(int M, int N, magma_desc_t **T) 
{
    magma_context_t *magma;
    int NB, IB, MT, NT;
    int64_t lm, ln;

    magma = magma_context_self();
    if (magma == NULL) {
        magma_error("MAGMA_Alloc_Workspace_zgeqrf", "MAGMA not initialized");
        return MAGMA_ERR_NOT_INITIALIZED;
    }

    NB = MAGMA_NB;
    IB = MAGMA_IB;

    MT = (M%NB==0) ? (M/NB) : (M/NB+1);
    NT = (N%NB==0) ? (N/NB) : (N/NB+1);

    *T = (magma_desc_t*)malloc(sizeof(magma_desc_t));

    if ( magma->householder == MAGMA_TREE_HOUSEHOLDER )
        NT *= 2;

    lm = IB * MT;
    ln = NB * NT;

    **T = magma_desc_init(
                          PlasmaComplexDouble, IB, NB, IB*NB, 
                          lm, ln, 0, 0, lm, ln );
    if ( magma_desc_mat_alloc( *T ) ) {
        magma_error( __func__, "magma_shared_alloc() failed");
        free(*T);
        return MAGMA_ERR_OUT_OF_RESOURCES;
    }

    return MAGMA_SUCCESS;
}


/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  MAGMA_Alloc_Workspace_zgeev - Allocates workspace for MAGMA_zgeev.
 *
 *******************************************************************************
 *
 * @param[in] N
 *          The order of the matrix A.  N >= 0.
 *
 * @param[out] T
 *          On exit, workspace handle for storage of the extra T factors.
 *
 *******************************************************************************
 *
 * @return
 *          \retval MAGMA_SUCCESS successful exit
 *
 ******************************************************************************/
int MAGMA_Alloc_Workspace_zgeev(int N, magma_desc_t **T) {
    return MAGMA_Alloc_Workspace_zgeqrf(N, N, T); 
}

/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  MAGMA_Alloc_Workspace_zgehrd - Allocates workspace for MAGMA_zgehrd.
 *
 *******************************************************************************
 *
 * @param[in] N
 *          The number of columns of the matrix A.  N >= 0.
 *
 * @param[out] T
 *          On exit, workspace handle for storage of the extra T factors.
 *
 *******************************************************************************
 *
 * @return
 *          \retval MAGMA_SUCCESS successful exit
 *
 ******************************************************************************/
int MAGMA_Alloc_Workspace_zgehrd(int N, magma_desc_t **T) {
    return MAGMA_Alloc_Workspace_zgeqrf(N, N, T); 
}

/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  MAGMA_Alloc_Workspace_zgebrd - Allocates workspace for MAGMA_zgebrd.
 *
 *******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the matrix A. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix A.  N >= 0.
 *
 * @param[out] T
 *          On exit, workspace handle for storage of the extra T factors.
 *
 *******************************************************************************
 *
 * @return
 *          \retval MAGMA_SUCCESS successful exit
 *
 ******************************************************************************/
int MAGMA_Alloc_Workspace_zgebrd(int M, int N, magma_desc_t **T) {
    return MAGMA_Alloc_Workspace_zgeqrf(M, N, T); 
}

/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  MAGMA_Alloc_Workspace_zgels - Allocates workspace for MAGMA_zgels.
 *
 *******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the matrix A. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix A.  N >= 0.
 *
 * @param[out] T
 *          On exit, workspace handle for storage of the extra T factors required by the tile QR
 *          or the tile LQ factorization.
 *
 *******************************************************************************
 *
 * @return
 *          \retval MAGMA_SUCCESS successful exit
 *
 ******************************************************************************/
int MAGMA_Alloc_Workspace_zgels(int M, int N, magma_desc_t **T) {
    return MAGMA_Alloc_Workspace_zgeqrf(M, N, T); 
}

/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  MAGMA_Alloc_Workspace_zgelqf - Allocates workspace for MAGMA_zgelqf.
 *
 *******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the matrix A. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix A.  N >= 0.
 *
 * @param[out] T
 *          On exit, workspace handle for storage of the extra T factors required by the tile LQ
 *          factorization.
 *
 *******************************************************************************
 *
 * @return
 *          \retval MAGMA_SUCCESS successful exit
 *
 ******************************************************************************/
int MAGMA_Alloc_Workspace_zgelqf(int M, int N, magma_desc_t **T) {
    return MAGMA_Alloc_Workspace_zgeqrf(M, N, T); 
}

/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  MAGMA_Alloc_Workspace_zgesvd - Allocates workspace for MAGMA_zgesvd.
 *
 *******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the matrix A. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix A.  N >= 0.
 *
 * @param[out] T
 *          On exit, workspace handle for storage of the extra T factors required by the tile BRD.
 *
 *******************************************************************************
 *
 * @return
 *          \retval MAGMA_SUCCESS successful exit
 *
 ******************************************************************************/
int MAGMA_Alloc_Workspace_zgesvd(int M, int N, magma_desc_t **T) {
    return MAGMA_Alloc_Workspace_zgeqrf(M, N, T); 
}

/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  MAGMA_Alloc_Workspace_zheev - Allocates workspace for MAGMA_zheev or MAGMA_zheev_Tile routine.
 *
 *******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the matrix A. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix A.  N >= 0.
 *
 * @param[out] T
 *          On exit, workspace handle for storage of the extra T factors required by the tile TRD.
 *
 *******************************************************************************
 *
 * @return
 *          \retval MAGMA_SUCCESS successful exit
 *
 ******************************************************************************/
int MAGMA_Alloc_Workspace_zheev(int M, int N, magma_desc_t **T) {
    return MAGMA_Alloc_Workspace_zgeqrf(M, N, T); 
}

/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  MAGMA_Alloc_Workspace_zhegv - Allocates workspace for MAGMA_zhegv.
 *
 *******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the matrix A. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix A.  N >= 0.
 *
 * @param[out] T
 *          On exit, workspace handle for storage of the extra T factors.
 *
 *******************************************************************************
 *
 * @return
 *          \retval MAGMA_SUCCESS successful exit
 *
 ******************************************************************************/
int MAGMA_Alloc_Workspace_zhegv(int M, int N, magma_desc_t **T) {
    return MAGMA_Alloc_Workspace_zgeqrf(M, N, T); 
}

 /***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  MAGMA_Alloc_Workspace_zhetrd - Allocates workspace for MAGMA_zhetrd or MAGMA_zhetrd_Tile routine.
 *
 *******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the matrix A. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix A.  N >= 0.
 *
 * @param[out] T
 *          On exit, workspace handle for storage of the extra T factors required by the tile TRD.
 *
 *******************************************************************************
 *
 * @return
 *          \retval MAGMA_SUCCESS successful exit
 *
 ******************************************************************************/
int MAGMA_Alloc_Workspace_zhetrd(int M, int N, magma_desc_t **T) {
    return MAGMA_Alloc_Workspace_zgeqrf(M, N, T); 
}


/***************************************************************************//**
 * @ingroup Auxiliary
 *
 *  MAGMA_Alloc_Workspace_zgetrf_incpiv - Allocates workspace for
 *  MAGMA_zgetrf_incpiv or MAGMA_zgetrf_incpiv_Tile or
 *  MAGMA_zgetrf_incpiv_Tile_Async routines.
 *
 *******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the matrix A. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix A.  N >= 0.
 *
 * @param[out] L
 *          On exit, workspace handle for storage of the extra L factors required by the tile LU
 *          factorization.
 *
 * @param[out] IPIV
 *          On exit, workspace handle for storage of pivot indexes required by the tile LU
 *          factorization (not equivalent to LAPACK).
 *
 *******************************************************************************
 *
 * @return
 *          \retval MAGMA_SUCCESS successful exit
 *
 ******************************************************************************
 *
 * @sa MAGMA_zgetrf_incpiv
 * @sa MAGMA_zgetrf_incpiv_Tile
 * @sa MAGMA_zgetrf_incpiv_Tile_Async
 *
 ******************************************************************************/
int MAGMA_Alloc_Workspace_zgetrf_incpiv(int M, int N, magma_desc_t **L, int **IPIV) {
    magma_context_t *magma;
    int NB, IB, MT, NT;
    int64_t lm, ln;

    magma = magma_context_self();
    if (magma == NULL) {
        magma_error("MAGMA_Alloc_Workspace_zgetrf_incpiv", "MAGMA not initialized");
        return MAGMA_ERR_NOT_INITIALIZED;
    }

    NB = MAGMA_NB;
    IB = MAGMA_IB * 2; /* TODO: check if *2 is not only required for GPU */

    MT = (M%NB==0) ? (M/NB) : (M/NB+1);
    NT = (N%NB==0) ? (N/NB) : (N/NB+1);

    lm = IB * MT;
    ln = NB * NT;

    *L    = (magma_desc_t*)malloc(sizeof(magma_desc_t));
    *IPIV = (int*)malloc( min(MT, NT) * NB * NT * sizeof(int) );

    **L = magma_desc_init(
                          PlasmaComplexDouble, IB, NB, IB*NB, 
                          lm, ln, 0, 0, lm, ln );
    if ( magma_desc_mat_alloc( *L ) ) {
        magma_error( __func__, "magma_shared_alloc() failed");
        free(*L);
        return MAGMA_ERR_OUT_OF_RESOURCES;
    }

    return MAGMA_SUCCESS;
}

/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  MAGMA_Alloc_Workspace_zgesv - Allocates workspace for MAGMA_zgesv or MAGMA_zgesv_Tile routines.
 *
 *******************************************************************************
 *
 * @param[in] N
 *          The number of linear equations, i.e., the order of the matrix A. N >= 0.
 *
 * @param[out] L
 *          On exit, workspace handle for storage of the extra L factors required by the tile LU
 *          factorization.
 *
 * @param[out] IPIV
 *          On exit, workspace handle for storage of pivot indexes required by the tile LU
 *          factorization (not equivalent to LAPACK).
 *
 *******************************************************************************
 *
 * @return
 *          \retval MAGMA_SUCCESS successful exit
 *
 ******************************************************************************/
int MAGMA_Alloc_Workspace_zgesv_incpiv(int N, magma_desc_t **L, int **IPIV) {
    return MAGMA_Alloc_Workspace_zgetrf_incpiv(N, N, L, IPIV);
}
