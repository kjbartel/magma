/**
 *
 * @file descriptor.c
 *
 *  MAGMA auxiliary routines
 *  MAGMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 1.1.0
 * @author Mathieu Faverge
 * @date 2010-11-15
 *
 **/
#include <stdlib.h>
#include "common.h"

/***************************************************************************//**
 *  Internal static descriptor initializer
 **/
magma_desc_t magma_desc_init(MAGMA_enum dtyp, int mb, int nb, int bsiz,
                           int lm, int ln, int i, int j, int m, int n)
{
    magma_desc_t desc;
    PLASMA_desc *pdesc = &(desc.desc);

    // Matrix address
    pdesc->mat = NULL;
    pdesc->A21 = (lm - lm%mb)*(ln - ln%nb);
    pdesc->A12 = (     lm%mb)*(ln - ln%nb) + pdesc->A21;
    pdesc->A22 = (lm - lm%mb)*(     ln%nb) + pdesc->A12;
    // Matrix properties
    pdesc->dtyp = dtyp;
    pdesc->mb = mb;
    pdesc->nb = nb;
    pdesc->bsiz = bsiz;
    // Large matrix parameters
    pdesc->lm = lm;
    pdesc->ln = ln;
    // Large matrix derived parameters
    pdesc->lm1 = (lm/mb);
    pdesc->ln1 = (ln/nb);
    pdesc->lmt = (lm%mb==0) ? (lm/mb) : (lm/mb+1);
    pdesc->lnt = (ln%nb==0) ? (ln/nb) : (ln/nb+1);
    // Submatrix parameters
    pdesc->i = i;
    pdesc->j = j;
    pdesc->m = m;
    pdesc->n = n;
    // Submatrix derived parameters
    pdesc->mt = (i+m-1)/mb - i/mb + 1;
    pdesc->nt = (j+n-1)/nb - j/nb + 1;

    desc.occurences = 0;
    morse_desc_init( &desc );

    return desc;
}

/***************************************************************************//**
 *  Internal static descriptor initializer for submatrices
 **/
magma_desc_t magma_desc_submatrix(magma_desc_t descA, int i, int j, int m, int n)
{
    magma_desc_t descB;
    PLASMA_desc *pdescB = &(descB.desc);
    int mb, nb;

    descB = descA;
    mb = descA.desc.mb;
    nb = descA.desc.nb;
    // Submatrix parameters
    pdescB->i = i;
    pdescB->j = j;
    pdescB->m = m;
    pdescB->n = n;
    // Submatrix derived parameters
    pdescB->mt = (i+m-1)/mb - i/mb + 1;
    pdescB->nt = (j+n-1)/nb - j/nb + 1;

    morse_desc_submatrix( &descB );

    return descB;
}

/***************************************************************************//**
 *  Check for descriptor correctness
 **/
int magma_desc_check(magma_desc_t *desc)
{
    if (desc->desc.mat == NULL) {
        magma_error("magma_desc_check", "NULL matrix pointer");
        return MAGMA_ERR_UNALLOCATED;
    }
    if (desc->desc.dtyp != PlasmaRealFloat &&
        desc->desc.dtyp != PlasmaRealDouble &&
        desc->desc.dtyp != PlasmaComplexFloat &&
        desc->desc.dtyp != PlasmaComplexDouble  ) {
        magma_error("magma_desc_check", "invalid matrix type");
        return MAGMA_ERR_ILLEGAL_VALUE;
    }
    if (desc->desc.mb <= 0 || desc->desc.nb <= 0) {
        magma_error("magma_desc_check", "negative tile dimension");
        return MAGMA_ERR_ILLEGAL_VALUE;
    }
    if (desc->desc.bsiz < desc->desc.mb*desc->desc.nb) {
        magma_error("magma_desc_check", "tile memory size smaller than the product of dimensions");
        return MAGMA_ERR_ILLEGAL_VALUE;
    }
    if (desc->desc.lm <= 0 || desc->desc.ln <= 0) {
        magma_error("magma_desc_check", "negative matrix dimension");
        return MAGMA_ERR_ILLEGAL_VALUE;
    }
    if (desc->desc.i >= desc->desc.lm || desc->desc.j >= desc->desc.ln) {
        magma_error("magma_desc_check", "beginning of the matrix out of scope");
        return MAGMA_ERR_ILLEGAL_VALUE;
    }
    if (desc->desc.i+desc->desc.m > desc->desc.lm || desc->desc.j+desc->desc.n > desc->desc.ln) {
        magma_error("magma_desc_check", "submatrix out of scope");
        return MAGMA_ERR_ILLEGAL_VALUE;
    }
    return MAGMA_SUCCESS;
}

/***************************************************************************//**
 *
 **/
int magma_desc_mat_alloc( magma_desc_t *desc )
{
    size_t size;

    size = (size_t)desc->desc.lm * (size_t)desc->desc.ln * (size_t)plasma_element_size(desc->desc.dtyp);
    if ((desc->desc.mat = malloc(size)) == NULL) {
        magma_error("magma_desc_mat_alloc", "malloc() failed");
        return MAGMA_ERR_OUT_OF_RESOURCES;
    }

    morse_desc_create( desc );

    return MAGMA_SUCCESS;
}

/***************************************************************************//**
 *
 **/
int magma_desc_mat_free( magma_desc_t *desc)
{

    morse_desc_destroy( desc );

    if (desc->desc.mat != NULL) {
        free(desc->desc.mat);
        desc->desc.mat = NULL;
    }
    return MAGMA_SUCCESS;
}

/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  MAGMA_Desc_Create - Create matrix descriptor.
 *
 *******************************************************************************
 *
 * @param[out] desc
 *          On exit, descriptor of the matrix.
 *
 * @param[in] mat
 *          Memory location of the matrix.
 *
 * @param[in] dtyp
 *          Data type of the matrix:
 *          @arg MagmaRealFloat:     single precision real (S),
 *          @arg MagmaRealDouble:    double precision real (D),
 *          @arg MagmaComplexFloat:  single precision complex (C),
 *          @arg MagmaComplexDouble: double precision complex (Z).
 *
 * @param[in] mb
 *          Number of rows in a tile.
 *
 * @param[in] nb
 *          Number of columns in a tile.
 *
 * @param[in] bsiz
 *          Size in bytes including padding.
 *
 * @param[in] lm
 *          Number of rows of the entire matrix.
 *
 * @param[in] ln
 *          Number of columns of the entire matrix.
 *
 * @param[in] i
 *          Row index to the beginning of the submatrix.
 *
 * @param[in] j
 *          Column indes to the beginning of the submatrix.
 *
 * @param[in] m
 *          Number of rows of the submatrix.
 *
 * @param[in] n
 *          Number of columns of the submatrix.
 *
 *******************************************************************************
 *
 * @return
 *          \retval MAGMA_SUCCESS successful exit
 *
 ******************************************************************************/
int MAGMA_Desc_Create(magma_desc_t **desc, void *mat, MAGMA_enum dtyp, int mb, int nb, int bsiz,
                       int lm, int ln, int i, int j, int m, int n)
{
    magma_context_t *magma;
    int status;

    magma = magma_context_self();
    if (magma == NULL) {
        magma_error("MAGMA_Desc_Create", "MAGMA not initialized");
        return MAGMA_ERR_NOT_INITIALIZED;
    }
    /* Allocate memory and initialize the descriptor */
    *desc = (magma_desc_t*)malloc(sizeof(magma_desc_t));
    if (*desc == NULL) {
        magma_error("MAGMA_Desc_Create", "malloc() failed");
        return MAGMA_ERR_OUT_OF_RESOURCES;
    }
    **desc = magma_desc_init(dtyp, mb, nb, bsiz, lm, ln, i, j, m, n);
    (**desc).desc.mat = mat;

    /* Create scheduler structure like registering data */
    morse_desc_create( *desc );

    status = magma_desc_check(*desc);
    if (status != MAGMA_SUCCESS) {
        magma_error("MAGMA_Desc_Create", "invalid descriptor");
        return status;
    }
   
    return MAGMA_SUCCESS;
}

/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  MAGMA_Desc_Destroy - Destroys matrix descriptor.
 *
 *******************************************************************************
 *
 * @param[in] desc
 *          Matrix descriptor.
 *
 *******************************************************************************
 *
 * @return
 *          \retval MAGMA_SUCCESS successful exit
 *
 ******************************************************************************/
int MAGMA_Desc_Destroy(magma_desc_t **desc)
{
    magma_context_t *magma;

    magma = magma_context_self();
    if (magma == NULL) {
        magma_error("MAGMA_Desc_Destroy", "MAGMA not initialized");
        return MAGMA_ERR_NOT_INITIALIZED;
    }
    if (*desc == NULL) {
        magma_error("MAGMA_Desc_Destroy", "attempting to destroy a NULL descriptor");
        return MAGMA_ERR_UNALLOCATED;
    }

    /* Clean at scheduler level like unregistering data */
    morse_desc_destroy( *desc );

    free(*desc);
    *desc = NULL;
    return MAGMA_SUCCESS;
}



