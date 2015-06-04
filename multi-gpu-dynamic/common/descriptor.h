/**
 *
 * @file descriptor.h
 *
 *  MAGMA auxiliary routines
 *  MAGMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 1.1.0
 * @author Jakub Kurzak
 * @author Mathieu Faverge
 * @date 2010-11-15
 *
 **/
#ifndef _MAGMA_DESCRIPTOR_H_
#define _MAGMA_DESCRIPTOR_H_

#ifdef __cplusplus
extern "C" {
#endif

/***************************************************************************//**
 *  Internal routines
 **/
inline static void *magma_geteltaddr( const magma_desc_t *A, int m, int n, int eltsize);
inline static void *magma_getaddr(magma_desc_t *A, int m, int n);
magma_desc_t magma_desc_init(MAGMA_enum dtyp, int mb, int nb, int bsiz, int lm, int ln, int i, int j, int m, int n);
magma_desc_t magma_desc_submatrix(magma_desc_t descA, int i, int j, int m, int n);
int magma_desc_check(magma_desc_t *desc);
int magma_desc_mat_alloc(magma_desc_t *desc);
int magma_desc_mat_free(magma_desc_t *desc);

extern int plasma_element_size(int type);

#define BLKLDD(A, k) ( ( (k) + (A)->desc.i/(A)->desc.mb) < (A)->desc.lm1 ? (A)->desc.mb : (A)->desc.lm%(A)->desc.mb )

/***************************************************************************//**
 *  Internal function to return adress of block (m,n)
 **/
inline static void *magma_getaddr(magma_desc_t *descA, int m, int n)
{
    PLASMA_desc *A = &(descA->desc);
    size_t mm = m+A->i/A->mb;
    size_t nn = n+A->j/A->nb;
    size_t eltsize = plasma_element_size(A->dtyp);
    size_t offset = 0;

    if (mm < A->lm1) {
        if (nn < A->ln1)
            offset = A->bsiz*(mm+A->lm1*nn);
        else
            offset = A->A12 + (A->mb*(A->ln%A->nb)*mm);
    }
    else {
        if (nn < A->ln1)
            offset = A->A21 + ((A->lm%A->mb)*A->nb*nn);
        else
            offset = A->A22;
    }

    return (void*)((intptr_t)A->mat + (offset*eltsize) );
}

/***************************************************************************//**
 *  Internal function to return adress of element A(m,n)
 **/
inline static void *magma_geteltaddr( const magma_desc_t *descA, int m, int n, int eltsize)
{
    const PLASMA_desc *A = &(descA->desc);
    size_t mm = m/A->mb;
    size_t nn = n/A->nb;
    size_t offset = 0;

    if (mm < A->lm1) {
        if (nn < A->ln1)
            offset = A->bsiz*(mm+A->lm1*nn) + m%A->mb + A->mb*(n%A->nb);
        else
            offset = A->A12 + (A->mb*(A->ln%A->nb)*mm) + m%A->mb + A->mb*(n%A->nb);
    }
    else {
        if (nn < A->ln1)
            offset = A->A21 + ((A->lm%A->mb)*A->nb*nn) + m%A->mb + (A->lm%A->mb)*(n%A->nb);
        else
            offset = A->A22 + m%A->mb  + (A->lm%A->mb)*(n%A->nb);
    }
    return (void*)((intptr_t)A->mat + (offset*eltsize) );
}

#ifdef __cplusplus
}
#endif

#endif
