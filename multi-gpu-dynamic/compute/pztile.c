/**
 *
 * @file ztile.c
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.3.1
 * @author Jakub Kurzak
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 **/
#include "common.h"

#define AF77(m, n) &(Af77[ ((int64_t)A.nb*(int64_t)lda*(int64_t)(n)) + (int64_t)(A.mb*(m)) ])

#undef BLKADDR
#define BLKADDR(A, type, m, n)  (type *)magma_getaddr(A, m, n)
#define ABDL(m, n) BLKADDR(dA, PLASMA_Complex64_t, m, n)

void CORE_ztile_zero(int X1, int X2, int Y1, int Y2, 
                     PLASMA_Complex64_t *A, int lda);

/***************************************************************************//**
 *  Conversion from LAPACK F77 matrix layout to tile layout - dynamic scheduling
 **/
void magma_pzlapack_to_tile(PLASMA_Complex64_t *Af77, int lda, magma_desc_t *dA,
                            magma_sequence_t *sequence, magma_request_t *request)
{
    PLASMA_Complex64_t *f77;
    PLASMA_Complex64_t *bdl;
    PLASMA_desc A = dA->desc;
    magma_context_t *magma;
    int X1, Y1;
    int X2, Y2;
    int n, m, ldt;

    magma = magma_context_self();
    if (sequence->status != MAGMA_SUCCESS)
        return;

    for (m = 0; m < A.mt; m++)
    {
        ldt = BLKLDD(dA, m);
        for (n = 0; n < A.nt; n++)
        {
            X1 = n == 0 ? A.j%A.nb : 0;
            Y1 = m == 0 ? A.i%A.mb : 0;
            X2 = n == A.nt-1 ? (A.j+A.n-1)%A.nb+1 : A.nb;
            Y2 = m == A.mt-1 ? (A.i+A.m-1)%A.mb+1 : A.mb;

            f77 = AF77(m, n);
            bdl = ABDL(m, n);
            CORE_zlacpy(
                PlasmaUpperLower, (Y2-Y1), (X2-X1),
                &(f77[X1*lda+Y1]), lda, 
                &(bdl[X1*lda+Y1]), ldt);
        }
    }
}

/***************************************************************************//**
 *  Conversion from LAPACK F77 matrix layout to tile layout - dynamic scheduling
 **/
void magma_pztile_to_lapack(magma_desc_t *dA, PLASMA_Complex64_t *Af77, int lda,
                            magma_sequence_t *sequence, magma_request_t *request)
{
    PLASMA_Complex64_t *f77;
    PLASMA_Complex64_t *bdl;
    PLASMA_desc A = dA->desc;
    magma_context_t *magma;
    int X1, Y1;
    int X2, Y2;
    int n, m, ldt;

    magma = magma_context_self();
    if (sequence->status != MAGMA_SUCCESS)
        return;

    for (m = 0; m < A.mt; m++)
    {
        ldt = BLKLDD(dA, m);
        for (n = 0; n < A.nt; n++)
        {
            X1 = n == 0 ? A.j%A.nb : 0;
            Y1 = m == 0 ? A.i%A.mb : 0;
            X2 = n == A.nt-1 ? (A.j+A.n-1)%A.nb+1 : A.nb;
            Y2 = m == A.mt-1 ? (A.i+A.m-1)%A.mb+1 : A.mb;

            f77 = AF77(m, n);
            bdl = ABDL(m, n);
            CORE_zlacpy(
                PlasmaUpperLower, (Y2-Y1), (X2-X1),
                &(bdl[X1*lda+Y1]), ldt,
                &(f77[X1*lda+Y1]), lda);
        }
    }
}

/***************************************************************************//**
 *  Zeroes a submatrix in tile layout - dynamic scheduling
 **/
void magma_pztile_zero(magma_desc_t *dA, magma_sequence_t *sequence, magma_request_t *request)
{
    PLASMA_Complex64_t *bdl;
    PLASMA_desc A = dA->desc;
    magma_context_t *magma;
    int X1, Y1;
    int X2, Y2;
    int n, m, ldt;

    magma = magma_context_self();
    if (sequence->status != MAGMA_SUCCESS)
        return;

    for (m = 0; m < A.mt; m++)
    {
        ldt = BLKLDD(dA, m);
        for (n = 0; n < A.nt; n++)
        {
            X1 = n == 0 ? A.j%A.nb : 0;
            Y1 = m == 0 ? A.i%A.mb : 0;
            X2 = n == A.nt-1 ? (A.j+A.n-1)%A.nb+1 : A.nb;
            Y2 = m == A.mt-1 ? (A.i+A.m-1)%A.mb+1 : A.mb;

            bdl = ABDL(m, n);
 
            CORE_ztile_zero(X1, X2, Y1, Y2, bdl, ldt);
        }
    }
}

/***************************************************************************//**
 *
 **/
void CORE_ztile_zero(int X1, int X2, int Y1, int Y2, 
                     PLASMA_Complex64_t *A, int lda)
{
    int x, y;

    for (x = X1; x < X2; x++)
        for (y = Y1; y < Y2; y++)
            A[lda*x+y] = 0.0;
}
