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

#define A(m, n) dA, m, n
#define B(m, n) dB, m, n

void CORE_ztile_zero(int X1, int X2, int Y1, int Y2, 
                     PLASMA_Complex64_t *A, int lda);

/***************************************************************************//**
 *  Conversion from LAPACK F77 matrix layout to tile layout - dynamic scheduling
 **/
void magma_pzlapack_to_tile(PLASMA_Complex64_t *Af77, int lda, magma_desc_t *dA,
                            magma_sequence_t *sequence, magma_request_t *request)
{
    magma_context_t *magma;
    MorseOption_t options;
    magma_desc_t *dB;
    int X1, Y1;
    int X2, Y2;
    int n, m;

    magma = magma_context_self();
    if (sequence->status != MAGMA_SUCCESS)
        return;

    morse_options_init( &options, magma, sequence, request );

    *dB = magma_desc_init(
        PlasmaComplexDouble, dA->mb, dA->nb, dA->bsiz,
        dA->m, dA->n, dA->i, dA->j, dA->m, dA->n);

    dB->get_blkaddr = magma_getaddr_cm;
    dB->get_blkldd  = magma_get_blkldd_cm;
    dB->styp = PlasmaCM;

    for (m = 0; m < dA->mt; m++)
    {
        for (n = 0; n < dA->nt; n++)
        {
            X1 = n == 0 ? dA->j%dA->nb : 0;
            Y1 = m == 0 ? dA->i%dA->mb : 0;
            X2 = n == dA->nt-1 ? (dA->j+dA->n-1)%dA->nb+1 : dA->nb;
            Y2 = m == dA->mt-1 ? (dA->i+dA->m-1)%dA->mb+1 : dA->mb;

            MORSE_zlacpy(
                &options,
                PlasmaUpperLower, (Y2-Y1), (X2-X1),
                A(m, n), B(m, n));
        }
    }
}

/***************************************************************************//**
 *  Conversion from LAPACK F77 matrix layout to tile layout - dynamic scheduling
 **/
void magma_pztile_to_lapack(magma_desc_t *dA, PLASMA_Complex64_t *Af77, int lda,
                            magma_sequence_t *sequence, magma_request_t *request)
{
    magma_context_t *magma;
    MorseOption_t options;
    magma_desc_t *dB;
    int X1, Y1;
    int X2, Y2;
    int n, m;

    magma = magma_context_self();
    if (sequence->status != MAGMA_SUCCESS)
        return;

    morse_options_init( &options, magma, sequence, request );

    *dB = magma_desc_init(
        PlasmaComplexDouble, dA->mb, dA->nb, dA->bsiz,
        dA->m, dA->n, dA->i, dA->j, dA->m, dA->n);

    dB->get_blkaddr = magma_getaddr_cm;
    dB->get_blkldd  = magma_get_blkldd_cm;
    dB->styp = PlasmaCM;

    for (m = 0; m < dA->mt; m++)
    {
        for (n = 0; n < dA->nt; n++)
        {
            X1 = n == 0 ? dA->j%dA->nb : 0;
            Y1 = m == 0 ? dA->i%dA->mb : 0;
            X2 = n == dA->nt-1 ? (dA->j+dA->n-1)%dA->nb+1 : dA->nb;
            Y2 = m == dA->mt-1 ? (dA->i+dA->m-1)%dA->mb+1 : dA->mb;

            MORSE_zlacpy(
                &options,
                PlasmaUpperLower, (Y2-Y1), (X2-X1),
                B(m, n), A(m, n));
        }
    }
}

#if 0
/***************************************************************************//**
 *  Zeroes a submatrix in tile layout - dynamic scheduling
 **/
void magma_pztile_zero(magma_desc_t *dA, magma_sequence_t *sequence, magma_request_t *request)
{
    magma_context_t *magma;
    MorseOption_t options;
    PLASMA_Complex64_t *bdl;
    int X1, Y1;
    int X2, Y2;
    int n, m, ldt;

    magma = magma_context_self();
    if (sequence->status != MAGMA_SUCCESS)
        return;

    morse_options_init( &options, magma, sequence, request );

    for (m = 0; m < dA->mt; m++)
    {
        ldt = BLKLDD(dA, m);
        for (n = 0; n < dA->nt; n++)
        {
            X1 = n == 0 ? dA->j%dA->nb : 0;
            Y1 = m == 0 ? dA->i%dA->mb : 0;
            X2 = n == dA->nt-1 ? (dA->j+dA->n-1)%dA->nb+1 : dA->nb;
            Y2 = m == dA->mt-1 ? (dA->i+dA->m-1)%dA->mb+1 : dA->mb;

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
#endif
