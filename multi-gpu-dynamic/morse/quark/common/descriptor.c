/**
 *
 * @file descriptor.c
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 
 * @author Vijay Joshi
 * @date 2011-10-29
 *
 **/
#include <stdlib.h>
#include "morse_quark.h"

void morse_desc_create( magma_desc_t *desc )
{
    return;
}

void morse_desc_destroy( magma_desc_t *desc )
{
    return;
}

void morse_desc_init( magma_desc_t *desc )
{
    return;
}

void morse_desc_submatrix( magma_desc_t *desc )
{
    return;
}


int morse_desc_acquire( magma_desc_t *desc )
{
    return MAGMA_SUCCESS;
}

int morse_desc_release( magma_desc_t *desc )
{
    return MAGMA_SUCCESS;
}

void *morse_desc_getaddr( magma_desc_t *desc, int m, int n )
{
    return magma_getaddr( desc, m, n );
}
