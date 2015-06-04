/**
 *
 * @file descriptor.c
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 1.1.0
 * @author Mathieu Faverge
 * @author Cedric Augonnet
 * @date 2010-11-15
 *
 **/
#include <stdlib.h>
#include "morse_starpu.h"

void morse_desc_create( magma_desc_t *desc )
{
    PLASMA_desc *pdesc = &(desc->desc);
    int64_t lm  = pdesc->lm;
    int64_t ln  = pdesc->ln;
    int64_t mb  = pdesc->mb;
    int64_t nb  = pdesc->nb;
    int64_t lmt = pdesc->lmt;
    int64_t lnt = pdesc->lnt;
    int64_t eltsze = plasma_element_size(pdesc->dtyp);
    int64_t m, n, tempmm, tempnn;
    int64_t block_ind = 0;
    starpu_data_handle *tiles;

    desc->occurences = 1;

#if defined(MORSE_USE_CUDA)
    /* Register the matrix as pinned memory */
/*     if ( cudaHostRegister( pdesc->mat, lm*ln*eltsze, cudaHostRegisterPortable ) != cudaSuccess )  */
/*     { */
/*         magma_warning("morse_desc_create(StarPU)", "cudaHostRegister failed to register the matrix as pinned memory"); */
/*     } */
#endif

    tiles = (starpu_data_handle*)malloc(lnt*lmt*sizeof(starpu_data_handle));
    STARPU_ASSERT(tiles);
    
    desc->schedopt.starpu_handles = tiles;
    
    for (n = 0; n < lnt; n++) {
        tempnn = n == lnt-1 ? ln - n * nb : nb;
        for (m = 0; m < lmt; m++) {
            tempmm = m == lmt-1 ? lm - m * mb : mb;
            starpu_matrix_data_register(&tiles[block_ind], 0,
                                        (uintptr_t)magma_getaddr(desc, m, n),
                                        BLKLDD(desc, m), tempmm, tempnn, eltsze);
            /* magma_set_reduction_methods(tiles[block_ind], pdesc->dtyp); */
            block_ind++;
        }
    }
}

void morse_desc_destroy( magma_desc_t *desc )
{
    int lmt = desc->desc.lmt;
    int lnt = desc->desc.lnt;
    
    starpu_data_handle *tiles = desc->schedopt.starpu_handles;
    
    int  m, n;
    int64_t block_ind = 0;

    desc->occurences--;
    if ( desc->occurences == 0 ) {
        for (n = 0; n < lnt; n++)
            for (m = 0; m < lmt; m++) {
                starpu_data_unregister(tiles[block_ind]);
                block_ind++;
            }
        free(tiles);
    }
}

void morse_desc_init( magma_desc_t *desc )
{
  return;
}

void morse_desc_submatrix( magma_desc_t *desc )
{
    desc->occurences++;
    return;
}


int morse_desc_acquire( magma_desc_t *desc )
{
    int lmt = desc->desc.lmt;
    int lnt = desc->desc.lnt;
    int  m, n;
    int64_t block_ind = 0;

    for (n = 0; n < lnt; n++)
    for (m = 0; m < lmt; m++)
    {
            starpu_data_acquire(desc->schedopt.starpu_handles[block_ind], STARPU_R);
        block_ind++;
    }
    return MAGMA_SUCCESS;
}

int morse_desc_release( magma_desc_t *desc )
{
    int lmt = desc->desc.lmt;
    int lnt = desc->desc.lnt;
    int  m, n;
    int64_t block_ind = 0;
    for (n = 0; n < lnt; n++)
    for (m = 0; m < lmt; m++)
    {
        starpu_data_release(desc->schedopt.starpu_handles[block_ind]);
        block_ind++;
    }
    return MAGMA_SUCCESS;
}

void *morse_desc_getaddr( magma_desc_t *desc, int m, int n )
{
    return (void *)(desc->schedopt.starpu_handles[(int64_t)(desc->desc.lmt) * (int64_t)n + (int64_t)m ]);
}
