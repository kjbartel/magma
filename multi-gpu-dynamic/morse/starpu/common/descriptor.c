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

#if defined(MORSE_USE_MPI)

/* Take 24 bits for the data identifier, and 7 bits for data name */
#define TAG_WIDTH_MIN 20
static int tag_width = 31;
static int tag_sep   = 24;

/* Returns the MPI node number where data indexes index is */
static int morse_mpi_mapping(int m, int n, int morse_mpi_comm_size) {
  /* TODO: more clever mapping */
  return (m+n) % morse_mpi_comm_size;
}
#endif

void morse_desc_create( magma_desc_t *desc )
{
    int64_t lm  = desc->lm;
    int64_t ln  = desc->ln;
    int64_t mb  = desc->mb;
    int64_t nb  = desc->nb;
    int64_t lmt = desc->lmt;
    int64_t lnt = desc->lnt;
    int64_t eltsze = plasma_element_size(desc->dtyp);
    int64_t m, n, tempmm, tempnn;
    int64_t block_ind = 0;
    starpu_data_handle_t *tiles;

    desc->occurences = 1;

#if defined(MORSE_USE_CUDA) && 0
    /* Register the matrix as pinned memory */
    if ( cudaHostRegister( desc->mat, lm*ln*eltsze, cudaHostRegisterPortable ) != cudaSuccess )
    {
        magma_warning("morse_desc_create(StarPU)", "cudaHostRegister failed to register the matrix as pinned memory");
    }
#endif

    tiles = (starpu_data_handle_t*)malloc(lnt*lmt*sizeof(starpu_data_handle_t));
    STARPU_ASSERT(tiles);

    desc->schedopt.starpu_handles = tiles;

#if defined(MORSE_USE_MPI)
    int mysize, myrank = desc->myrank;
    magma_context_t *magma;
    static unsigned num_descr;
    unsigned mynum;
    int tag_up = 0, ok;

    magma = magma_context_self();
    if (magma == NULL) {
        magma_error("morse_desc_create", "MAGMA not initialized");
        return;
    }
    mysize = MAGMA_MPI_SIZE;

    MPI_Comm_get_attr(MPI_COMM_WORLD, MPI_TAG_UB, &tag_up, &ok);

    if ( !ok ) {
        magma_error("morse_desc_create", "MPI_TAG_UB not known by MPI");
    }

    while ( ((int)(1UL<<tag_width)) >= tag_up 
            && (tag_width >= TAG_WIDTH_MIN) ) {
      tag_width--;
      tag_sep--;
    }

    if ( tag_width < TAG_WIDTH_MIN ) {
        magma_error("morse_desc_create", "MPI_TAG_UB too small to identify all the data");
        return;
    }

    if ( lnt*lmt > ((int)(1UL<<tag_sep)) ) {
        magma_error("morse_desc_create", "Too many tiles in the descriptor for MPI tags");
        return;
    }

    /* TODO: manage to reuse MPI tag space */
    mynum = num_descr++;
    if (mynum >= 1UL<<(tag_width-tag_sep)) {
        magma_error("morse_desc_create", "Number of descriptor available in MPI mode out of stock");
        return;
    }

    for (n = 0; n < lnt; n++) {
        tempnn = n == lnt-1 ? ln - n * nb : nb;

        for (m = 0; m < lmt; m++) {

            int owner = morse_mpi_mapping(m, n, mysize);
            tempmm = m == lmt-1 ? lm - m * mb : mb;

            if ( myrank == owner 
                    || myrank == 0 /* XXX: For checks on node 0 */
                    ) {
                starpu_matrix_data_register(&tiles[block_ind], 0,
                                            (uintptr_t)desc->get_blkaddr(desc, m, n),
                                            BLKLDD(desc, m), tempmm, tempnn, eltsze);
            }
            else {
                starpu_matrix_data_register(&tiles[block_ind], -1,
                                            (uintptr_t) NULL,
                                            BLKLDD(desc, m), tempmm, tempnn, eltsze);
            }
            if (tiles[block_ind])
            {
                starpu_data_set_rank(tiles[block_ind], owner);
                starpu_data_set_tag(tiles[block_ind], (mynum<<tag_sep) | (m*lnt+n));
            }
            /* splagma_set_reduction_methods(tiles[block_ind], desc->dtyp); */
            block_ind++;
        }
    }
#else
    for (n = 0; n < lnt; n++) {
        tempnn = n == lnt-1 ? ln - n * nb : nb;
        for (m = 0; m < lmt; m++) {
            tempmm = m == lmt-1 ? lm - m * mb : mb;
                starpu_matrix_data_register(&tiles[block_ind], 0,
                                            (uintptr_t)desc->get_blkaddr(desc, m, n),
                                            BLKLDD(desc, m), tempmm, tempnn, eltsze);
            /* magma_set_reduction_methods(tiles[block_ind], desc->dtyp); */
            block_ind++;
        }
    }
#endif
}

void morse_desc_destroy( magma_desc_t *desc )
{
    int lmt = desc->lmt;
    int lnt = desc->lnt;

    starpu_data_handle_t *tiles = desc->schedopt.starpu_handles;

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
  (void)desc;
  return;
}

void morse_desc_submatrix( magma_desc_t *desc )
{
    desc->occurences++;
    return;
}


#if defined(MORSE_USE_MPI)
int morse_desc_getoncpu( magma_desc_t *desc )
{
    int lmt = desc->lmt;
    int lnt = desc->lnt;
    int  m, n;
    int64_t block_ind = 0;

    for (n = 0; n < lnt; n++)
    for (m = 0; m < lmt; m++)
    {
        /* "CPU" in MPI means node 0.  First make sure we have a copy on node 0.  */
        starpu_mpi_get_data_on_node_detached(MPI_COMM_WORLD, desc->schedopt.starpu_handles[block_ind], 0, NULL, NULL);

        /* And make sure it's not in on some GPU on node 0.  */
        if (MAGMA_my_mpi_rank() == 0) {
            starpu_data_acquire(desc->schedopt.starpu_handles[block_ind], STARPU_R);
            starpu_data_release(desc->schedopt.starpu_handles[block_ind]);
        }
        block_ind++;
    }
    return MAGMA_SUCCESS;
}

#else
int morse_desc_acquire( magma_desc_t *desc )
{
    int lmt = desc->lmt;
    int lnt = desc->lnt;
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
    int lmt = desc->lmt;
    int lnt = desc->lnt;
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
#endif

void *morse_desc_getaddr( magma_desc_t *desc, int m, int n )
{
    return (void *)(desc->schedopt.starpu_handles[(int64_t)(desc->lmt) * (int64_t)n + (int64_t)m ]);
}
