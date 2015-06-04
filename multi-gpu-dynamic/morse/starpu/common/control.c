/**
 *
 * @file control.c
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
#include <stdio.h>
#include <stdlib.h>
#include "morse_starpu.h"

static void morse_starpu_init_plasma_on_worker(void *_magma)
{
    magma_context_t *magma = (magma_context_t*)_magma;
    int i;
    int id = starpu_worker_get_id();
    int nthreads = magma->nthreads_per_worker;
    int bindtab[nthreads];

    for (i = 0; i < nthreads - 1; i++)
    {
        bindtab[i] = i + id*nthreads;
    }

    PLASMA_Init_Affinity(nthreads - 1, bindtab);

    PLASMA_Set(PLASMA_SCHEDULING_MODE, PLASMA_STATIC_SCHEDULING);
    /* PLASMA_Set(PLASMA_SCHEDULING_MODE, PLASMA_DYNAMIC_SCHEDULING); */
}

static void morse_starpu_finalize_plasma_on_worker(void *arg __attribute__((unused)))
{
    PLASMA_Finalize();
}


/***************************************************************************//**
 *  Busy-waiting barrier
 **/
void morse_barrier( magma_context_t *magma )
{
    (void)magma;
#if defined(MORSE_USE_MPI)
    starpu_mpi_barrier(MPI_COMM_WORLD);
#else
    starpu_task_wait_for_all();
#endif
}

int morse_init_scheduler( magma_context_t *magma, int nworkers, int ncudas, int nthreads_per_worker)
{
    int hres = -1;
    if ((nworkers == -1)||(nthreads_per_worker == -1))
    {
        magma->parallel_enabled = MAGMA_FALSE;
        magma->schedopt.starpu->ncpus = nworkers;
        magma->schedopt.starpu->ncuda = ncudas;

        hres = starpu_init( magma->schedopt.starpu );
    }
    else {
        int worker;

        magma->parallel_enabled = MAGMA_TRUE;
        magma->schedopt.starpu->ncpus = nworkers;
        magma->schedopt.starpu->ncuda = ncudas;

        for (worker = 0; worker < nworkers; worker++)
            magma->schedopt.starpu->workers_bindid[worker] = (worker+1)*nthreads_per_worker - 1;

        for (worker = 0; worker < nworkers; worker++)
            magma->schedopt.starpu->workers_bindid[worker + ncudas] = worker*nthreads_per_worker;

        magma->schedopt.starpu->use_explicit_workers_bindid = 1;

        hres = starpu_init( magma->schedopt.starpu );

        magma->nworkers = nworkers;
        magma->nthreads_per_worker = nthreads_per_worker;

        starpu_execute_on_each_worker(morse_starpu_init_plasma_on_worker, &magma, STARPU_CPU);
    }

#if defined(MORSE_USE_MPI)
    starpu_mpi_initialize_extended(&(magma->my_mpi_rank), &(magma->mpi_comm_size));
#endif // MORSE_USE_MPI

    return hres;
}

/***************************************************************************//**
 *
 */
void morse_finalize_scheduler( magma_context_t *magma )
{
    if ( magma->parallel_enabled )
        starpu_execute_on_each_worker(morse_starpu_finalize_plasma_on_worker, NULL, STARPU_CPU);

    starpu_shutdown();
    return;
}

