/**
 *
 * @file context.h
 *
 *  PLAGMA auxiliary routines
 *  PLAGMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.3.1
 * @author Cedric Augonnet
 * @author Mathieu Faverge
 * @author Jakub Kurzak
 * @date 2010-11-15
 *
 **/
#ifndef _MAGMA_CONTEXT_H_
#define _MAGMA_CONTEXT_H_

enum sched_e {
  MAGMA_SCHED_QUARK,
  MAGMA_SCHED_STARPU,
};
typedef enum sched_e sched_t;


/* Quark */
#if !defined(QUARK_H)
struct quark_s;
typedef struct quark_s Quark;
#endif

/* StarPU */
#if !defined(__STARPU_H__)
struct starpu_conf;
#endif
typedef struct starpu_conf Starpu_Conf_t; 

/***************************************************************************//**
 *  PLAGMA context
 **/
struct magma_context_s {
    sched_t scheduler;
    int nworkers;
    int ncudas;
    int nthreads_per_worker;
    int world_size;
    int group_size;

    /* Boolean flags */
    MAGMA_bool errors_enabled;
    MAGMA_bool warnings_enabled;
    MAGMA_bool autotuning_enabled;
    MAGMA_bool parallel_enabled;
    MAGMA_bool profiling_enabled;

    MAGMA_enum householder;    // "domino" (flat) or tree-based (reduction) Householder
    MAGMA_enum translation;    // In place or Out of place layout conversion

    int nb;
    int ib;
    int rhblock;

    union {
        Quark         *quark;
        Starpu_Conf_t *starpu;
    } schedopt;

};

#ifdef __cplusplus
extern "C" {
#endif

/***************************************************************************//**
 *  Internal routines
 **/
magma_context_t *magma_context_create();
magma_context_t *magma_context_self();

int   magma_context_destroy();

#ifdef __cplusplus
}
#endif

#endif
