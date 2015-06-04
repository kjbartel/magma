/**
 *
 * @file zpotrf_extra.c
 *
 *  PLASMA computational routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.3.1
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @precisions normal z -> s d c
 *
 **/
#include "common.h"

void MAGMA_locality_restrict_zpotrf(uint32_t potrf, uint32_t trsm, uint32_t herk, uint32_t gemm )
{
#ifdef MORSE_USE_CUDA
    cl_zpotrf_restrict_where(potrf);
    cl_ztrsm_restrict_where(trsm);
    cl_zgemm_restrict_where(gemm);
    cl_zherk_restrict_where(herk);
#endif
}

void MAGMA_locality_restore_zpotrf(void)
{
#ifdef MORSE_USE_CUDA
    cl_zpotrf_restore_where();
    cl_ztrsm_restore_where();
    cl_zgemm_restore_where();
    cl_zherk_restore_where();
#endif
}

void MAGMA_profile_zpotrf(void)
{
#ifdef MORSE_SCHEDULER_QUARK
#else
    profiling_display_zpotrf_info();
    profiling_display_ztrsm_info();
    profiling_display_zherk_info();
    profiling_display_zgemm_info();
#endif
}

/* void MAGMA_cload_potrf_FakeModel(void) */
/* { */
 /*     cl_cpotrf_load_fake_model(); */
/*     cl_cherk_load_fake_model(); */
/*     cl_ctrsm_load_fake_model(); */
/*     cl_cgemm_load_fake_model(); */
/* } */

/* void MAGMA_crestore_potrf_Model(void) */
/* { */
/*     cl_cpotrf_restore_model(); */
/*     cl_cherk_restore_model(); */
/*     cl_ctrsm_restore_model(); */
/*     cl_cgemm_restore_model(); */
/* } */
