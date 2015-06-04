/**
 *
 * @file ztrtri_extra.c
 *
 *  PLASMA computational routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *	TODO	
 * @version 2.3.1
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @precisions normal z -> s d c
 *
 **/
#include "common.h"

void SPLAGMA_locality_restrict_ztrtri(uint32_t trsm, uint32_t gemm, uint32_t trtri )
{
#ifdef MORSE_USE_CUDA
    cl_ztrsm_restrict_where(trsm);
    cl_zgemm_restrict_where(gemm);
    cl_zherk_restrict_where(trtri);
#endif
}

void SPLAGMA_locality_restore_ztrtri(void)
{
#ifdef MORSE_USE_CUDA
    cl_ztrsm_restore_where();
    cl_zgemm_restore_where();
    cl_ztrtri_restore_where();
#endif
}

void SPLAGMA_profile_ztrtri(void)
{
    profiling_display_ztrsm_info();
    profiling_display_zgemm_info();
    profiling_display_ztrtri_info();
}

