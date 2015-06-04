/**	
 *
 * @file zpotri_extra.c
 *
 *  PLASMA computational routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 * //TODO
 * @version 2.3.1
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @precisions normal z -> s d c
 *
 **/
#include "common.h"

void SPLAGMA_locality_restrict_zpotri(uint32_t trtri, uint32_t lauum)
{
#ifdef MORSE_USE_CUDA
    cl_ztrtri_restrict_where(trtri);
    cl_zlauum_restrict_where(lauum);
#endif
}

void SPLAGMA_locality_restore_zpotri(void)
{
#ifdef MORSE_USE_CUDA
    cl_ztrtri_restore_where();
    cl_zlauum_restore_where();
#endif
}

void SPLAGMA_profile_zpotri(void)
{
    profiling_display_ztrtri_info();
    profiling_display_zlauum_info();
}


