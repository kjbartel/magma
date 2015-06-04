/*
    -- MAGMA (version 1.4.0-beta2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       June 2013

       @author Azzam Haidar
*/

#ifndef _MAGMA_THREADSETTING_H_
#define _MAGMA_THREADSETTING_H_

#ifdef __cplusplus
extern "C" {
#endif
/***************************************************************************//**
 *  Internal routines
 **/
void magma_setlapack_multithreads(magma_int_t numthreads);
void magma_setlapack_sequential();
void magma_setlapack_numthreads(magma_int_t numthreads);
magma_int_t magma_get_numthreads();
/***************************************************************************/
#ifdef __cplusplus
}
#endif

#endif
