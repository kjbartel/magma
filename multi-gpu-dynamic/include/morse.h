/**
 *
 *  @file morse.h
 *
 *  MAGMA codelets kernel
 *  MAGMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver,
 *  and INRIA Bordeaux Sud-Ouest
 *
 *  @version 2.3.1
 *  @author Mathieu Faverge
 *  @author Cedric Augonnet
 *  @date 2011-06-01
 *
 **/
#ifndef _MORSE_H_
#define _MORSE_H_

#if defined( _WIN32 )
  /* This must be included before INPUT is defined below, otherwise we
     have a name clash/problem  */
  #include <windows.h>
  #include <limits.h>
#else
  #include <inttypes.h>
#endif

#include "morse_kernels.h"

typedef struct MorseOption_s {
    magma_sequence_t *sequence;
    magma_request_t  *request;
    int    profiling;
    int    parallel;
    int    priority; 
    Quark *quark;
    Quark_Task_Flags *task_flags;
    int    nb;
    size_t ws_hsize;
    size_t ws_dsize;
    void  *ws_host;
    void  *ws_device;
} MorseOption_t;

/*********************
 * Context
 */
int morse_init_scheduler    ( magma_context_t*, int, int, int );
void morse_finalize_scheduler( magma_context_t* );

void morse_context_create ( magma_context_t* );
void morse_context_destroy( magma_context_t* );

void morse_enable ( MAGMA_enum );
void morse_disable( MAGMA_enum );

void morse_barrier( magma_context_t* );

/*********************
 * Descriptor
 */
void morse_desc_create   ( magma_desc_t* );
void morse_desc_destroy  ( magma_desc_t* );
void morse_desc_init     ( magma_desc_t* );
void morse_desc_submatrix( magma_desc_t* );
void*morse_desc_getaddr  ( magma_desc_t*, int, int );

#if defined(MORSE_USE_MPI)
int  morse_desc_getoncpu  ( magma_desc_t* );
#else
int  morse_desc_acquire  ( magma_desc_t* );
int  morse_desc_release  ( magma_desc_t* );
#define morse_desc_getoncpu( desc ) { morse_desc_acquire( desc ); morse_desc_release( desc ); }
#endif

/*********************
 * Async 
 */
int  morse_sequence_create ( magma_context_t*, magma_sequence_t* );
int  morse_sequence_destroy( magma_context_t*, magma_sequence_t* );
int  morse_sequence_wait   ( magma_context_t*, magma_sequence_t* );
void morse_sequence_flush  ( void*, magma_sequence_t*, magma_request_t*, int);

/*********************
 * Kernels options
 */
#define MORSE_PRIORITY_MIN  0
#define MORSE_PRIORITY_MAX  INT_MAX

void  morse_options_init        ( MorseOption_t*, magma_context_t*, magma_sequence_t*, magma_request_t* );
void  morse_options_finalize    ( MorseOption_t*, magma_context_t* );
int   morse_options_ws_alloc    ( MorseOption_t*, size_t, size_t );
int   morse_options_ws_free     ( MorseOption_t* );
/* void *morse_options_ws_gethost  ( MorseOption_t* ); */
/* void *morse_options_ws_getdevice( MorseOption_t* ); */

void  morse_schedprofile_display(void);

#ifdef MORSE_USE_CUDA
void  morse_zlocality_allrestrict( uint32_t );
void  morse_zlocality_onerestrict( morse_kernel_t, uint32_t );
void  morse_zlocality_onerestore ( morse_kernel_t );
#else
void  morse_zlocality_allrestrict();
void  morse_zlocality_onerestrict();
void  morse_zlocality_onerestore( );
#endif
void  morse_zlocality_allrestore ( );
void  morse_zdisplay_allprofile  ( );
void  morse_zdisplay_oneprofile  ( morse_kernel_t );

#ifdef MORSE_USE_CUDA
void  morse_clocality_allrestrict( uint32_t );
void  morse_clocality_onerestrict( morse_kernel_t, uint32_t );
void  morse_clocality_onerestore ( morse_kernel_t );
#else
void  morse_clocality_allrestrict();
void  morse_clocality_onerestrict();
void  morse_clocality_onerestore( );
#endif
void  morse_clocality_allrestore ( );
void  morse_cdisplay_allprofile  ( );
void  morse_cdisplay_oneprofile  ( morse_kernel_t );

#ifdef MORSE_USE_CUDA
void  morse_dlocality_allrestrict( uint32_t );
void  morse_dlocality_onerestrict( morse_kernel_t, uint32_t );
void  morse_dlocality_onerestore ( morse_kernel_t );
#else
void  morse_dlocality_allrestrict();
void  morse_dlocality_onerestrict();
void  morse_dlocality_onerestore( );
#endif
void  morse_dlocality_allrestore ( );
void  morse_ddisplay_allprofile  ( );
void  morse_ddisplay_oneprofile  ( morse_kernel_t );

#ifdef MORSE_USE_CUDA
void  morse_slocality_allrestrict( uint32_t );
void  morse_slocality_onerestrict( morse_kernel_t, uint32_t );
void  morse_slocality_onerestore ( morse_kernel_t );
#else
void  morse_slocality_allrestrict();
void  morse_slocality_onerestrict();
void  morse_slocality_onerestore( );
#endif
void  morse_slocality_allrestore ( );
void  morse_sdisplay_allprofile  ( );
void  morse_sdisplay_oneprofile  ( morse_kernel_t );

/*********************
 * Kernels
 */
#include "morse_z.h"
#include "morse_d.h"
#include "morse_c.h"
#include "morse_s.h"

#endif
