/*
    -- MAGMA (version 1.3.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2012
*/

/*#include <quark.h>*/

#ifndef _MAGMA_
#define _MAGMA_

/* ------------------------------------------------------------
 * MAGMA Blas Functions
 * --------------------------------------------------------- */
#include "magmablas.h"

/* ------------------------------------------------------------
 * MAGMA functions
 * --------------------------------------------------------- */
#include "magma_z.h"
#include "magma_c.h"
#include "magma_d.h"
#include "magma_s.h"
#include "magma_zc.h"
#include "magma_ds.h"

/* ------------------------------------------------------------
 * MAGMA constants
 * --------------------------------------------------------- */
#define MAGMA_VERSION_MAJOR 1
#define MAGMA_VERSION_MINOR 3
#define MAGMA_VERSION_MICRO 0

#define MagmaNoTrans       'N'
#define MagmaTrans         'T'
#define MagmaConjTrans     'C'

#define MagmaUpper         'U'
#define MagmaLower         'L'
#define MagmaUpperLower    'A'

#define MagmaNonUnit       'N'
#define MagmaUnit          'U'

#define MagmaLeft          'L'
#define MagmaRight         'R'

#define MagmaForward       'F'
#define MagmaBackward      'B'

#define MagmaColumnwise    'C'
#define MagmaRowwise       'R'

#define MagmaNoVectors     'N'
#define MagmaVectors       'V'

#define MagmaNoTransStr    "NonTrans"
#define MagmaTransStr      "Trans"
#define MagmaConjTransStr  "Conj"

#define MagmaUpperStr      "Upper"
#define MagmaLowerStr      "Lower"
#define MagmaUpperLowerStr "All"

#define MagmaNonUnitStr    "NonUnit"
#define MagmaUnitStr       "Unit"

#define MagmaLeftStr       "Left"
#define MagmaRightStr      "Right"

#define MagmaForwardStr    "Forward"
#define MagmaBackwardStr   "Backward"

#define MagmaColumnwiseStr "Columnwise"
#define MagmaRowwiseStr    "Rowwise"

#define MagmaNoVectorsStr  "NoVectors"
#define MagmaVectorsStr    "Vectors"

/* ------------------------------------------------------------
 * Return codes
 * LAPACK argument errors are < 0 but > MAGMA_ERR.
 * MAGMA errors are < MAGMA_ERR.
 * --------------------------------------------------------- */
#define MAGMA_SUCCESS                 0
#define MAGMA_ERR                  -100
#define MAGMA_ERR_NOT_INITIALIZED  -101
#define MAGMA_ERR_REINITIALIZED    -102
#define MAGMA_ERR_NOT_SUPPORTED    -103
#define MAGMA_ERR_ILLEGAL_VALUE    -104
#define MAGMA_ERR_NOT_FOUND        -105
#define MAGMA_ERR_OUT_OF_RESOURCES -106
#define MAGMA_ERR_ALLOCATION       -106
#define MAGMA_ERR_INTERNAL_LIMIT   -107
#define MAGMA_ERR_UNALLOCATED      -108
#define MAGMA_ERR_FILESYSTEM       -109
#define MAGMA_ERR_UNEXPECTED       -110
#define MAGMA_ERR_SEQUENCE_FLUSHED -111
#define MAGMA_ERR_HOST_ALLOC       -112
#define MAGMA_ERR_DEVICE_ALLOC     -113
#define MAGMA_ERR_CUDASTREAM       -114
#define MAGMA_ERR_INVALID_PTR      -115

/* ------------------------------------------------------------
 *   Define new type that will not be changed by the generator
 * --------------------------------------------------------- */
typedef double real_Double_t;

/* ------------------------------------------------------------
 *   Macros to deal with cuda complex
 * --------------------------------------------------------- */
#define MAGMA_Z_SET2REAL(v, t)    {(v).x = (t); (v).y = 0.0;}
#define MAGMA_Z_EQUAL(u,v)        (((u).x == (v).x) && ((u).y == (v).y))
#define MAGMA_Z_DSCALE(v, t, s)   {(v).x = (t).x/(s); (v).y = (t).y/(s);}
#define MAGMA_Z_MAKE(r, i)        make_cuDoubleComplex((r), (i))
#define MAGMA_Z_REAL(a)           cuCreal(a)
#define MAGMA_Z_IMAG(a)           cuCimag(a)
#define MAGMA_Z_ADD(a, b)         cuCadd((a), (b))
#define MAGMA_Z_SUB(a, b)         cuCsub((a), (b))
#define MAGMA_Z_MUL(a, b)         cuCmul((a), (b))
#define MAGMA_Z_DIV(a, b)         cuCdiv((a), (b))
#define MAGMA_Z_ABS(a)            cuCabs((a))
#define MAGMA_Z_CNJG(a)           cuConj(a)
#define MAGMA_Z_NEGATE(a)         make_cuDoubleComplex( -(a).x, -(a).y )
#define MAGMA_Z_ZERO              make_cuDoubleComplex(0.0, 0.0)
#define MAGMA_Z_ONE               make_cuDoubleComplex(1.0, 0.0)
#define MAGMA_Z_HALF              make_cuDoubleComplex(0.5, 0.0)
#define MAGMA_Z_NEG_ONE           make_cuDoubleComplex(-1.0, 0.0)
#define MAGMA_Z_NEG_HALF          make_cuDoubleComplex(-0.5, 0.0)

#define MAGMA_C_SET2REAL(v, t)    {(v).x = (t); (v).y = 0.0;}
#define MAGMA_C_EQUAL(u,v)        (((u).x == (v).x) && ((u).y == (v).y))
#define MAGMA_C_SSCALE(v, t, s)   {(v).x = (t).x/(s); (v).y = (t).y/(s);}
#define MAGMA_C_MAKE(r, i)        make_cuFloatComplex((r), (i))
#define MAGMA_C_REAL(a)           cuCrealf(a)
#define MAGMA_C_IMAG(a)           cuCimagf(a)
#define MAGMA_C_ADD(a, b)         cuCaddf((a), (b))
#define MAGMA_C_SUB(a, b)         cuCsubf((a), (b))
#define MAGMA_C_MUL(a, b)         cuCmulf((a), (b))
#define MAGMA_C_DIV(a, b)         cuCdivf((a), (b))
#define MAGMA_C_ABS(a)            cuCabsf((a))
#define MAGMA_C_CNJG(a)           cuConjf(a)
#define MAGMA_C_NEGATE(a)         make_cuFloatComplex( -(a).x, -(a).y )
#define MAGMA_C_ZERO              make_cuFloatComplex(0.0, 0.0)
#define MAGMA_C_ONE               make_cuFloatComplex(1.0, 0.0)
#define MAGMA_C_HALF              make_cuFloatComplex(0.5, 0.0)
#define MAGMA_C_NEG_ONE           make_cuFloatComplex(-1.0, 0.0)
#define MAGMA_C_NEG_HALF          make_cuFloatComplex(-0.5, 0.0)

#define MAGMA_D_SET2REAL(v, t)    (v) = (t)
#define MAGMA_D_OP_NEG_ASGN(t, z) (t) = -(z)
#define MAGMA_D_EQUAL(u,v)        ((u) == (v))
#define MAGMA_D_DSCALE(v, t, s)   (v) = (t)/(s)
#define MAGMA_D_MAKE(r, i)        (r)
#define MAGMA_D_REAL(a)           (a)
#define MAGMA_D_IMAG(a)           (a)
#define MAGMA_D_ADD(a, b)         ( (a) + (b) )
#define MAGMA_D_SUB(a, b)         ( (a) - (b) )
#define MAGMA_D_MUL(a, b)         ( (a) * (b) )
#define MAGMA_D_DIV(a, b)         ( (a) / (b) )
#define MAGMA_D_ABS(a)            ((a)>0?(a):-(a))
#define MAGMA_D_CNJG(a)           (a)
#define MAGMA_D_NEGATE(a)         (-(a))
#define MAGMA_D_ZERO              (0.0)
#define MAGMA_D_ONE               (1.0)
#define MAGMA_D_HALF              (0.5)
#define MAGMA_D_NEG_ONE           (-1.0)
#define MAGMA_D_NEG_HALF          (-0.5)

#define MAGMA_S_SET2REAL(v, t)    (v) = (t)
#define MAGMA_S_OP_NEG_ASGN(t, z) (t) = -(z)
#define MAGMA_S_EQUAL(u,v)        ((u) == (v))
#define MAGMA_S_SSCALE(v, t, s)   (v) = (t)/(s)
#define MAGMA_S_MAKE(r, i)        (r)
#define MAGMA_S_REAL(a)           (a)
#define MAGMA_S_IMAG(a)           (a)
#define MAGMA_S_ADD(a, b)         ( (a) + (b) )
#define MAGMA_S_SUB(a, b)         ( (a) - (b) )
#define MAGMA_S_MUL(a, b)         ( (a) * (b) )
#define MAGMA_S_DIV(a, b)         ( (a) / (b) )
#define MAGMA_S_ABS(a)            ((a)>0?(a):-(a))
#define MAGMA_S_CNJG(a)           (a)
#define MAGMA_S_NEGATE(a)         (-(a))
#define MAGMA_S_ZERO              (0.0)
#define MAGMA_S_ONE               (1.0)
#define MAGMA_S_HALF              (0.5)
#define MAGMA_S_NEG_ONE           (-1.0)
#define MAGMA_S_NEG_HALF          (-0.5)

#ifndef CBLAS_SADDR
#define CBLAS_SADDR(a)  &(a)
#endif

#ifdef __cplusplus
extern "C" {
#endif

/* ------------------------------------------------------------
 *   -- MAGMA function definitions
 * --------------------------------------------------------- */

// ========================================
// initialization
void magma_init( void );
void magma_finalize( void );
void magma_version( int* major, int* minor, int* micro );

// ========================================
// memory allocation
magma_err_t magma_malloc( magma_devptr *ptrPtr, size_t bytes );
magma_err_t magma_free  ( magma_devptr ptr );

magma_err_t magma_malloc_cpu( void **ptrPtr, size_t bytes );
magma_err_t magma_free_cpu  ( void *ptr );

magma_err_t magma_malloc_pinned( void **ptrPtr, size_t bytes );
magma_err_t magma_free_pinned  ( void *ptr );

// type-safe convenience functions to avoid using (void**) cast and sizeof(...)
// here n is the number of elements (floats, doubles, etc.) not the number of bytes.
static inline magma_err_t magma_smalloc( float           **ptrPtr, size_t n ) { return magma_malloc( (void**) ptrPtr, n*sizeof(float)           ); }
static inline magma_err_t magma_dmalloc( double          **ptrPtr, size_t n ) { return magma_malloc( (void**) ptrPtr, n*sizeof(double)          ); }
static inline magma_err_t magma_cmalloc( cuFloatComplex  **ptrPtr, size_t n ) { return magma_malloc( (void**) ptrPtr, n*sizeof(cuFloatComplex)  ); }
static inline magma_err_t magma_zmalloc( cuDoubleComplex **ptrPtr, size_t n ) { return magma_malloc( (void**) ptrPtr, n*sizeof(cuDoubleComplex) ); }

static inline magma_err_t magma_smalloc_cpu( float           **ptrPtr, size_t n ) { return magma_malloc_cpu( (void**) ptrPtr, n*sizeof(float)           ); }
static inline magma_err_t magma_dmalloc_cpu( double          **ptrPtr, size_t n ) { return magma_malloc_cpu( (void**) ptrPtr, n*sizeof(double)          ); }
static inline magma_err_t magma_cmalloc_cpu( cuFloatComplex  **ptrPtr, size_t n ) { return magma_malloc_cpu( (void**) ptrPtr, n*sizeof(cuFloatComplex)  ); }
static inline magma_err_t magma_zmalloc_cpu( cuDoubleComplex **ptrPtr, size_t n ) { return magma_malloc_cpu( (void**) ptrPtr, n*sizeof(cuDoubleComplex) ); }

static inline magma_err_t magma_smalloc_pinned( float           **ptrPtr, size_t n ) { return magma_malloc_pinned( (void**) ptrPtr, n*sizeof(float)           ); }
static inline magma_err_t magma_dmalloc_pinned( double          **ptrPtr, size_t n ) { return magma_malloc_pinned( (void**) ptrPtr, n*sizeof(double)          ); }
static inline magma_err_t magma_cmalloc_pinned( cuFloatComplex  **ptrPtr, size_t n ) { return magma_malloc_pinned( (void**) ptrPtr, n*sizeof(cuFloatComplex)  ); }
static inline magma_err_t magma_zmalloc_pinned( cuDoubleComplex **ptrPtr, size_t n ) { return magma_malloc_pinned( (void**) ptrPtr, n*sizeof(cuDoubleComplex) ); }


// ========================================
// device & queue support
void magma_getdevices(
    magma_device_t* devices,
    magma_int_t     size,
    magma_int_t*    numPtr );

void magma_getdevice( magma_device_t* dev );

void magma_setdevice( magma_device_t dev );

void magma_device_sync();


// ========================================
// queue support
void magma_queue_create( /*magma_device_t device,*/ magma_queue_t* queuePtr );

void magma_queue_destroy( magma_queue_t queue );

void magma_queue_sync( magma_queue_t queue );


// ========================================
// event support
void magma_event_create( magma_event_t* eventPtr );

void magma_event_destroy( magma_event_t event );

void magma_event_record( magma_event_t event, magma_queue_t queue );

// blocks CPU until event occurs
void magma_event_sync( magma_event_t event );

// blocks queue (but not CPU) until event occurs
void magma_queue_wait_event( magma_queue_t queue, magma_event_t event );

// magma GPU-complex connection
int magma_buildconnection_mgpu(  magma_int_t gnode[MagmaMaxGPUs+2][MagmaMaxGPUs+2], magma_int_t *nbcmplx, magma_int_t ngpu);



// ========================================
// generic, type-independent routines to copy data.
// type-safe versions which avoid the user needing sizeof(...) are in [sdcz]set_get.cpp
void magma_setvector(
    magma_int_t n, size_t elemSize,
    void const *hx_src, magma_int_t incx,
    void       *dy_dst, magma_int_t incy );

void magma_getvector(
    magma_int_t n, size_t elemSize,
    void const *dx_src, magma_int_t incx,
    void       *hy_dst, magma_int_t incy );

void magma_setvector_async(
    magma_int_t n, size_t elemSize,
    void const *hx_src, magma_int_t incx,
    void       *dy_dst, magma_int_t incy,
    magma_stream_t stream );

void magma_getvector_async(
    magma_int_t n, size_t elemSize,
    void const *dx_src, magma_int_t incx,
    void       *hy_dst, magma_int_t incy,
    magma_stream_t stream );


// ========================================
// error handler
void magma_xerbla( const char *name, magma_int_t info );

/* ------------------------------------------------------------
 *   -- MAGMA Auxiliary structures and functions
 * --------------------------------------------------------- */
typedef struct magma_timestr_s
{
  unsigned int sec;
  unsigned int usec;
} magma_timestr_t;

magma_timestr_t get_current_time(void);
double GetTimerValue(magma_timestr_t time_1, magma_timestr_t time_2);
void printout_devices();
void swp2pswp(char trans, magma_int_t n, magma_int_t *ipiv, magma_int_t *newipiv);
float getv(float *da);

double magma_wtime( void );
double magma_sync_wtime( magma_queue_t queue );
size_t magma_strlcpy(char *dst, const char *src, size_t siz);
int magma_num_gpus( void );
int magma_is_devptr( const void* A );

#ifdef __cplusplus
}
#endif

#endif
