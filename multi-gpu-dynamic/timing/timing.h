#ifndef TIMING_H
#define TIMING_H

typedef double morse_time_t;

enum iparam_timing {
    IPARAM_THRDNBR,        /* Number of cores                            */
    IPARAM_THRDNBR_SUBGRP, /* Number of cores in a subgroup (NUMA node)  */
    IPARAM_SCHEDULER,      /* What scheduler do we choose (dyn, stat)    */
    IPARAM_M,              /* Number of rows of the matrix               */
    IPARAM_N,              /* Number of columns of the matrix            */
    IPARAM_K,              /* RHS or K                                   */
    IPARAM_LDA,            /* Leading dimension of A                     */
    IPARAM_LDB,            /* Leading dimension of B                     */
    IPARAM_LDC,            /* Leading dimension of C                     */
    IPARAM_IB,             /* Inner-blocking size                        */
    IPARAM_NB,             /* Number of columns in a tile                */
    IPARAM_MB,             /* Number of rows in a tile                   */
    IPARAM_NITER,          /* Number of iteration of each test           */
    IPARAM_WARMUP,         /* Run one test to load dynamic libraries     */
    IPARAM_CHECK,          /* Checking activated or not                  */
    IPARAM_VERBOSE,        /* How much noise do we want?                 */
    IPARAM_AUTOTUNING,     /* Disable/enable autotuning                  */
    IPARAM_INPUTFMT,       /* Input format (Use only for getmi/gecfi)    */
    IPARAM_OUTPUTFMT,      /* Output format (Use only for getmi/gecfi)   */
    IPARAM_TRACE,          /* Generate trace on the first non warmup run */
    IPARAM_DAG,            /* Do we require to output the DOT file?      */
    IPARAM_ASYNC,          /* Asynchronous calls                         */
    IPARAM_MX,             /* */
    IPARAM_NX,             /* */
    IPARAM_RHBLK,          /* Householder reduction parameter for QR/LQ  */
    IPARAM_INPLACE,        /* InPlace/OutOfPlace translation mode        */

    IPARAM_INVERSE,
    IPARAM_NCUDAS,
    IPARAM_NDOM,
    /* Added for StarPU version */
    IPARAM_PROFILE,
    IPARAM_PEAK,
    IPARAM_PARALLEL_TASKS,
    IPARAM_NO_CPU,
    IPARAM_BOUND,
    IPARAM_BOUNDDEPS,
    IPARAM_BOUNDDEPSPRIO,
    /* End */
    IPARAM_SIZEOF
};

enum dparam_timing {
  IPARAM_TIME,
  IPARAM_ANORM,
  IPARAM_BNORM,
  IPARAM_XNORM,
  IPARAM_RNORM,
  IPARAM_AinvNORM,
  IPARAM_ESTIMATED_PEAK,
  IPARAM_RES,
  /* Begin section for hydra integration tool */
  IPARAM_THRESHOLD_CHECK, /* Maximum value accepted for: |Ax-b||/N/eps/(||A||||x||+||b||) */
  /* End section for hydra integration tool  */
  IPARAM_DNBPARAM
};

#define PASTE_CODE_IPARAM_LOCALS(iparam)           \
    double  t;                                     \
    int64_t M     = iparam[IPARAM_M];              \
    int64_t N     = iparam[IPARAM_N];              \
    int64_t K     = iparam[IPARAM_K];              \
    int64_t NRHS  = K;                             \
    int64_t LDA   = max(M, iparam[IPARAM_LDA]);    \
    int64_t LDB   = max(N, iparam[IPARAM_LDB]);    \
    int64_t LDC   = max(K, iparam[IPARAM_LDC]);    \
    int64_t IB    = iparam[IPARAM_IB];             \
    int64_t MB    = iparam[IPARAM_MB];             \
    int64_t NB    = iparam[IPARAM_NB];             \
    int64_t MT    = (M%MB==0) ? (M/MB) : (M/MB+1); \
    int64_t NT    = (N%NB==0) ? (N/NB) : (N/NB+1); \
    int check     = iparam[IPARAM_CHECK];          \
    int loud      = iparam[IPARAM_VERBOSE];        \
    (void)M;(void)N;(void)K;(void)NRHS;            \
    (void)LDA;(void)LDB;(void)LDC;                 \
    (void)IB;(void)MB;(void)NB;(void)MT;(void)NT;  \
    (void)check;(void)loud;

/* Paste code to allocate a matrix in desc if cond_init is true */
#define PASTE_CODE_ALLOCATE_MATRIX_TILE(_desc_, _cond_, _type_, _type2_, _lda_, _m_, _n_) \
    magma_desc_t *_desc_ = NULL;                                         \
    if( _cond_ ) {                                                      \
        _type_ *ptr = NULL;                                             \
        ptr = (_type_*)malloc( (_lda_) * (_n_) * sizeof(_type_) );      \
        if ( ! ptr ) {                                                  \
            fprintf(stderr, "Our of Memory for %s\n", #_desc_);         \
            return -1;                                                  \
        }                                                               \
        MAGMA_Desc_Create(&(_desc_), ptr, _type2_, MB, NB, MB*NB, _lda_, _n_, 0, 0, _m_, _n_); \
    }

#define PASTE_CODE_FREE_MATRIX(_desc_)                                  \
    if ( _desc_ != NULL ) {                                             \
        free(_desc_->mat);                                              \
    }                                                                   \
    MAGMA_Desc_Destroy( &_desc_ );

#define PASTE_TILE_TO_LAPACK(_desc_, _name_, _cond_, _type_, _lda_, _n_) \
    _type_ *_name_ = NULL;                                              \
    if ( _cond_ ) {                                                     \
        _name_ = (_type_*)malloc( (_lda_) * (_n_) * sizeof(_type_));    \
        if ( ! _name_ ) {                                               \
            fprintf(stderr, "Our of Memory for %s\n", #_name_);         \
            return -1;                                                  \
        }                                                               \
        switch( (_desc_)->dtyp ) {                                      \
        case PlasmaComplexDouble:                                       \
            MAGMA_zTile_to_Lapack(_desc_, (void*)_name_, _lda_);        \
            break;                                                      \
        case PlasmaComplexFloat:                                        \
            MAGMA_cTile_to_Lapack(_desc_, (void*)_name_, _lda_);        \
            break;                                                      \
        case PlasmaRealDouble:                                          \
            MAGMA_dTile_to_Lapack(_desc_, (void*)_name_, _lda_);        \
            break;                                                      \
        case PlasmaRealFloat:                                           \
            MAGMA_sTile_to_Lapack(_desc_, (void*)_name_, _lda_);        \
            break;                                                      \
        }                                                               \
    }

#define PASTE_CODE_ALLOCATE_MATRIX(_name_, _cond_, _type_, _lda_, _n_)  \
    _type_ *_name_ = NULL;                                              \
    if( _cond_ ) {                                                      \
        _name_ = (_type_*)malloc( (_lda_) * (_n_) * sizeof(_type_) );   \
        if ( ! _name_ ) {                                               \
            fprintf(stderr, "Our of Memory for %s\n", #_name_);         \
            return -1;                                                  \
        }                                                               \
    }

#define PASTE_CODE_ALLOCATE_COPY(_name_, _cond_, _type_, _orig_, _lda_, _n_) \
    _type_ *_name_ = NULL;                                              \
    if( _cond_ ) {                                                      \
        _name_ = (_type_*)malloc( (_lda_) * (_n_) * sizeof(_type_) );   \
        if ( ! _name_ ) {                                               \
            fprintf(stderr, "Our of Memory for %s\n", #_name_);         \
            return -1;                                                  \
        }                                                               \
        memcpy(_name_, _orig_, (_lda_) * (_n_) * sizeof(_type_) );      \
    }

/*********************
 *
 * Macro for trace generation
 *
 */
#if defined(MORSE_TRACE)
#define START_TRACING()                         \
    if(iparam[IPARAM_TRACE] == 2)               \
      eztrace_start();

#define STOP_TRACING()                         \
    if(iparam[IPARAM_TRACE] == 2)              \
        eztrace_stop();

#else /* defined(MORSE_TRACE) */

#define START_TRACING()                         \
  if( 0 ) {};

#define STOP_TRACING()                          \
  if( 0 ) {};

#endif /* defined(MORSE_TRACE) */

/*********************
 *
 * Macro for DAG generation
 *
 */
#if 0
#define START_DAG()                   \
    if ( iparam[IPARAM_DAG] == 2 )    \
        MAGMA_Enable(MAGMA_DAG);

#define STOP_DAG()                             \
    if ( iparam[IPARAM_DAG] == 2 )             \
        MAGMA_Disable(MAGMA_DAG);
#else
#define START_DAG()  do {} while(0);
#define STOP_DAG()  do {} while(0);
#endif

/*********************
 *
 * General Macros for timing
 *
 */
#define START_TIMING()                          \
  START_DAG();                                  \
  START_TRACING();                              \
  t = -cWtime();

#define STOP_TIMING()                           \
  t += cWtime();                                \
  STOP_TRACING();                               \
  STOP_DAG();                                   \
  *t_ = t;

#endif /* TIMING_H */
