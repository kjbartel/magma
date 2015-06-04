/*
 *   -- MAGMA (version 1.2.1) --
 *      Univ. of Tennessee, Knoxville
 *      Univ. of California, Berkeley
 *      Univ. of Colorado, Denver
 *      June 2012
 *
 * @generated c Thu Jun 28 12:30:03 2012
 */

#ifndef _MAGMA_CBULGEINC_H_
#define _MAGMA_CBULGEINC_H_

#define PRECISION_c
#ifdef __cplusplus
extern "C" {
#endif


/***************************************************************************//**
 *  Configuration
 **/

 // maximum contexts
#define MAX_THREADS_BLG         256

real_Double_t get_time_azz(void);
void findVTpos(magma_int_t N, magma_int_t NB, magma_int_t Vblksiz, magma_int_t sweep, magma_int_t st, magma_int_t *Vpos, magma_int_t *TAUpos, magma_int_t *Tpos, magma_int_t *myblkid);
void findVTsiz(magma_int_t N, magma_int_t NB, magma_int_t Vblksiz, magma_int_t *blkcnt, magma_int_t *LDV);
magma_int_t plasma_ceildiv(magma_int_t a, magma_int_t b);


/*
extern volatile magma_int_t barrier_in[MAX_THREADS_BLG];
extern volatile magma_int_t barrier_out[MAX_THREADS_BLG];
extern volatile magma_int_t *ss_prog;
*/

 /***************************************************************************//**
 *  Static scheduler
 **/
/*
#define ssched_init(nbtiles) \
{ \
        volatile int   prog_ol[2*nbtiles+10];\
                 int   iamdone[MAX_THREADS_BLG]; \
                 int   thread_num[MAX_THREADS_BLG];\
        pthread_t      thread_id[MAX_THREADS_BLG];\
        pthread_attr_t thread_attr;\
}
*/
////////////////////////////////////////////////////////////////////////////////////////////////////



 struct gbstrct_blg {
    cuFloatComplex *dQ1;
    cuFloatComplex *dT1;
    cuFloatComplex *dT2;
    cuFloatComplex *dV2;
    cuFloatComplex *dE;
    cuFloatComplex *T;
    cuFloatComplex *A;
    cuFloatComplex *V;
    cuFloatComplex *TAU;
    cuFloatComplex *E;
    cuFloatComplex *E_CPU;
    int cores_num;
    int locores_num;
    int overlapQ1;
    int usemulticpu;
    int NB;
    int NBTILES;
    int N;
    int NE;
    int N_CPU;
    int N_GPU;
    int LDA;
    int LDE;
    int BAND;
    int grsiz;
    int Vblksiz;
    int WANTZ;
    char SIDE;
    real_Double_t *timeblg;
    real_Double_t *timeaplQ;
    volatile int *ss_prog;
} ;

// declare globals here; defined in chetrd_bhe2trc.cpp
extern struct gbstrct_blg core_in_all;





////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////
#define MAX_EVENTSBLG 163840
//#define MAX_EVENTSBLG 1048576

// declare globals here; defined in chetrd_bhe2trc.cpp
extern int           event_numblg        [MAX_THREADS_BLG]                 __attribute__ ((aligned (128)));
extern real_Double_t event_start_timeblg [MAX_THREADS_BLG]                 __attribute__ ((aligned (128)));
extern real_Double_t event_end_timeblg   [MAX_THREADS_BLG]                 __attribute__ ((aligned (128)));
extern real_Double_t event_logblg        [MAX_THREADS_BLG][MAX_EVENTSBLG]  __attribute__ ((aligned (128)));
extern int           log_eventsblg;

#define core_event_startblg(my_core_id)\
    event_start_timeblg[my_core_id] = get_time_azz();\

#define core_event_endblg(my_core_id)\
    event_end_timeblg[my_core_id] = get_time_azz();\

#define core_log_eventblg(event, my_core_id)\
    event_logblg[my_core_id][event_numblg[my_core_id]+0] = my_core_id;\
    event_logblg[my_core_id][event_numblg[my_core_id]+1] = event_start_timeblg[my_core_id];\
    event_logblg[my_core_id][event_numblg[my_core_id]+2] = event_end_timeblg[my_core_id];\
    event_logblg[my_core_id][event_numblg[my_core_id]+3] = (event);\
    event_numblg[my_core_id] += (log_eventsblg << 2);\
    event_numblg[my_core_id] &= (MAX_EVENTSBLG-1);

////////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef __cplusplus
}
#endif

#undef PRECISION_c
#endif




