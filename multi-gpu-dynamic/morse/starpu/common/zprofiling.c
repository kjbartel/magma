/**
 *
 *  @file zprofiling.c
 *
 *  MAGMA codelets kernel
 *  MAGMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver,
 *  and INRIA Bordeaux Sud-Ouest
 *
 *  @version 1.1.0
 *  @author Mathieu Faverge
 *  @author Cedric Augonnet
 *  @date 2011-06-01
 *  @precisions normal z -> s d c
 *
 **/
#include "morse_starpu.h"

void morse_zdisplay_allprofile()
{
    profiling_display_zgemm_info();
    profiling_display_ztrsm_info();
    /* profiling_display_ztrmm_info(); */
    /* profiling_display_zhemm_info(); */
    /* profiling_display_zsymm_info(); */
    profiling_display_zherk_info();
    /* profiling_display_zsyrk_info(); */
    /* profiling_display_zher2k_info(); */
    /* profiling_display_zsyr2k_info(); */

    profiling_display_zlacpy_info();
#if defined(PRECISION_z) || defined(PRECISION_c)
    profiling_display_zplghe_info();
#endif
    profiling_display_zplgsy_info();
    profiling_display_zplrnt_info();
    /* profiling_display_zlange_info(); */
    /* profiling_display_zlanhe_info(); */
    /* profiling_display_zlansy_info(); */

    profiling_display_zpotrf_info();

    profiling_display_zgetrl_info();
    profiling_display_zgessm_info();
    profiling_display_ztstrf_info();
    profiling_display_zssssm_info();

    profiling_display_zgeqrt_info();
    profiling_display_zunmqr_info();
    profiling_display_ztsqrt_info();
    profiling_display_ztsmqr_info();
    /* profiling_display_zttqrt_info(); */
    /* profiling_display_zttmqr_info(); */

    /* profiling_display_zgelqt_info(); */
    /* profiling_display_zunmlq_info(); */
    /* profiling_display_ztslqt_info(); */
    /* profiling_display_ztsmlq_info(); */
    /* profiling_display_zttlqt_info(); */
    /* profiling_display_zttmlq_info(); */
}

void morse_zdisplay_oneprofile( morse_kernel_t kernel )
{
    switch( kernel ) {
    case MORSE_GEMM:  profiling_display_zgemm_info();  break;
    case MORSE_TRSM:  profiling_display_ztrsm_info();  break;
    /* case MORSE_TRMM:  profiling_display_ztrmm_info();  break; */
    /* case MORSE_HEMM:  profiling_display_zhemm_info();  break; */
    /* case MORSE_SYMM:  profiling_display_zsymm_info();  break; */
    case MORSE_HERK:  profiling_display_zherk_info();  break;
    /* case MORSE_SYRK:  profiling_display_zsyrk_info();  break; */
    /* case MORSE_HER2K: profiling_display_zher2k_info(); break; */
    /* case MORSE_SYR2K: profiling_display_zsyr2k_info(); break; */

    case MORSE_LACPY: profiling_display_zlacpy_info(); break;
#if defined(PRECISION_z) || defined(PRECISION_c)
    case MORSE_PLGHE: profiling_display_zplghe_info(); break;
#endif
    case MORSE_PLGSY: profiling_display_zplgsy_info(); break;
    case MORSE_PLRNT: profiling_display_zplrnt_info(); break;
    /* case MORSE_LANGE: profiling_display_zlange_info(); break; */
    /* case MORSE_LANHE: profiling_display_zlanhe_info(); break; */
    /* case MORSE_LANSY: profiling_display_zlansy_info(); break; */

    case MORSE_POTRF: profiling_display_zpotrf_info(); break;

    case MORSE_GETRL: profiling_display_zgetrl_info(); break;
    case MORSE_GESSM: profiling_display_zgessm_info(); break;
    case MORSE_TSTRF: profiling_display_ztstrf_info(); break;
    case MORSE_SSSSM: profiling_display_zssssm_info(); break;

    case MORSE_GEQRT: profiling_display_zgeqrt_info(); break;
    case MORSE_UNMQR: profiling_display_zunmqr_info(); break;
    case MORSE_TSQRT: profiling_display_ztsqrt_info(); break;
    case MORSE_TSMQR: profiling_display_ztsmqr_info(); break;
    /* case MORSE_TTQRT: profiling_display_zttqrt_info(); break; */
    /* case MORSE_TTMQR: profiling_display_zttmqr_info(); break; */

    /* case MORSE_GELQT: profiling_display_zgelqt_info(); break; */
    /* case MORSE_UNMLQ: profiling_display_zunmlq_info(); break; */
    /* case MORSE_TSLQT: profiling_display_ztslqt_info(); break; */
    /* case MORSE_TSMLQ: profiling_display_ztsmlq_info(); break; */
    /* case MORSE_TTLQT: profiling_display_zttlqt_info(); break; */
    /* case MORSE_TTMLQ: profiling_display_zttmlq_info(); break; */
    default:
        return;
    }
}

