/**
 *
 *  @file zlocality.c
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

void morse_zlocality_allrestrict( uint32_t where )
{
#ifdef MORSE_USE_CUDA
    cl_zgemm_restrict_where( where );
    cl_ztrsm_restrict_where( where );
    /* cl_ztrmm_restrict_where( where ); */
    /* cl_zhemm_restrict_where( where ); */
    /* cl_zsymm_restrict_where( where ); */
    cl_zherk_restrict_where( where );
    /* cl_zsyrk_restrict_where( where ); */
    /* cl_zher2k_restrict_where( where ); */
    /* cl_zsyr2k_restrict_where( where ); */

    cl_zlacpy_restrict_where( where );
    /* cl_zplghe_restrict_where( where ); */
    /* cl_zplgsy_restrict_where( where ); */
    /* cl_zplrnt_restrict_where( where ); */
    /* cl_zlange_restrict_where( where ); */
    /* cl_zlanhe_restrict_where( where ); */
    /* cl_zlansy_restrict_where( where ); */

    cl_zpotrf_restrict_where( where );

    cl_zgetrl_restrict_where( where );
    cl_zgessm_restrict_where( where );
    cl_ztstrf_restrict_where( where );
    cl_zssssm_restrict_where( where );

    /* cl_zgeqrt_restrict_where( where ); */
    /* cl_zunmqr_restrict_where( where ); */
    /* cl_ztsqrt_restrict_where( where ); */
    /* cl_ztsmqr_restrict_where( where ); */
    /* cl_zttqrt_restrict_where( where ); */
    /* cl_zttmqr_restrict_where( where ); */

    /* cl_zgelqt_restrict_where( where ); */
    /* cl_zunmlq_restrict_where( where ); */
    /* cl_ztslqt_restrict_where( where ); */
    /* cl_ztsmlq_restrict_where( where ); */
    /* cl_zttlqt_restrict_where( where ); */
    /* cl_zttmlq_restrict_where( where ); */
#endif
}

void morse_zlocality_onerestrict( morse_kernel_t kernel, uint32_t where )
{
#ifdef MORSE_USE_CUDA
    switch( kernel ) {
    case MORSE_GEMM:  cl_zgemm_restrict_where( where );  break;
    case MORSE_TRSM:  cl_ztrsm_restrict_where( where );  break;
    /* case MORSE_TRMM:  cl_ztrmm_restrict_where( where );  break; */
    /* case MORSE_HEMM:  cl_zhemm_restrict_where( where );  break; */
    /* case MORSE_SYMM:  cl_zsymm_restrict_where( where );  break; */
    case MORSE_HERK:  cl_zherk_restrict_where( where );  break;
    /* case MORSE_SYRK:  cl_zsyrk_restrict_where( where );  break; */
    /* case MORSE_HER2K: cl_zher2k_restrict_where( where ); break; */
    /* case MORSE_SYR2K: cl_zsyr2k_restrict_where( where ); break; */

    case MORSE_LACPY: cl_zlacpy_restrict_where( where ); break;
    /* case MORSE_PLGHE: cl_zplghe_restrict_where( where ); break; */
    /* case MORSE_PLGSY: cl_zplgsy_restrict_where( where ); break; */
    /* case MORSE_PLRNT: cl_zplrnt_restrict_where( where ); break; */
    /* case MORSE_LANGE: cl_zlange_restrict_where( where ); break; */
    /* case MORSE_LANHE: cl_zlanhe_restrict_where( where ); break; */
    /* case MORSE_LANSY: cl_zlansy_restrict_where( where ); break; */

    case MORSE_POTRF: cl_zpotrf_restrict_where( where ); break;

    case MORSE_GETRL: cl_zgetrl_restrict_where( where ); break;
    case MORSE_GESSM: cl_zgessm_restrict_where( where ); break;
    case MORSE_TSTRF: cl_ztstrf_restrict_where( where ); break;
    case MORSE_SSSSM: cl_zssssm_restrict_where( where ); break;

    /* case MORSE_GEQRT: cl_zgeqrt_restrict_where( where ); break; */
    /* case MORSE_UNMQR: cl_zunmqr_restrict_where( where ); break; */
    /* case MORSE_TSQRT: cl_ztsqrt_restrict_where( where ); break; */
    /* case MORSE_TSMQR: cl_ztsmqr_restrict_where( where ); break; */
    /* case MORSE_TTQRT: cl_zttqrt_restrict_where( where ); break; */
    /* case MORSE_TTMQR: cl_zttmqr_restrict_where( where ); break; */

    /* case MORSE_GELQT: cl_zgelqt_restrict_where( where ); break; */
    /* case MORSE_UNMLQ: cl_zunmlq_restrict_where( where ); break; */
    /* case MORSE_TSLQT: cl_ztslqt_restrict_where( where ); break; */
    /* case MORSE_TSMLQ: cl_ztsmlq_restrict_where( where ); break; */
    /* case MORSE_TTLQT: cl_zttlqt_restrict_where( where ); break; */
    /* case MORSE_TTMLQ: cl_zttmlq_restrict_where( where ); break; */
    default:
      return;
    }
#endif
}

void morse_zlocality_allrestore( )
{
#ifdef MORSE_USE_CUDA
    cl_zgemm_restore_where();
    cl_ztrsm_restore_where();
    /* cl_ztrmm_restore_where(); */
    /* cl_zhemm_restore_where(); */
    /* cl_zsymm_restore_where(); */
    cl_zherk_restore_where();
    /* cl_zsyrk_restore_where(); */
    /* cl_zher2k_restore_where(); */
    /* cl_zsyr2k_restore_where(); */

    cl_zlacpy_restore_where();
    /* cl_zplghe_restore_where(); */
    /* cl_zplgsy_restore_where(); */
    /* cl_zplrnt_restore_where(); */
    /* cl_zlange_restore_where(); */
    /* cl_zlanhe_restore_where(); */
    /* cl_zlansy_restore_where(); */

    cl_zpotrf_restore_where();

    cl_zgetrl_restore_where();
    cl_zgessm_restore_where();
    cl_ztstrf_restore_where();
    cl_zssssm_restore_where();

    /* cl_zgeqrt_restore_where(); */
    /* cl_zunmqr_restore_where(); */
    /* cl_ztsqrt_restore_where(); */
    /* cl_ztsmqr_restore_where(); */
    /* cl_zttqrt_restore_where(); */
    /* cl_zttmqr_restore_where(); */

    /* cl_zgelqt_restore_where(); */
    /* cl_zunmlq_restore_where(); */
    /* cl_ztslqt_restore_where(); */
    /* cl_ztsmlq_restore_where(); */
    /* cl_zttlqt_restore_where(); */
    /* cl_zttmlq_restore_where(); */
#endif
}

void morse_zlocality_onerestore( morse_kernel_t kernel )
{
#ifdef MORSE_USE_CUDA
    switch( kernel ) {
    case MORSE_GEMM:  cl_zgemm_restore_where();  break;
    case MORSE_TRSM:  cl_ztrsm_restore_where();  break;
    /* case MORSE_TRMM:  cl_ztrmm_restore_where();  break; */
    /* case MORSE_HEMM:  cl_zhemm_restore_where();  break; */
    /* case MORSE_SYMM:  cl_zsymm_restore_where();  break; */
    case MORSE_HERK:  cl_zherk_restore_where();  break;
    /* case MORSE_SYRK:  cl_zsyrk_restore_where();  break; */
    /* case MORSE_HER2K: cl_zher2k_restore_where(); break; */
    /* case MORSE_SYR2K: cl_zsyr2k_restore_where(); break; */

    case MORSE_LACPY: cl_zlacpy_restore_where(); break;
    /* case MORSE_PLGHE: cl_zplghe_restore_where(); break; */
    /* case MORSE_PLGSY: cl_zplgsy_restore_where(); break; */
    /* case MORSE_PLRNT: cl_zplrnt_restore_where(); break; */
    /* case MORSE_LANGE: cl_zlange_restore_where(); break; */
    /* case MORSE_LANHE: cl_zlanhe_restore_where(); break; */
    /* case MORSE_LANSY: cl_zlansy_restore_where(); break; */

    case MORSE_POTRF: cl_zpotrf_restore_where(); break;

    case MORSE_GETRL: cl_zgetrl_restore_where(); break;
    case MORSE_GESSM: cl_zgessm_restore_where(); break;
    case MORSE_TSTRF: cl_ztstrf_restore_where(); break;
    case MORSE_SSSSM: cl_zssssm_restore_where(); break;

    /* case MORSE_GEQRT: cl_zgeqrt_restore_where(); break; */
    /* case MORSE_UNMQR: cl_zunmqr_restore_where(); break; */
    /* case MORSE_TSQRT: cl_ztsqrt_restore_where(); break; */
    /* case MORSE_TSMQR: cl_ztsmqr_restore_where(); break; */
    /* case MORSE_TTQRT: cl_zttqrt_restore_where(); break; */
    /* case MORSE_TTMQR: cl_zttmqr_restore_where(); break; */

    /* case MORSE_GELQT: cl_zgelqt_restore_where(); break; */
    /* case MORSE_UNMLQ: cl_zunmlq_restore_where(); break; */
    /* case MORSE_TSLQT: cl_ztslqt_restore_where(); break; */
    /* case MORSE_TSMLQ: cl_ztsmlq_restore_where(); break; */
    /* case MORSE_TTLQT: cl_zttlqt_restore_where(); break; */
    /* case MORSE_TTMLQ: cl_zttmlq_restore_where(); break; */
    default:
      return;
    }
#endif
}

