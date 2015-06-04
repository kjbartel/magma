/**
 *
 * @file core_blas_dag.h
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.4.2
 * @author Mathieu Faverge
 * @date 2010-11-15
 *
 **/
#ifndef _PLASMA_CORE_BLAS_DAG_H_
#define _PLASMA_CORE_BLAS_DAG_H_

#if defined(QUARK_DOT_DAG_ENABLE) /* || 1 */
#define DAG_SET_PROPERTIES( _name, _color )                            \
  QUARK_Task_Flag_Set(task_flags, TASK_LABEL, (intptr_t)(_name));      \
  QUARK_Task_Flag_Set(task_flags, TASK_COLOR, (intptr_t)(_color));
#else
#define DAG_SET_PROPERTIES( _name, _color )
#endif

#define DAG_CORE_ASUM   DAG_SET_PROPERTIES( "ASUM"  , "white"   )
#define DAG_CORE_AXPY   DAG_SET_PROPERTIES( "AXPY"  , "white"   )
#define DAG_CORE_GELQT  DAG_SET_PROPERTIES( "GELQT" , "green"   )
#define DAG_CORE_GEMM   DAG_SET_PROPERTIES( "GEMM"  , "yellow"  )
#define DAG_CORE_GEQRT  DAG_SET_PROPERTIES( "GEQRT" , "green"   )
#define DAG_CORE_GESSM  DAG_SET_PROPERTIES( "GESSM" , "cyan"    )
#define DAG_CORE_GETRF  DAG_SET_PROPERTIES( "GETRF" , "green"   )
#define DAG_CORE_GETRIP DAG_SET_PROPERTIES( "GETRIP", "white"   )
#define DAG_CORE_GETRO  DAG_SET_PROPERTIES( "GETRO" , "white"   )
#define DAG_CORE_HEMM   DAG_SET_PROPERTIES( "HEMM"  , "white"   )
#define DAG_CORE_HER2K  DAG_SET_PROPERTIES( "HER2K" , "white"   )
#define DAG_CORE_HERK   DAG_SET_PROPERTIES( "HERK"  , "yellow"  )
#define DAG_CORE_LACPY  DAG_SET_PROPERTIES( "LACPY" , "white"   )
#define DAG_CORE_LAG2C  DAG_SET_PROPERTIES( "LAG2C" , "white"   )
#define DAG_CORE_LAG2Z  DAG_SET_PROPERTIES( "LAG2Z" , "white"   )
#define DAG_CORE_LANGE  DAG_SET_PROPERTIES( "LANGE" , "white"   )
#define DAG_CORE_LANHE  DAG_SET_PROPERTIES( "LANHE" , "white"   )
#define DAG_CORE_LANSY  DAG_SET_PROPERTIES( "LANSY" , "white"   )
#define DAG_CORE_LASET  DAG_SET_PROPERTIES( "LASET" , "orange"  )
#define DAG_CORE_LASWP  DAG_SET_PROPERTIES( "LASWP" , "orange"  )
#define DAG_CORE_LAUUM  DAG_SET_PROPERTIES( "LAUUM" , "white"   )
#define DAG_CORE_PLGHE  DAG_SET_PROPERTIES( "PLGHE" , "white"   )
#define DAG_CORE_PLGSY  DAG_SET_PROPERTIES( "PLGSY" , "white"   )
#define DAG_CORE_PLRNT  DAG_SET_PROPERTIES( "PLRNT" , "white"   )
#define DAG_CORE_POTRF  DAG_SET_PROPERTIES( "POTRF" , "green"   )
#define DAG_CORE_SHIFT  DAG_SET_PROPERTIES( "SHIFT" , "white"   )
#define DAG_CORE_SHIFTW DAG_SET_PROPERTIES( "SHIFTW", "white"   )
#define DAG_CORE_SSSSM  DAG_SET_PROPERTIES( "SSSSM" , "yellow"  )
#define DAG_CORE_SWPAB  DAG_SET_PROPERTIES( "SWPAB" , "white"   )
#define DAG_CORE_SYMM   DAG_SET_PROPERTIES( "SYMM"  , "white"   )
#define DAG_CORE_SYR2K  DAG_SET_PROPERTIES( "SYR2K" , "white"   )
#define DAG_CORE_SYRK   DAG_SET_PROPERTIES( "SYRK"  , "red"     )
#define DAG_CORE_TRMM   DAG_SET_PROPERTIES( "TRMM"  , "cyan"    )
#define DAG_CORE_TRSM   DAG_SET_PROPERTIES( "TRSM"  , "cyan"    )
#define DAG_CORE_TRTRI  DAG_SET_PROPERTIES( "TRTRI" , "white"   )
#define DAG_CORE_TSLQT  DAG_SET_PROPERTIES( "TSLQT" , "red"     )
#define DAG_CORE_TSMLQ  DAG_SET_PROPERTIES( "TSMLQ" , "yellow"  )
#define DAG_CORE_TSMQR  DAG_SET_PROPERTIES( "TSMQR" , "yellow"  )
#define DAG_CORE_TSQRT  DAG_SET_PROPERTIES( "TSQRT" , "red"     )
#define DAG_CORE_TSTRF  DAG_SET_PROPERTIES( "TSTRF" , "red"     )
#define DAG_CORE_TTLQT  DAG_SET_PROPERTIES( "TTLQT" , "pink"    )
#define DAG_CORE_TTMLQ  DAG_SET_PROPERTIES( "TTMLQ" , "magenta" )
#define DAG_CORE_TTMQR  DAG_SET_PROPERTIES( "TTMQR" , "magenta" )
#define DAG_CORE_TTQRT  DAG_SET_PROPERTIES( "TTQRT" , "pink"    )
#define DAG_CORE_UNMLQ  DAG_SET_PROPERTIES( "UNMLQ" , "cyan"    )
#define DAG_CORE_UNMQR  DAG_SET_PROPERTIES( "UNMQR" , "cyan"    )

#endif /* _PLASMA_CORE_BLAS_DAG_H_ */
