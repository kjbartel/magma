/**
 *
 *  @file morse_kernels.h
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
#ifndef _MORSE_KERNELS_H_
#define _MORSE_KERNELS_H_

typedef enum morse_kernel_e {
  MORSE_GEMM,
  MORSE_TRSM,
  MORSE_TRMM,
  MORSE_HEMM,
  MORSE_SYMM,
  MORSE_HERK,
  MORSE_SYRK,
  MORSE_HER2K,
  MORSE_SYR2K,

  MORSE_LACPY,
  MORSE_PLGHE,
  MORSE_PLGSY,
  MORSE_PLRNT,
  MORSE_LANGE,
  MORSE_LANHE,
  MORSE_LANSY,

  MORSE_POTRF,

  MORSE_GETRL,
  MORSE_GESSM,
  MORSE_TSTRF,
  MORSE_SSSSM,

  MORSE_GEQRT,
  MORSE_UNMQR,
  MORSE_ORMQR,
  MORSE_TSQRT,
  MORSE_TSMQR,
  MORSE_TTQRT,
  MORSE_TTMQR,

  MORSE_GELQT,
  MORSE_UNMLQ,
  MORSE_ORMLQ,
  MORSE_TSLQT,
  MORSE_TSMLQ,
  MORSE_TTLQT,
  MORSE_TTMLQ,

  MORSE_NBKERNELS
} morse_kernel_t;

#endif /* _MORSE_KERNELS_H_ */
