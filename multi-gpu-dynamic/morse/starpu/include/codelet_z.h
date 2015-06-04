/**
 *
 *  @file codelet_z.h
 *
 *  MAGMA codelets kernel
 *  MAGMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver,
 *  and INRIA Bordeaux Sud-Ouest
 *
 *  @version 2.3.1
 *  @author Cedric Augonnet
 *  @author Mathieu Faverge
 *  @date 2011-06-01
 *  @precisions normal z -> c d s
 *
 **/

#ifndef _CODELETS_Z_H_
#define _CODELETS_Z_H_

#define PRECISION_z

/*
 * Alphabetical order
 */
ZCODELETS_HEADER(gemm)
ZCODELETS_HEADER(geqrt)
ZCODELETS_HEADER(getrl)
ZCODELETS_HEADER(gessm)
ZCODELETS_HEADER(herk)
ZCODELETS_HEADER(lacpy)
ZCODELETS_HEADER(potrf)
ZCODELETS_HEADER(ssssm)
ZCODELETS_HEADER(trsm)
ZCODELETS_HEADER(tsmqr)
ZCODELETS_HEADER(tsqrt)
ZCODELETS_HEADER(tstrf)
ZCODELETS_HEADER(unmqr)

ZCODELETS_HEADER(trtri)
ZCODELETS_HEADER(lauum)
ZCODELETS_HEADER(trmm)

/*
 * CPU only functions
 */
ZCODELETS_CPU_HEADER(plrnt)
#if defined(PRECISION_z) || defined(PRECISION_c)
ZCODELETS_CPU_HEADER(plghe)
#endif
ZCODELETS_CPU_HEADER(plgsy)

#undef PRECISION_z

#endif /* _CODELETS_Z_H_ */
