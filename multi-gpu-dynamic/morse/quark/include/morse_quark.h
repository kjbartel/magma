/**
 *
 * @file morse_quark.h
 *
 *  MAGMA codelets kernel
 *  MAGMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver,
 *  and INRIA Bordeaux Sud-Ouest
 *
 * @version 2.3.1
 * @author Mathieu Faverge
 * @date 2011-06-01
 *
 **/

/***************************************************************************//**
 *  MAGMA facilities of interest to both src and magmablas directories
 **/
#ifndef _MORSE_QUARK_H_
#define _MORSE_QUARK_H_

#include <quark.h>
#include <plasma.h>
#include <core_blas.h>
#include <core_blas_dag.h>

#include "common.h"

/* 
 * Access to block pointer and leading dimension
 */
#define BLKADDR( desc, type, m, n ) ( (type*)morse_desc_getaddr( desc, m, n ) )

#endif /* _MORSE_QUARK_H_ */
