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
#include "morse_quark.h"

void morse_zdisplay_allprofile()
{
    magma_warning("morse_zdisplay_allprofile(quark)", "Profiling is not available with Quark");
}

void morse_zdisplay_oneprofile( morse_kernel_t kernel )
{
    magma_warning("morse_zdisplay_oneprofile(quark)", "Profiling is not available with Quark\n");
}

