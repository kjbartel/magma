/**
 *
 * @file profiling.c
 *
 *  MAGMA auxiliary routines
 *  MAGMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 
 * @author Vijay Joshi
 * @date 2011-10-29 
 *
 **/
#include "morse_quark.h"

void morse_schedprofile_display(void)
{
    magma_warning("morse_schedprofile_display(quark)", "Scheduler profiling is not available with Quark\n");
}
