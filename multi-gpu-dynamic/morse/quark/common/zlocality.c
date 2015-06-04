/**
 *
 *  @file zlocality.c
 *
 *  MAGMA codelets kernel
 *  MAGMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver,
 *  and INRIA Bordeaux Sud-Ouest
 *
 *  @version
 *  @author Vijay Joshi
 *  @date 2011-10-29
 *  @precisions normal z -> s d c
 *
 **/
#include "morse_quark.h"

void morse_zlocality_allrestrict( uint32_t where )
{
    magma_warning("morse_zlocality_allrestrict(quark)", "Kernel locality cannot be specified with Quark");
}

void morse_zlocality_onerestrict( morse_kernel_t kernel, uint32_t where )
{
    magma_warning("morse_zlocality_onerestrict(quark)", "Kernel locality cannot be specified with Quark");
}

void morse_zlocality_allrestore( )
{
    magma_warning("morse_zlocality_allrestore(quark)", "Kernel locality cannot be specified with Quark");
}

void morse_zlocality_onerestore( morse_kernel_t kernel )
{
    magma_warning("morse_zlocality_onerestore(quark)", "Kernel locality cannot be specified with Quark");
}
