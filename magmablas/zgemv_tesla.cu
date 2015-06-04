/*
    -- MAGMA (version 1.3.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2012

       @precisions normal z -> z
       
*/
#include "common_magma.h"

#define magmablas_zgemv_tesla magmablas_zgemv

extern "C" void
magmablas_zgemv_tesla(char trans, magma_int_t m, magma_int_t n, 
                      cuDoubleComplex alpha, const cuDoubleComplex *A, magma_int_t lda, 
                                             const cuDoubleComplex *x, magma_int_t incx, 
                      cuDoubleComplex beta,  cuDoubleComplex       *y, magma_int_t incy) 
{
    cublasZgemv(trans, m, n, alpha, A, lda, x, incx, beta, y, incy);
}
