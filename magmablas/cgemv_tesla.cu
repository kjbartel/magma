/*
    -- MAGMA (version 1.3.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2012

       @precisions normal d

*/
#include "common_magma.h"

#define magmablas_cgemv_tesla magmablas_cgemv

extern "C" void
magmablas_cgemv_tesla(char trans, magma_int_t m, magma_int_t n, 
                      cuFloatComplex alpha, const cuFloatComplex *A, magma_int_t lda, 
                                            const cuFloatComplex *x, magma_int_t incx, 
                      cuFloatComplex beta,  cuFloatComplex       *y, magma_int_t incy) 
{
    cublasCgemv(trans, m, n, alpha, A, lda, x, incx, beta, y, incy);
}
