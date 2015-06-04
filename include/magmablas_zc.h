/*
 *   -- MAGMA (version 1.3.0) --
 *      Univ. of Tennessee, Knoxville
 *      Univ. of California, Berkeley
 *      Univ. of Colorado, Denver
 *      November 2012
 *
 * @precisions mixed zc -> ds
 */

#ifndef MAGMABLAS_ZC_H
#define MAGMABLAS_ZC_H

#ifdef __cplusplus
extern "C" {
#endif

  /* Mixed precision */
void magmablas_zcaxpycp(
    cuFloatComplex *R, cuDoubleComplex *X,
    magma_int_t m, cuDoubleComplex *B, cuDoubleComplex *W );

void magmablas_zaxpycp(
    cuDoubleComplex *R, cuDoubleComplex *X,
    magma_int_t m, cuDoubleComplex *B );

void magmablas_zclaswp(
    magma_int_t n,
    cuDoubleComplex *A, magma_int_t lda,
    cuFloatComplex *SA, magma_int_t m,
    const magma_int_t *ipiv, magma_int_t incx );

void magmablas_zlag2c(
    magma_int_t m, magma_int_t n,
    const cuDoubleComplex *A,  magma_int_t lda,
    cuFloatComplex        *SA, magma_int_t ldsa,
    magma_int_t *info );

void magmablas_clag2z(
    magma_int_t m, magma_int_t n, 
    const cuFloatComplex  *SA, magma_int_t ldsa, 
    cuDoubleComplex       *A,  magma_int_t lda, 
    magma_int_t *info );

void magmablas_zlat2c(
    char uplo, magma_int_t n, 
    const cuDoubleComplex *A,  magma_int_t lda, 
    cuFloatComplex        *SA, magma_int_t ldsa, 
    magma_int_t *info );

#ifdef __cplusplus
}
#endif

#endif // MAGMABLAS_ZC_H
