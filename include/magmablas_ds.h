/*
 *   -- MAGMA (version 1.3.0) --
 *      Univ. of Tennessee, Knoxville
 *      Univ. of California, Berkeley
 *      Univ. of Colorado, Denver
 *      November 2012
 *
 * @generated ds Wed Nov 14 22:52:26 2012
 */

#ifndef MAGMABLAS_DS_H
#define MAGMABLAS_DS_H

#ifdef __cplusplus
extern "C" {
#endif

  /* Mixed precision */
void magmablas_dsaxpycp(
    float *R, double *X,
    magma_int_t m, double *B, double *W );

void magmablas_daxpycp(
    double *R, double *X,
    magma_int_t m, double *B );

void magmablas_dslaswp(
    magma_int_t n,
    double *A, magma_int_t lda,
    float *SA, magma_int_t m,
    const magma_int_t *ipiv, magma_int_t incx );

void magmablas_dlag2s(
    magma_int_t m, magma_int_t n,
    const double *A,  magma_int_t lda,
    float        *SA, magma_int_t ldsa,
    magma_int_t *info );

void magmablas_slag2d(
    magma_int_t m, magma_int_t n, 
    const float  *SA, magma_int_t ldsa, 
    double       *A,  magma_int_t lda, 
    magma_int_t *info );

void magmablas_dlat2s(
    char uplo, magma_int_t n, 
    const double *A,  magma_int_t lda, 
    float        *SA, magma_int_t ldsa, 
    magma_int_t *info );

#ifdef __cplusplus
}
#endif

#endif // MAGMABLAS_DS_H
