/*
    -- MAGMA (version 1.6.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date November 2014

       @generated from commonblas_z.h normal z -> s, Sat Nov 15 19:53:59 2014
*/

#ifndef COMMONBLAS_S_H
#define COMMONBLAS_S_H

#ifdef __cplusplus
extern "C" {
#endif

/* ======================================================================
 * Internal prototypes
 */

// Tesla GEMM kernels
#define MAGMABLAS_SGEMM( name ) \
void magmablas_sgemm_##name( \
    float *C, const float *A, const float *B, \
    magma_int_t m, magma_int_t n, magma_int_t k, \
    magma_int_t lda, magma_int_t ldb, magma_int_t ldc, \
    float alpha, float beta )

MAGMABLAS_SGEMM( a_0  );
MAGMABLAS_SGEMM( ab_0 );
MAGMABLAS_SGEMM( N_N_64_16_16_16_4_special );
MAGMABLAS_SGEMM( N_N_64_16_16_16_4         );
MAGMABLAS_SGEMM( N_T_64_16_4_16_4          );
MAGMABLAS_SGEMM( T_N_32_32_8_8_8           );
MAGMABLAS_SGEMM( T_T_64_16_16_16_4_special );
MAGMABLAS_SGEMM( T_T_64_16_16_16_4         );
                   
void magmablas_sgemm_tesla(
    magma_trans_t transA, magma_trans_t transB, magma_int_t m, magma_int_t n, magma_int_t k,
    float alpha,
    const float *A, magma_int_t lda,
    const float *B, magma_int_t ldb,
    float beta,
    float *C, magma_int_t ldc );

void magmablas_sgemv_tesla(
    magma_trans_t trans, magma_int_t m, magma_int_t n,
    float alpha,
    const float *A, magma_int_t lda,
    const float *x, magma_int_t incx,
    float beta,
    float *y, magma_int_t incy );


// kernels used in snrm2, sgeqr2x-v4, laqps2_gpu, laqps3_gpu, slarfbx, slarfgx-v2, slarfx
__global__ void
magma_sgemv_kernel1(int m, const float * __restrict__ V, int ldv,
                    const float * __restrict__ c,
                    float *dwork);

__global__ void
magma_sgemv_kernel2(int m, int n, const float * __restrict__ V, int ldv,
                    const float * __restrict__ x, float *c);

__global__ void
magma_sgemv_kernel3(int m, const float * __restrict__ V, int ldv,
                    float *c, float *dwork,
                    float *tau);

__global__ void
magma_strmv_tkernel(float *T, int ldt, float *v,
                                    float *y);

__global__ void
magma_strmv_kernel2(const float *T, int ldt,
                    float *v, float *y, float *tau);

__global__ void
magma_snrm2_adjust_kernel(float *xnorm, float *c);

#ifdef __cplusplus
}
#endif

#endif /* COMMONBLAS_S_H */
