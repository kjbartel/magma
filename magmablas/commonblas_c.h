/*
    -- MAGMA (version 1.5.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date September 2014

       @generated from commonblas_z.h normal z -> c, Wed Sep 17 15:08:23 2014
*/

#ifndef COMMONBLAS_C_H
#define COMMONBLAS_C_H

#ifdef __cplusplus
extern "C" {
#endif

/* ======================================================================
 * Internal prototypes
 */

// Tesla GEMM kernels
#define MAGMABLAS_CGEMM( name ) \
void magmablas_cgemm_##name( \
    magmaFloatComplex *C, const magmaFloatComplex *A, const magmaFloatComplex *B, \
    magma_int_t m, magma_int_t n, magma_int_t k, \
    magma_int_t lda, magma_int_t ldb, magma_int_t ldc, \
    magmaFloatComplex alpha, magmaFloatComplex beta )

MAGMABLAS_CGEMM( a_0  );
MAGMABLAS_CGEMM( ab_0 );
MAGMABLAS_CGEMM( N_N_64_16_16_16_4_special );
MAGMABLAS_CGEMM( N_N_64_16_16_16_4         );
MAGMABLAS_CGEMM( N_T_64_16_4_16_4          );
MAGMABLAS_CGEMM( T_N_32_32_8_8_8           );
MAGMABLAS_CGEMM( T_T_64_16_16_16_4_special );
MAGMABLAS_CGEMM( T_T_64_16_16_16_4         );
                   
void magmablas_cgemm_tesla(
    magma_trans_t transA, magma_trans_t transB, magma_int_t m, magma_int_t n, magma_int_t k,
    magmaFloatComplex alpha,
    const magmaFloatComplex *A, magma_int_t lda,
    const magmaFloatComplex *B, magma_int_t ldb,
    magmaFloatComplex beta,
    magmaFloatComplex *C, magma_int_t ldc );

void magmablas_cgemv_tesla(
    magma_trans_t trans, magma_int_t m, magma_int_t n,
    magmaFloatComplex alpha,
    const magmaFloatComplex *A, magma_int_t lda,
    const magmaFloatComplex *x, magma_int_t incx,
    magmaFloatComplex beta,
    magmaFloatComplex *y, magma_int_t incy );


// kernels used in scnrm2, cgeqr2x-v4, laqps2_gpu, laqps3_gpu, clarfbx, clarfgx-v2, clarfx
__global__ void
magma_cgemv_kernel1(int m, const magmaFloatComplex * __restrict__ V, int ldv,
                    const magmaFloatComplex * __restrict__ c,
                    magmaFloatComplex *dwork);

__global__ void
magma_cgemv_kernel2(int m, int n, const magmaFloatComplex * __restrict__ V, int ldv,
                    const magmaFloatComplex * __restrict__ x, magmaFloatComplex *c);

__global__ void
magma_cgemv_kernel3(int m, const magmaFloatComplex * __restrict__ V, int ldv,
                    magmaFloatComplex *c, magmaFloatComplex *dwork,
                    magmaFloatComplex *tau);

__global__ void
magma_ctrmv_tkernel(magmaFloatComplex *T, int ldt, magmaFloatComplex *v,
                                    magmaFloatComplex *y);

__global__ void
magma_ctrmv_kernel2(const magmaFloatComplex *T, int ldt,
                    magmaFloatComplex *v, magmaFloatComplex *y, magmaFloatComplex *tau);

__global__ void
magma_scnrm2_adjust_kernel(float *xnorm, magmaFloatComplex *c);

#ifdef __cplusplus
}
#endif

#endif /* COMMONBLAS_C_H */
