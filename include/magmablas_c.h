/*
 *   -- MAGMA (version 1.2.1) --
 *      Univ. of Tennessee, Knoxville
 *      Univ. of California, Berkeley
 *      Univ. of Colorado, Denver
 *      June 2012
 *
 * @generated c Thu Jun 28 12:30:02 2012
 */

#ifndef _MAGMABLAS_C_H_
#define _MAGMABLAS_C_H_

#define PRECISION_c

#ifdef __cplusplus
extern "C" {
#endif

  /*
   * Interface to clean
   */
float cpu_gpu_cdiff(             magma_int_t M, magma_int_t N, 
                  cuFloatComplex * a, magma_int_t lda, 
                  cuFloatComplex *da, magma_int_t ldda);

void czero_32x32_block(           cuFloatComplex *, magma_int_t);

void czero_nbxnb_block(           magma_int_t, cuFloatComplex *, magma_int_t);

void magmablas_cpermute_long(     cuFloatComplex *, magma_int_t, 
                  magma_int_t *, magma_int_t, magma_int_t);

void magmablas_cpermute_long2(magma_int_t n, cuFloatComplex *, magma_int_t, 
                  magma_int_t *, magma_int_t, magma_int_t);

void magmablas_cpermute_long3( cuFloatComplex *dAT, magma_int_t lda, 
                               magma_int_t *ipiv, magma_int_t nb, magma_int_t ind );

  /*
   * Transpose functions
   */
void magmablas_cinplace_transpose(cuFloatComplex *, magma_int_t, magma_int_t);

void magmablas_ctranspose(        cuFloatComplex *, magma_int_t, 
                  cuFloatComplex *, magma_int_t, 
                  magma_int_t, magma_int_t);

void magmablas_ctranspose2(       cuFloatComplex *, magma_int_t, 
                  cuFloatComplex *, magma_int_t, 
                  magma_int_t, magma_int_t);

void magmablas_ctranspose2s(cuFloatComplex *odata, magma_int_t ldo,
                       cuFloatComplex *idata, magma_int_t ldi,
                       magma_int_t m, magma_int_t n, cudaStream_t *stream );

void magmablas_cgetmatrix_transpose(  magma_int_t m, magma_int_t n,
                                      cuFloatComplex *dat, magma_int_t ldda,
                                      cuFloatComplex  *ha, magma_int_t lda,
                                      cuFloatComplex  *dB, magma_int_t lddb, magma_int_t nb );
void magmablas_csetmatrix_transpose(  magma_int_t m, magma_int_t n,
                                      cuFloatComplex  *ha, magma_int_t lda, 
                                      cuFloatComplex *dat, magma_int_t ldda,
                                      cuFloatComplex  *dB, magma_int_t lddb, magma_int_t nb );

  /*
   * Multi-GPU functions
   */
void magmablas_cgetmatrix_transpose_mgpu(
                  magma_int_t num_gpus, cudaStream_t **stream0,
                  cuFloatComplex **dat, magma_int_t ldda,
                  cuFloatComplex   *ha, magma_int_t lda,
                  cuFloatComplex  **dB, magma_int_t lddb,
                  magma_int_t m, magma_int_t n, magma_int_t nb);

void magmablas_csetmatrix_transpose_mgpu(
                  magma_int_t num_gpus, cudaStream_t **stream0,
                  cuFloatComplex  *ha,  magma_int_t lda, 
                  cuFloatComplex **dat, magma_int_t ldda, magma_int_t starti,
                  cuFloatComplex **dB,  magma_int_t lddb,
                  magma_int_t m, magma_int_t n, magma_int_t nb);

void magmablas_cgetmatrix_1D_bcyclic( magma_int_t m, magma_int_t n,
                                      cuFloatComplex  *da[], magma_int_t ldda,
                                      cuFloatComplex  *ha, magma_int_t lda,
                                      magma_int_t num_gpus, magma_int_t nb );

void magmablas_csetmatrix_1D_bcyclic( magma_int_t m, magma_int_t n,
                                      cuFloatComplex  *ha, magma_int_t lda,
                                      cuFloatComplex  *da[], magma_int_t ldda,
                                      magma_int_t num_gpus, magma_int_t nb );

void magmablas_chemm_mgpu(
    char side, char uplo, magma_int_t m, magma_int_t n,
    cuFloatComplex alpha, cuFloatComplex *dA[], magma_int_t ldda,  magma_int_t offset,
                           cuFloatComplex *dB[], magma_int_t lddb,
    cuFloatComplex beta,  cuFloatComplex *dC[], magma_int_t lddc,
                           cuFloatComplex *C,    magma_int_t ldc,
    magma_int_t ngpu, magma_int_t nb, cudaStream_t streams[][20], magma_int_t nstream );

void magmablas_cher2k_mgpu2(
    char uplo, char trans, magma_int_t n, magma_int_t k,
    cuFloatComplex alpha, cuFloatComplex *dA[], magma_int_t lda,
                           cuFloatComplex *dB[], magma_int_t ldb,
    float beta,           cuFloatComplex *dC[], magma_int_t ldc,  magma_int_t offset,
    magma_int_t ngpu, magma_int_t nb, cudaStream_t streams[][10], magma_int_t nstream );

  /*
   * LAPACK auxiliary functions
   */
void   magmablas_clacpy( char uplo, 
             magma_int_t m, magma_int_t n, 
             cuFloatComplex *A, magma_int_t lda, 
             cuFloatComplex *B, magma_int_t ldb);

float magmablas_clange( char norm, 
             magma_int_t m, magma_int_t n, 
             cuFloatComplex *A, magma_int_t lda, float *WORK);

float magmablas_clanhe( char norm, char uplo, 
             magma_int_t n,
             cuFloatComplex *A, magma_int_t lda, float *WORK);

float magmablas_clansy( char norm, char uplo,
             magma_int_t n, 
             cuFloatComplex *A, magma_int_t lda, float *WORK);

void   magmablas_clascl( char type, magma_int_t kl, magma_int_t ku,
             float cfrom, float cto,
             magma_int_t m, magma_int_t n,
             cuFloatComplex *A, magma_int_t lda, magma_int_t *info );

void   magmablas_claset( char uplo, magma_int_t m, magma_int_t n,
             cuFloatComplex *A, magma_int_t lda);

void   magmablas_claset_identity(
             magma_int_t m, magma_int_t n,
             cuFloatComplex *A, magma_int_t lda);

void   magmablas_claswp( magma_int_t N, 
             cuFloatComplex *dAT, magma_int_t lda, 
             magma_int_t i1,  magma_int_t i2, 
             magma_int_t *ipiv, magma_int_t inci );

void   magmablas_claswpx(magma_int_t N, 
             cuFloatComplex *dAT, magma_int_t ldx, magma_int_t ldy, 
             magma_int_t i1, magma_int_t i2,
             magma_int_t *ipiv, magma_int_t inci );

  /*
   * Level 1 BLAS
   */
void   magmablas_cswap(   magma_int_t N, 
              cuFloatComplex *dA1, magma_int_t lda1, 
              cuFloatComplex *dA2, magma_int_t lda2 );

void   magmablas_cswapblk(char storev, 
              magma_int_t N, 
              cuFloatComplex *dA1, magma_int_t lda1, 
              cuFloatComplex *dA2, magma_int_t lda2,
              magma_int_t i1, magma_int_t i2, 
              magma_int_t *ipiv, magma_int_t inci, 
              magma_int_t offset);

void magmablas_cswapdblk(magma_int_t n, magma_int_t nb,
             cuFloatComplex *dA1, magma_int_t ldda1, magma_int_t inca1,
             cuFloatComplex *dA2, magma_int_t ldda2, magma_int_t inca2 );

  /*
   * Level 2 BLAS
   */
void magmablas_cgemv(char t, magma_int_t M, magma_int_t N, 
             cuFloatComplex alpha,
             cuFloatComplex *A, magma_int_t lda, 
             cuFloatComplex * X, magma_int_t incX, 
             cuFloatComplex beta, 
             cuFloatComplex *Y, magma_int_t incY);

#if defined(PRECISION_z) || defined(PRECISION_c)
magma_int_t magmablas_chemv(char u, magma_int_t N, 
                            cuFloatComplex alpha, 
                            cuFloatComplex *A, magma_int_t lda, 
                            cuFloatComplex *X, magma_int_t incX, 
                            cuFloatComplex beta, 
                            cuFloatComplex *Y, magma_int_t incY);
#endif

magma_int_t magmablas_csymv(char u, magma_int_t N, 
                            cuFloatComplex alpha, 
                            cuFloatComplex *A, magma_int_t lda, 
                            cuFloatComplex *X, magma_int_t incX, 
                            cuFloatComplex beta, 
                            cuFloatComplex *Y, magma_int_t incY);

  /*
   * Level 3 BLAS
   */
void magmablas_cgemm(char tA, char tB,
             magma_int_t m, magma_int_t n, magma_int_t k, 
             cuFloatComplex alpha,
             const cuFloatComplex *A, magma_int_t lda, 
             const cuFloatComplex *B, magma_int_t ldb, 
             cuFloatComplex beta,
             cuFloatComplex *C, magma_int_t ldc);

void magmablas_cgemm_fermi80(char tA, char tB, 
                 magma_int_t m, magma_int_t n, magma_int_t k,
                 cuFloatComplex alpha, 
                 const cuFloatComplex *A, magma_int_t lda, 
                 const cuFloatComplex *B, magma_int_t ldb,
                 cuFloatComplex beta, 
                 cuFloatComplex *C, magma_int_t ldc);

void magmablas_cgemm_fermi64(char tA, char tB, 
                 magma_int_t m, magma_int_t n, magma_int_t k,
                 cuFloatComplex alpha, 
                 const cuFloatComplex *A, magma_int_t lda, 
                 const cuFloatComplex *B, magma_int_t ldb, 
                 cuFloatComplex beta, 
                 cuFloatComplex *C, magma_int_t ldc);

void magmablas_chemm(char s, char u,          
             magma_int_t m, magma_int_t n,
             cuFloatComplex alpha, 
             const cuFloatComplex *A, magma_int_t lda,
             const cuFloatComplex *B, magma_int_t ldb,
             cuFloatComplex beta, 
             cuFloatComplex *C, magma_int_t ldc);

void magmablas_csymm(char s, char u,
             magma_int_t m, magma_int_t n,
             cuFloatComplex alpha, 
             const cuFloatComplex *A, magma_int_t lda, 
             const cuFloatComplex *B, magma_int_t ldb,
             cuFloatComplex beta,
             cuFloatComplex *C, magma_int_t ldc);

void magmablas_csyrk(char u, char t,
             magma_int_t n, magma_int_t k, 
             cuFloatComplex alpha, 
             const cuFloatComplex *A, magma_int_t lda,
             cuFloatComplex beta,
             cuFloatComplex *C, magma_int_t ldc);

void magmablas_cherk(char u, char t,
             magma_int_t n, magma_int_t k, 
             float  alpha, 
             const cuFloatComplex *A, magma_int_t lda,
             float  beta, 
             cuFloatComplex *C, magma_int_t ldc);

void magmablas_csyr2k(char u, char t,
              magma_int_t n, magma_int_t k,
              cuFloatComplex alpha, 
              const cuFloatComplex *A, magma_int_t lda,
              const cuFloatComplex *B, magma_int_t ldb, 
              cuFloatComplex beta, 
              cuFloatComplex *C, magma_int_t ldc);

void magmablas_cher2k(char u, char t,
              magma_int_t n, magma_int_t k, 
              cuFloatComplex alpha, 
              const cuFloatComplex *A, magma_int_t lda, 
              const cuFloatComplex *B, magma_int_t ldb,
              float  beta,
              cuFloatComplex *C, magma_int_t ldc);

void magmablas_ctrmm(char s, char u, char t,  char d, 
             magma_int_t m, magma_int_t n,
             cuFloatComplex alpha,
             const cuFloatComplex *A, magma_int_t lda,
             cuFloatComplex *B, magma_int_t ldb);

void magmablas_ctrsm(char s, char u, char t, char d,
             magma_int_t m, magma_int_t n,
             cuFloatComplex alpha,
             /*const*/ cuFloatComplex *A, magma_int_t lda,
             cuFloatComplex *B, magma_int_t ldb);


  /*
   * Wrappers for platform independence.
   * These wrap CUBLAS or AMD OpenCL BLAS functions.
   */

// ========================================
// copying vectors
// set copies host to device
// get copies device to host

void magma_csetvector(
    magma_int_t n,
    cuFloatComplex const *hx_src, magma_int_t incx,
    cuFloatComplex       *dy_dst, magma_int_t incy );

void magma_cgetvector(
    magma_int_t n,
    cuFloatComplex const *dx_src, magma_int_t incx,
    cuFloatComplex       *hy_dst, magma_int_t incy );

void magma_csetvector_async(
    magma_int_t n,
    cuFloatComplex const *hx_src, magma_int_t incx,
    cuFloatComplex       *dy_dst, magma_int_t incy,
    magma_stream_t stream );

void magma_cgetvector_async(
    magma_int_t n,
    cuFloatComplex const *dx_src, magma_int_t incx,
    cuFloatComplex       *hy_dst, magma_int_t incy,
    magma_stream_t stream );


// ========================================
// copying sub-matrices (contiguous columns)
// set copies host to device
// get copies device to host
// cpy copies device to device (with CUDA unified addressing, can be same or different devices)

void magma_csetmatrix(
    magma_int_t m, magma_int_t n,
    cuFloatComplex const *hA_src, magma_int_t lda,
    cuFloatComplex       *dB_dst, magma_int_t ldb );

void magma_cgetmatrix(
    magma_int_t m, magma_int_t n,
    cuFloatComplex const *dA_src, magma_int_t lda,
    cuFloatComplex       *hB_dst, magma_int_t ldb );

void magma_csetmatrix_async(
    magma_int_t m, magma_int_t n,
    cuFloatComplex const *hA_src, magma_int_t lda,
    cuFloatComplex       *dB_dst, magma_int_t ldb,
    magma_stream_t stream );

void magma_cgetmatrix_async(
    magma_int_t m, magma_int_t n,
    cuFloatComplex const *dA_src, magma_int_t lda,
    cuFloatComplex       *hB_dst, magma_int_t ldb,
    magma_stream_t stream );

void magma_ccopymatrix(
    magma_int_t m, magma_int_t n,
    cuFloatComplex const *dA_src, magma_int_t lda,
    cuFloatComplex       *dB_dst, magma_int_t ldb );

void magma_ccopymatrix_async(
    magma_int_t m, magma_int_t n,
    cuFloatComplex const *dA_src, magma_int_t lda,
    cuFloatComplex       *dB_dst, magma_int_t ldb,
    magma_stream_t stream );


// ========================================
// Level 1 BLAS

void magma_cswap(
    magma_int_t n,
    cuFloatComplex *dx, magma_int_t incx,
    cuFloatComplex *dy, magma_int_t incy );

magma_int_t magma_icamax(
    magma_int_t n,
    cuFloatComplex *dx, magma_int_t incx );

// ========================================
// Level 2 BLAS

void magma_cgemv(
    magma_trans_t transA,
    magma_int_t m, magma_int_t n,
    cuFloatComplex alpha, cuFloatComplex const *dA, magma_int_t lda,
                           cuFloatComplex const *dx, magma_int_t incx,
    cuFloatComplex beta,  cuFloatComplex       *dy, magma_int_t incy );

void magma_chemv(
    magma_uplo_t uplo,
    magma_int_t n,
    cuFloatComplex alpha, cuFloatComplex const *dA, magma_int_t lda,
                           cuFloatComplex const *dx, magma_int_t incx,
    cuFloatComplex beta,  cuFloatComplex       *dy, magma_int_t incy );

void magma_ctrsv(
    magma_uplo_t uplo, magma_trans_t trans, magma_diag_t diag, 
    magma_int_t n, 
    cuFloatComplex const *dA, magma_int_t lda, 
    cuFloatComplex       *dx, magma_int_t incx );

// ========================================
// Level 3 BLAS

void magma_cgemm(
    magma_trans_t transA, magma_trans_t transB,
    magma_int_t m, magma_int_t n, magma_int_t k,
    cuFloatComplex alpha, cuFloatComplex const *dA, magma_int_t lda,
                           cuFloatComplex const *dB, magma_int_t ldb,
    cuFloatComplex beta,  cuFloatComplex       *dC, magma_int_t ldc );

void magma_chemm(
    magma_side_t side, magma_uplo_t uplo,
    magma_int_t m, magma_int_t n,
    cuFloatComplex alpha, cuFloatComplex const *dA, magma_int_t lda,
                           cuFloatComplex const *dB, magma_int_t ldb,
    cuFloatComplex beta,  cuFloatComplex       *dC, magma_int_t ldc );

void magma_cherk(
    magma_uplo_t uplo, magma_trans_t trans,
    magma_int_t n, magma_int_t k,
    float alpha, cuFloatComplex const *dA, magma_int_t lda,
    float beta,  cuFloatComplex       *dC, magma_int_t ldc );

void magma_cher2k(
    magma_uplo_t uplo, magma_trans_t trans,
    magma_int_t n, magma_int_t k,
    cuFloatComplex alpha, cuFloatComplex const *dA, magma_int_t lda,
                           cuFloatComplex const *dB, magma_int_t ldb,
    float beta,           cuFloatComplex       *dC, magma_int_t ldc );

void magma_ctrmm(
    magma_side_t side, magma_uplo_t uplo, magma_trans_t trans, magma_diag_t diag,
    magma_int_t m, magma_int_t n,
    cuFloatComplex alpha, cuFloatComplex const *dA, magma_int_t lda,
                           cuFloatComplex       *dB, magma_int_t ldb );

void magma_ctrsm(
    magma_side_t side, magma_uplo_t uplo, magma_trans_t trans, magma_diag_t diag,
    magma_int_t m, magma_int_t n,
    cuFloatComplex alpha, cuFloatComplex const *dA, magma_int_t lda,
                           cuFloatComplex       *dB, magma_int_t ldb );

#ifdef __cplusplus
}
#endif

#undef PRECISION_c
#endif
