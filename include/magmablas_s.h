/*
 *   -- MAGMA (version 1.2.1) --
 *      Univ. of Tennessee, Knoxville
 *      Univ. of California, Berkeley
 *      Univ. of Colorado, Denver
 *      June 2012
 *
 * @generated s Thu Jun 28 12:30:02 2012
 */

#ifndef _MAGMABLAS_S_H_
#define _MAGMABLAS_S_H_

#define PRECISION_s

#ifdef __cplusplus
extern "C" {
#endif

  /*
   * Interface to clean
   */
float cpu_gpu_sdiff(             magma_int_t M, magma_int_t N, 
                  float * a, magma_int_t lda, 
                  float *da, magma_int_t ldda);

void szero_32x32_block(           float *, magma_int_t);

void szero_nbxnb_block(           magma_int_t, float *, magma_int_t);

void magmablas_spermute_long(     float *, magma_int_t, 
                  magma_int_t *, magma_int_t, magma_int_t);

void magmablas_spermute_long2(magma_int_t n, float *, magma_int_t, 
                  magma_int_t *, magma_int_t, magma_int_t);

void magmablas_spermute_long3( float *dAT, magma_int_t lda, 
                               magma_int_t *ipiv, magma_int_t nb, magma_int_t ind );

  /*
   * Transpose functions
   */
void magmablas_sinplace_transpose(float *, magma_int_t, magma_int_t);

void magmablas_stranspose(        float *, magma_int_t, 
                  float *, magma_int_t, 
                  magma_int_t, magma_int_t);

void magmablas_stranspose2(       float *, magma_int_t, 
                  float *, magma_int_t, 
                  magma_int_t, magma_int_t);

void magmablas_stranspose2s(float *odata, magma_int_t ldo,
                       float *idata, magma_int_t ldi,
                       magma_int_t m, magma_int_t n, cudaStream_t *stream );

void magmablas_sgetmatrix_transpose(  magma_int_t m, magma_int_t n,
                                      float *dat, magma_int_t ldda,
                                      float  *ha, magma_int_t lda,
                                      float  *dB, magma_int_t lddb, magma_int_t nb );
void magmablas_ssetmatrix_transpose(  magma_int_t m, magma_int_t n,
                                      float  *ha, magma_int_t lda, 
                                      float *dat, magma_int_t ldda,
                                      float  *dB, magma_int_t lddb, magma_int_t nb );

  /*
   * Multi-GPU functions
   */
void magmablas_sgetmatrix_transpose_mgpu(
                  magma_int_t num_gpus, cudaStream_t **stream0,
                  float **dat, magma_int_t ldda,
                  float   *ha, magma_int_t lda,
                  float  **dB, magma_int_t lddb,
                  magma_int_t m, magma_int_t n, magma_int_t nb);

void magmablas_ssetmatrix_transpose_mgpu(
                  magma_int_t num_gpus, cudaStream_t **stream0,
                  float  *ha,  magma_int_t lda, 
                  float **dat, magma_int_t ldda, magma_int_t starti,
                  float **dB,  magma_int_t lddb,
                  magma_int_t m, magma_int_t n, magma_int_t nb);

void magmablas_sgetmatrix_1D_bcyclic( magma_int_t m, magma_int_t n,
                                      float  *da[], magma_int_t ldda,
                                      float  *ha, magma_int_t lda,
                                      magma_int_t num_gpus, magma_int_t nb );

void magmablas_ssetmatrix_1D_bcyclic( magma_int_t m, magma_int_t n,
                                      float  *ha, magma_int_t lda,
                                      float  *da[], magma_int_t ldda,
                                      magma_int_t num_gpus, magma_int_t nb );

void magmablas_ssymm_mgpu(
    char side, char uplo, magma_int_t m, magma_int_t n,
    float alpha, float *dA[], magma_int_t ldda,  magma_int_t offset,
                           float *dB[], magma_int_t lddb,
    float beta,  float *dC[], magma_int_t lddc,
                           float *C,    magma_int_t ldc,
    magma_int_t ngpu, magma_int_t nb, cudaStream_t streams[][20], magma_int_t nstream );

void magmablas_ssyr2k_mgpu2(
    char uplo, char trans, magma_int_t n, magma_int_t k,
    float alpha, float *dA[], magma_int_t lda,
                           float *dB[], magma_int_t ldb,
    float beta,           float *dC[], magma_int_t ldc,  magma_int_t offset,
    magma_int_t ngpu, magma_int_t nb, cudaStream_t streams[][10], magma_int_t nstream );

  /*
   * LAPACK auxiliary functions
   */
void   magmablas_slacpy( char uplo, 
             magma_int_t m, magma_int_t n, 
             float *A, magma_int_t lda, 
             float *B, magma_int_t ldb);

float magmablas_slange( char norm, 
             magma_int_t m, magma_int_t n, 
             float *A, magma_int_t lda, float *WORK);

float magmablas_slansy( char norm, char uplo, 
             magma_int_t n,
             float *A, magma_int_t lda, float *WORK);

float magmablas_slansy( char norm, char uplo,
             magma_int_t n, 
             float *A, magma_int_t lda, float *WORK);

void   magmablas_slascl( char type, magma_int_t kl, magma_int_t ku,
             float cfrom, float cto,
             magma_int_t m, magma_int_t n,
             float *A, magma_int_t lda, magma_int_t *info );

void   magmablas_slaset( char uplo, magma_int_t m, magma_int_t n,
             float *A, magma_int_t lda);

void   magmablas_slaset_identity(
             magma_int_t m, magma_int_t n,
             float *A, magma_int_t lda);

void   magmablas_slaswp( magma_int_t N, 
             float *dAT, magma_int_t lda, 
             magma_int_t i1,  magma_int_t i2, 
             magma_int_t *ipiv, magma_int_t inci );

void   magmablas_slaswpx(magma_int_t N, 
             float *dAT, magma_int_t ldx, magma_int_t ldy, 
             magma_int_t i1, magma_int_t i2,
             magma_int_t *ipiv, magma_int_t inci );

  /*
   * Level 1 BLAS
   */
void   magmablas_sswap(   magma_int_t N, 
              float *dA1, magma_int_t lda1, 
              float *dA2, magma_int_t lda2 );

void   magmablas_sswapblk(char storev, 
              magma_int_t N, 
              float *dA1, magma_int_t lda1, 
              float *dA2, magma_int_t lda2,
              magma_int_t i1, magma_int_t i2, 
              magma_int_t *ipiv, magma_int_t inci, 
              magma_int_t offset);

void magmablas_sswapdblk(magma_int_t n, magma_int_t nb,
             float *dA1, magma_int_t ldda1, magma_int_t inca1,
             float *dA2, magma_int_t ldda2, magma_int_t inca2 );

  /*
   * Level 2 BLAS
   */
void magmablas_sgemv(char t, magma_int_t M, magma_int_t N, 
             float alpha,
             float *A, magma_int_t lda, 
             float * X, magma_int_t incX, 
             float beta, 
             float *Y, magma_int_t incY);

#if defined(PRECISION_z) || defined(PRECISION_c)
magma_int_t magmablas_ssymv(char u, magma_int_t N, 
                            float alpha, 
                            float *A, magma_int_t lda, 
                            float *X, magma_int_t incX, 
                            float beta, 
                            float *Y, magma_int_t incY);
#endif

magma_int_t magmablas_ssymv(char u, magma_int_t N, 
                            float alpha, 
                            float *A, magma_int_t lda, 
                            float *X, magma_int_t incX, 
                            float beta, 
                            float *Y, magma_int_t incY);

  /*
   * Level 3 BLAS
   */
void magmablas_sgemm(char tA, char tB,
             magma_int_t m, magma_int_t n, magma_int_t k, 
             float alpha,
             const float *A, magma_int_t lda, 
             const float *B, magma_int_t ldb, 
             float beta,
             float *C, magma_int_t ldc);

void magmablas_sgemm_fermi80(char tA, char tB, 
                 magma_int_t m, magma_int_t n, magma_int_t k,
                 float alpha, 
                 const float *A, magma_int_t lda, 
                 const float *B, magma_int_t ldb,
                 float beta, 
                 float *C, magma_int_t ldc);

void magmablas_sgemm_fermi64(char tA, char tB, 
                 magma_int_t m, magma_int_t n, magma_int_t k,
                 float alpha, 
                 const float *A, magma_int_t lda, 
                 const float *B, magma_int_t ldb, 
                 float beta, 
                 float *C, magma_int_t ldc);

void magmablas_ssymm(char s, char u,          
             magma_int_t m, magma_int_t n,
             float alpha, 
             const float *A, magma_int_t lda,
             const float *B, magma_int_t ldb,
             float beta, 
             float *C, magma_int_t ldc);

void magmablas_ssymm(char s, char u,
             magma_int_t m, magma_int_t n,
             float alpha, 
             const float *A, magma_int_t lda, 
             const float *B, magma_int_t ldb,
             float beta,
             float *C, magma_int_t ldc);

void magmablas_ssyrk(char u, char t,
             magma_int_t n, magma_int_t k, 
             float alpha, 
             const float *A, magma_int_t lda,
             float beta,
             float *C, magma_int_t ldc);

void magmablas_ssyrk(char u, char t,
             magma_int_t n, magma_int_t k, 
             float  alpha, 
             const float *A, magma_int_t lda,
             float  beta, 
             float *C, magma_int_t ldc);

void magmablas_ssyr2k(char u, char t,
              magma_int_t n, magma_int_t k,
              float alpha, 
              const float *A, magma_int_t lda,
              const float *B, magma_int_t ldb, 
              float beta, 
              float *C, magma_int_t ldc);

void magmablas_ssyr2k(char u, char t,
              magma_int_t n, magma_int_t k, 
              float alpha, 
              const float *A, magma_int_t lda, 
              const float *B, magma_int_t ldb,
              float  beta,
              float *C, magma_int_t ldc);

void magmablas_strmm(char s, char u, char t,  char d, 
             magma_int_t m, magma_int_t n,
             float alpha,
             const float *A, magma_int_t lda,
             float *B, magma_int_t ldb);

void magmablas_strsm(char s, char u, char t, char d,
             magma_int_t m, magma_int_t n,
             float alpha,
             /*const*/ float *A, magma_int_t lda,
             float *B, magma_int_t ldb);


  /*
   * Wrappers for platform independence.
   * These wrap CUBLAS or AMD OpenCL BLAS functions.
   */

// ========================================
// copying vectors
// set copies host to device
// get copies device to host

void magma_ssetvector(
    magma_int_t n,
    float const *hx_src, magma_int_t incx,
    float       *dy_dst, magma_int_t incy );

void magma_sgetvector(
    magma_int_t n,
    float const *dx_src, magma_int_t incx,
    float       *hy_dst, magma_int_t incy );

void magma_ssetvector_async(
    magma_int_t n,
    float const *hx_src, magma_int_t incx,
    float       *dy_dst, magma_int_t incy,
    magma_stream_t stream );

void magma_sgetvector_async(
    magma_int_t n,
    float const *dx_src, magma_int_t incx,
    float       *hy_dst, magma_int_t incy,
    magma_stream_t stream );


// ========================================
// copying sub-matrices (contiguous columns)
// set copies host to device
// get copies device to host
// cpy copies device to device (with CUDA unified addressing, can be same or different devices)

void magma_ssetmatrix(
    magma_int_t m, magma_int_t n,
    float const *hA_src, magma_int_t lda,
    float       *dB_dst, magma_int_t ldb );

void magma_sgetmatrix(
    magma_int_t m, magma_int_t n,
    float const *dA_src, magma_int_t lda,
    float       *hB_dst, magma_int_t ldb );

void magma_ssetmatrix_async(
    magma_int_t m, magma_int_t n,
    float const *hA_src, magma_int_t lda,
    float       *dB_dst, magma_int_t ldb,
    magma_stream_t stream );

void magma_sgetmatrix_async(
    magma_int_t m, magma_int_t n,
    float const *dA_src, magma_int_t lda,
    float       *hB_dst, magma_int_t ldb,
    magma_stream_t stream );

void magma_scopymatrix(
    magma_int_t m, magma_int_t n,
    float const *dA_src, magma_int_t lda,
    float       *dB_dst, magma_int_t ldb );

void magma_scopymatrix_async(
    magma_int_t m, magma_int_t n,
    float const *dA_src, magma_int_t lda,
    float       *dB_dst, magma_int_t ldb,
    magma_stream_t stream );


// ========================================
// Level 1 BLAS

void magma_sswap(
    magma_int_t n,
    float *dx, magma_int_t incx,
    float *dy, magma_int_t incy );

magma_int_t magma_isamax(
    magma_int_t n,
    float *dx, magma_int_t incx );

// ========================================
// Level 2 BLAS

void magma_sgemv(
    magma_trans_t transA,
    magma_int_t m, magma_int_t n,
    float alpha, float const *dA, magma_int_t lda,
                           float const *dx, magma_int_t incx,
    float beta,  float       *dy, magma_int_t incy );

void magma_ssymv(
    magma_uplo_t uplo,
    magma_int_t n,
    float alpha, float const *dA, magma_int_t lda,
                           float const *dx, magma_int_t incx,
    float beta,  float       *dy, magma_int_t incy );

void magma_strsv(
    magma_uplo_t uplo, magma_trans_t trans, magma_diag_t diag, 
    magma_int_t n, 
    float const *dA, magma_int_t lda, 
    float       *dx, magma_int_t incx );

// ========================================
// Level 3 BLAS

void magma_sgemm(
    magma_trans_t transA, magma_trans_t transB,
    magma_int_t m, magma_int_t n, magma_int_t k,
    float alpha, float const *dA, magma_int_t lda,
                           float const *dB, magma_int_t ldb,
    float beta,  float       *dC, magma_int_t ldc );

void magma_ssymm(
    magma_side_t side, magma_uplo_t uplo,
    magma_int_t m, magma_int_t n,
    float alpha, float const *dA, magma_int_t lda,
                           float const *dB, magma_int_t ldb,
    float beta,  float       *dC, magma_int_t ldc );

void magma_ssyrk(
    magma_uplo_t uplo, magma_trans_t trans,
    magma_int_t n, magma_int_t k,
    float alpha, float const *dA, magma_int_t lda,
    float beta,  float       *dC, magma_int_t ldc );

void magma_ssyr2k(
    magma_uplo_t uplo, magma_trans_t trans,
    magma_int_t n, magma_int_t k,
    float alpha, float const *dA, magma_int_t lda,
                           float const *dB, magma_int_t ldb,
    float beta,           float       *dC, magma_int_t ldc );

void magma_strmm(
    magma_side_t side, magma_uplo_t uplo, magma_trans_t trans, magma_diag_t diag,
    magma_int_t m, magma_int_t n,
    float alpha, float const *dA, magma_int_t lda,
                           float       *dB, magma_int_t ldb );

void magma_strsm(
    magma_side_t side, magma_uplo_t uplo, magma_trans_t trans, magma_diag_t diag,
    magma_int_t m, magma_int_t n,
    float alpha, float const *dA, magma_int_t lda,
                           float       *dB, magma_int_t ldb );

#ifdef __cplusplus
}
#endif

#undef PRECISION_s
#endif
