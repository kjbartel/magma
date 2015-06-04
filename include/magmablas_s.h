/*
 *   -- MAGMA (version 1.3.0) --
 *      Univ. of Tennessee, Knoxville
 *      Univ. of California, Berkeley
 *      Univ. of Colorado, Denver
 *      November 2012
 *
 * @generated s Wed Nov 14 22:52:26 2012
 */

#ifndef MAGMABLAS_S_H
#define MAGMABLAS_S_H

#define PRECISION_s

#ifdef __cplusplus
extern "C" {
#endif

  /*
   * Interface to clean
   */
float cpu_gpu_sdiff(
    magma_int_t m, magma_int_t n,
    const float *hA, magma_int_t lda,
    const float *dA, magma_int_t ldda );

// see also slaset
void szero_32x32_block(
    float *dA, magma_int_t ldda );

void szero_nbxnb_block(
    magma_int_t nb, float *dA, magma_int_t ldda );

// see also slaswp
void magmablas_spermute_long2(
    magma_int_t n, float *dAT, magma_int_t ldda,
    magma_int_t *ipiv, magma_int_t nb, magma_int_t ind );

void magmablas_spermute_long3(
    /*magma_int_t n,*/ float *dAT, magma_int_t ldda,
    const magma_int_t *ipiv, magma_int_t nb, magma_int_t ind );

  /*
   * Transpose functions
   */
void magmablas_sinplace_transpose(
    float *dA, magma_int_t ldda, magma_int_t n );

void magmablas_stranspose(
    float       *odata, magma_int_t ldo,
    const float *idata, magma_int_t ldi,
    magma_int_t m, magma_int_t n );

void magmablas_stranspose2(
    float       *odata, magma_int_t ldo,
    const float *idata, magma_int_t ldi,
    magma_int_t m, magma_int_t n );

void magmablas_stranspose2s(
    float       *odata, magma_int_t ldo,
    const float *idata, magma_int_t ldi,
    magma_int_t m, magma_int_t n,
    cudaStream_t *stream );

void magmablas_sgetmatrix_transpose(
    magma_int_t m, magma_int_t n,
    const float *dAT,   magma_int_t ldda,
    float       *hA,    magma_int_t lda,
    float       *dwork, magma_int_t lddwork, magma_int_t nb );

void magmablas_ssetmatrix_transpose(
    magma_int_t m, magma_int_t n,
    const float *hA,    magma_int_t lda,
    float       *dAT,   magma_int_t ldda,
    float       *dwork, magma_int_t lddwork, magma_int_t nb );

  /*
   * Multi-GPU functions
   */
void magmablas_sgetmatrix_transpose_mgpu(
    magma_int_t ngpu, cudaStream_t stream[][2],
    float **dAT, magma_int_t ldda,
    float  *hA,  magma_int_t lda,
    float **dB,  magma_int_t lddb,
    magma_int_t m, magma_int_t n, magma_int_t nb );

void magmablas_ssetmatrix_transpose_mgpu(
    magma_int_t ngpu, cudaStream_t stream[][2],
    const float  *hA,  magma_int_t lda,
    float       **dAT, magma_int_t ldda,
    float       **dB,  magma_int_t lddb,
    magma_int_t m, magma_int_t n, magma_int_t nb );

void magmablas_sgetmatrix_1D_bcyclic(
    magma_int_t m, magma_int_t n,
    float *dA[], magma_int_t ldda,
    float       *hA,   magma_int_t lda,
    magma_int_t ngpu, magma_int_t nb );

void magmablas_ssetmatrix_1D_bcyclic(
    magma_int_t m, magma_int_t n,
    const float *hA,   magma_int_t lda,
    float       *dA[], magma_int_t ldda,
    magma_int_t ngpu, magma_int_t nb );

void magmablas_ssymm_1gpu_old(
    char side, char uplo, magma_int_t m, magma_int_t n,
    float alpha,
    float *dA[], magma_int_t ldda,  magma_int_t offset,
    float *dB[], magma_int_t lddb,
    float beta,
    float *dC[], magma_int_t lddc,
    float *C,    magma_int_t ldc,
    magma_int_t ngpu, magma_int_t nb,
    cudaStream_t streams[][20], magma_int_t nstream );

void magmablas_ssymm_1gpu(
    char side, char uplo, magma_int_t m, magma_int_t n,
    float alpha,
    float *dA[], magma_int_t ldda,  magma_int_t offset,
    float *dB[], magma_int_t lddb,
    float beta,
    float *dC[], magma_int_t lddc,
    float *C,    magma_int_t ldc,
    magma_int_t ngpu, magma_int_t nb,
    cudaStream_t streams[][20], magma_int_t nstream );


void magmablas_ssymm_mgpu(
    char side, char uplo, magma_int_t m, magma_int_t n,
    float alpha, float *dA[], magma_int_t ldda,  magma_int_t offset,
                           float *dB[], magma_int_t lddb,
    float beta,  float *dC[], magma_int_t lddc,
                           float *dwork[],    magma_int_t lddwork,
                           float *C,    magma_int_t ldc,
                           float *work[], magma_int_t ldwork,
                           magma_int_t ngpu, magma_int_t nb, 
                           cudaStream_t streams[][20], magma_int_t nstream, 
                           cudaEvent_t redevents[][20],magma_int_t nbevents );

void magmablas_ssymm_mgpu_com(
    char side, char uplo, magma_int_t m, magma_int_t n,
    float alpha, float *dA[], magma_int_t ldda,  magma_int_t offset,
                           float *dB[], magma_int_t lddb,
    float beta,  float *dC[], magma_int_t lddc,
                           float *dwork[],    magma_int_t lddwork,
                           float *C,    magma_int_t ldc,
                           float *work[], magma_int_t ldwork,
                           magma_int_t ngpu, magma_int_t nb, 
                           cudaStream_t streams[][20], magma_int_t nstream, 
                           cudaEvent_t redevents[][20],magma_int_t nbevents, 
                           magma_int_t gnode[MagmaMaxGPUs][MagmaMaxGPUs+2], magma_int_t nbcmplx );

void magmablas_ssymm_mgpu_spec(
    char side, char uplo, magma_int_t m, magma_int_t n,
    float alpha, float *dA[], magma_int_t ldda,  magma_int_t offset,
                           float *dB[], magma_int_t lddb,
    float beta,  float *dC[], magma_int_t lddc,
                           float *dwork[],    magma_int_t lddwork,
                           float *C,    magma_int_t ldc,
                           float *work[], magma_int_t ldwork,
                           magma_int_t ngpu, magma_int_t nb, 
                           cudaStream_t streams[][20], magma_int_t nstream, 
                           cudaEvent_t redevents[][MagmaMaxGPUs*MagmaMaxGPUs],magma_int_t nbevents, 
                           magma_int_t gnode[MagmaMaxGPUs][MagmaMaxGPUs+2], magma_int_t nbcmplx );



void magmablas_ssyr2k_mgpu2(
    char uplo, char trans, magma_int_t n, magma_int_t k,
    float alpha,
    float *dA[], magma_int_t ldda, magma_int_t aoff,
    float *dB[], magma_int_t lddb, magma_int_t boff,
    float beta,
    float *dC[], magma_int_t lddc, magma_int_t offset,
    magma_int_t ngpu, magma_int_t nb,
    cudaStream_t streams[][20], magma_int_t nstream );

void magmablas_ssyr2k_mgpu_spec(
    char uplo, char trans, magma_int_t n, magma_int_t k,
    float alpha, float *dA[], magma_int_t lda, magma_int_t aoff,
                           float *dB[], magma_int_t ldb, magma_int_t boff,
    float beta,           float *dC[], magma_int_t ldc,  magma_int_t offset,
    magma_int_t ngpu, magma_int_t nb, cudaStream_t streams[][20], magma_int_t nstream );

  /*
   * LAPACK auxiliary functions
   */
void magmablas_sgeadd(
    magma_int_t m, magma_int_t n,
    float alpha,
    const float *dA, magma_int_t ldda,
    float       *dB, magma_int_t lddb );

void magmablas_slacpy(
    char uplo,
    magma_int_t m, magma_int_t n,
    const float *dA, magma_int_t ldda,
    float       *dB, magma_int_t lddb );

float magmablas_slange(
    char norm,
    magma_int_t m, magma_int_t n,
    const float *dA, magma_int_t ldda, float *dwork );

float magmablas_slansy(
    char norm, char uplo,
    magma_int_t n,
    const float *dA, magma_int_t ldda, float *dwork );

float magmablas_slansy(
    char norm, char uplo,
    magma_int_t n,
    const float *dA, magma_int_t ldda, float *dwork );

void magmablas_slascl(
    char type, magma_int_t kl, magma_int_t ku,
    float cfrom, float cto,
    magma_int_t m, magma_int_t n,
    float *dA, magma_int_t ldda, magma_int_t *info );

void magmablas_slaset(
    char uplo, magma_int_t m, magma_int_t n,
    float *dA, magma_int_t ldda );

void magmablas_slaset_identity(
    magma_int_t m, magma_int_t n,
    float *dA, magma_int_t ldda );

void magmablas_slaswp(
    magma_int_t n,
    float *dAT, magma_int_t ldda,
    magma_int_t i1,  magma_int_t i2,
    const magma_int_t *ipiv, magma_int_t inci );

void magmablas_slaswpx(
    magma_int_t n,
    float *dAT, magma_int_t ldx, magma_int_t ldy,
    magma_int_t i1, magma_int_t i2,
    const magma_int_t *ipiv, magma_int_t inci );

void magmablas_slaswp2(
    magma_int_t n,
    float* dAT, magma_int_t ldda,
    magma_int_t i1, magma_int_t i2,
    const magma_int_t *d_ipiv );

void magmablas_ssymmetrize(
    char uplo, magma_int_t m, float *dA, magma_int_t ldda );

void magmablas_ssymmetrize_tiles(
    char uplo, magma_int_t m, float *dA, magma_int_t ldda,
    magma_int_t ntile, magma_int_t mstride, magma_int_t nstride );

  /*
   * Level 1 BLAS
   */
void magmablas_sswap(
    magma_int_t n,
    float *dA, magma_int_t ldda,
    float *dB, magma_int_t lddb );

void magmablas_sswapblk(
    char storev,
    magma_int_t n,
    float *dA, magma_int_t ldda,
    float *dB, magma_int_t lddb,
    magma_int_t i1, magma_int_t i2,
    const magma_int_t *ipiv, magma_int_t inci,
    magma_int_t offset );

void magmablas_sswapdblk(
    magma_int_t n, magma_int_t nb,
    float *dA, magma_int_t ldda, magma_int_t inca,
    float *dB, magma_int_t lddb, magma_int_t incb );

  /*
   * Level 2 BLAS
   */
void magmablas_sgemv(
    char t, magma_int_t m, magma_int_t n,
    float alpha,
    const float *dA, magma_int_t ldda,
    const float *dx, magma_int_t incx,
    float beta,
    float       *dy, magma_int_t incy );

#if defined(PRECISION_z ) || defined(PRECISION_c )
magma_int_t magmablas_ssymv(
    char u, magma_int_t n,
    float alpha,
    const float *dA, magma_int_t ldda,
    const float *dx, magma_int_t incx,
    float beta,
    float       *dy, magma_int_t incy );
#endif

magma_int_t magmablas_ssymv(
    char u, magma_int_t n,
    float alpha,
    const float *dA, magma_int_t ldda,
    const float *dx, magma_int_t incx,
    float beta,
    float       *dy, magma_int_t incy );

  /*
   * Level 3 BLAS
   */
void magmablas_sgemm(
    char tA, char tB,
    magma_int_t m, magma_int_t n, magma_int_t k,
    float alpha,
    const float *dA, magma_int_t ldda,
    const float *dB, magma_int_t lddb,
    float beta,
    float       *dC, magma_int_t lddc );

void magmablas_sgemm_fermi80(
    char tA, char tB,
    magma_int_t m, magma_int_t n, magma_int_t k,
    float alpha,
    const float *dA, magma_int_t ldda,
    const float *dB, magma_int_t lddb,
    float beta,
    float       *dC, magma_int_t lddc );

void magmablas_sgemm_fermi64(
    char tA, char tB,
    magma_int_t m, magma_int_t n, magma_int_t k,
    float alpha,
    const float *dA, magma_int_t ldda,
    const float *dB, magma_int_t lddb,
    float beta,
    float       *dC, magma_int_t lddc );

void magmablas_ssymm(
    char s, char u,
    magma_int_t m, magma_int_t n,
    float alpha,
    const float *dA, magma_int_t ldda,
    const float *dB, magma_int_t lddb,
    float beta,
    float       *dC, magma_int_t lddc );

void magmablas_ssymm(
    char s, char u,
    magma_int_t m, magma_int_t n,
    float alpha,
    const float *dA, magma_int_t ldda,
    const float *dB, magma_int_t lddb,
    float beta,
    float       *dC, magma_int_t lddc );

void magmablas_ssyrk(
    char u, char t,
    magma_int_t n, magma_int_t k,
    float alpha,
    const float *dA, magma_int_t ldda,
    float beta,
    float       *dC, magma_int_t lddc );

void magmablas_ssyrk(
    char u, char t,
    magma_int_t n, magma_int_t k,
    float  alpha,
    const float *dA, magma_int_t ldda,
    float  beta,
    float       *dC, magma_int_t lddc );

void magmablas_ssyr2k(
    char u, char t,
    magma_int_t n, magma_int_t k,
    float alpha,
    const float *dA, magma_int_t ldda,
    const float *dB, magma_int_t lddb,
    float beta,
    float       *dC, magma_int_t lddc );

void magmablas_ssyr2k(
    char u, char t,
    magma_int_t n, magma_int_t k,
    float alpha,
    const float *dA, magma_int_t ldda,
    const float *dB, magma_int_t lddb,
    float  beta,
    float       *dC, magma_int_t lddc );

void magmablas_strmm(
    char s, char u, char t,  char d,
    magma_int_t m, magma_int_t n,
    float alpha,
    const float *dA, magma_int_t ldda,
    float       *dB, magma_int_t lddb );

void magmablas_strsm(
    char s, char u, char t, char d,
    magma_int_t m, magma_int_t n,
    float alpha,
    const float *dA, magma_int_t ldda,
    float       *dB, magma_int_t lddb );


  /*
   * Wrappers for platform independence.
   * These wrap CUBLAS or AMD OpenCL BLAS functions.
   */

// ========================================
// copying vectors
// set copies host to device
// get copies device to host
// Add the function, file, and line for error-reporting purposes.

#define magma_ssetvector(           n, hx_src, incx, dy_dst, incy ) \
        magma_ssetvector_internal(  n, hx_src, incx, dy_dst, incy, __func__, __FILE__, __LINE__ )

#define magma_sgetvector(           n, dx_src, incx, hy_dst, incy ) \
        magma_sgetvector_internal(  n, dx_src, incx, hy_dst, incy, __func__, __FILE__, __LINE__ )

#define magma_ssetvector_async(           n, hx_src, incx, dy_dst, incy, stream ) \
        magma_ssetvector_async_internal(  n, hx_src, incx, dy_dst, incy, stream, __func__, __FILE__, __LINE__ )

#define magma_sgetvector_async(           n, dx_src, incx, hy_dst, incy, stream ) \
        magma_sgetvector_async_internal(  n, dx_src, incx, hy_dst, incy, stream, __func__, __FILE__, __LINE__ )

void magma_ssetvector_internal(
    magma_int_t n,
    const float *hx_src, magma_int_t incx,
    float       *dy_dst, magma_int_t incy,
    const char* func, const char* file, int line );

void magma_sgetvector_internal(
    magma_int_t n,
    const float *dx_src, magma_int_t incx,
    float       *hy_dst, magma_int_t incy,
    const char* func, const char* file, int line );

void magma_ssetvector_async_internal(
    magma_int_t n,
    const float *hx_src, magma_int_t incx,
    float       *dy_dst, magma_int_t incy,
    magma_stream_t stream,
    const char* func, const char* file, int line );

void magma_sgetvector_async_internal(
    magma_int_t n,
    const float *dx_src, magma_int_t incx,
    float       *hy_dst, magma_int_t incy,
    magma_stream_t stream,
    const char* func, const char* file, int line );


// ========================================
// copying sub-matrices (contiguous columns )
// set  copies host to device
// get  copies device to host
// copy copies device to device
// (with CUDA unified addressing, copy can be between same or different devices )
// Add the function, file, and line for error-reporting purposes.

#define magma_ssetmatrix(           m, n, hA_src, lda, dB_dst, lddb ) \
        magma_ssetmatrix_internal(  m, n, hA_src, lda, dB_dst, lddb, __func__, __FILE__, __LINE__ )

#define magma_sgetmatrix(           m, n, dA_src, ldda, hB_dst, ldb ) \
        magma_sgetmatrix_internal(  m, n, dA_src, ldda, hB_dst, ldb, __func__, __FILE__, __LINE__ )

#define magma_scopymatrix(          m, n, dA_src, ldda, dB_dst, lddb ) \
        magma_scopymatrix_internal( m, n, dA_src, ldda, dB_dst, lddb, __func__, __FILE__, __LINE__ )

#define magma_ssetmatrix_async(           m, n, hA_src, lda, dB_dst, lddb, stream ) \
        magma_ssetmatrix_async_internal(  m, n, hA_src, lda, dB_dst, lddb, stream, __func__, __FILE__, __LINE__ )

#define magma_sgetmatrix_async(           m, n, dA_src, ldda, hB_dst, ldb, stream ) \
        magma_sgetmatrix_async_internal(  m, n, dA_src, ldda, hB_dst, ldb, stream, __func__, __FILE__, __LINE__ )

#define magma_scopymatrix_async(          m, n, dA_src, ldda, dB_dst, lddb, stream ) \
        magma_scopymatrix_async_internal( m, n, dA_src, ldda, dB_dst, lddb, stream, __func__, __FILE__, __LINE__ )

void magma_ssetmatrix_internal(
    magma_int_t m, magma_int_t n,
    const float *hA_src, magma_int_t lda,
    float       *dB_dst, magma_int_t lddb,
    const char* func, const char* file, int line );

void magma_sgetmatrix_internal(
    magma_int_t m, magma_int_t n,
    const float *dA_src, magma_int_t ldda,
    float       *hB_dst, magma_int_t ldb,
    const char* func, const char* file, int line );

void magma_scopymatrix_internal(
    magma_int_t m, magma_int_t n,
    const float *dA_src, magma_int_t ldda,
    float       *dB_dst, magma_int_t lddb,
    const char* func, const char* file, int line );

void magma_ssetmatrix_async_internal(
    magma_int_t m, magma_int_t n,
    const float *hA_src, magma_int_t lda,
    float       *dB_dst, magma_int_t lddb,
    magma_stream_t stream,
    const char* func, const char* file, int line );

void magma_sgetmatrix_async_internal(
    magma_int_t m, magma_int_t n,
    const float *dA_src, magma_int_t ldda,
    float       *hB_dst, magma_int_t ldb,
    magma_stream_t stream,
    const char* func, const char* file, int line );

void magma_scopymatrix_async_internal(
    magma_int_t m, magma_int_t n,
    const float *dA_src, magma_int_t ldda,
    float       *dB_dst, magma_int_t lddb,
    magma_stream_t stream,
    const char* func, const char* file, int line );


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
    float alpha,
    const float *dA, magma_int_t ldda,
    const float *dx, magma_int_t incx,
    float beta,
    float       *dy, magma_int_t incy );

void magma_ssymv(
    magma_uplo_t uplo,
    magma_int_t n,
    float alpha,
    const float *dA, magma_int_t ldda,
    const float *dx, magma_int_t incx,
    float beta,
    float       *dy, magma_int_t incy );

void magma_strsv(
    magma_uplo_t uplo, magma_trans_t trans, magma_diag_t diag,
    magma_int_t n,
    const float *dA, magma_int_t ldda,
    float       *dx, magma_int_t incx );

// ========================================
// Level 3 BLAS

void magma_sgemm(
    magma_trans_t transA, magma_trans_t transB,
    magma_int_t m, magma_int_t n, magma_int_t k,
    float alpha,
    const float *dA, magma_int_t ldda,
    const float *dB, magma_int_t lddb,
    float beta,
    float       *dC, magma_int_t lddc );

void magma_ssymm(
    magma_side_t side, magma_uplo_t uplo,
    magma_int_t m, magma_int_t n,
    float alpha,
    const float *dA, magma_int_t ldda,
    const float *dB, magma_int_t lddb,
    float beta,
    float       *dC, magma_int_t lddc );

void magma_ssyrk(
    magma_uplo_t uplo, magma_trans_t trans,
    magma_int_t n, magma_int_t k,
    float alpha,
    const float *dA, magma_int_t ldda,
    float beta,
    float       *dC, magma_int_t lddc );

void magma_ssyr2k(
    magma_uplo_t uplo, magma_trans_t trans,
    magma_int_t n, magma_int_t k,
    float alpha,
    const float *dA, magma_int_t ldda,
    const float *dB, magma_int_t lddb,
    float beta,
    float       *dC, magma_int_t lddc );

void magma_strmm(
    magma_side_t side, magma_uplo_t uplo, magma_trans_t trans, magma_diag_t diag,
    magma_int_t m, magma_int_t n,
    float alpha,
    const float *dA, magma_int_t ldda,
    float       *dB, magma_int_t lddb );

void magma_strsm(
    magma_side_t side, magma_uplo_t uplo, magma_trans_t trans, magma_diag_t diag,
    magma_int_t m, magma_int_t n,
    float alpha,
    const float *dA, magma_int_t ldda,
    float       *dB, magma_int_t lddb );

#ifdef __cplusplus
}
#endif

#undef PRECISION_s

#endif  // MAGMABLAS_S_H
