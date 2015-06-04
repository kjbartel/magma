/*
 *   -- MAGMA (version 1.3.0) --
 *      Univ. of Tennessee, Knoxville
 *      Univ. of California, Berkeley
 *      Univ. of Colorado, Denver
 *      November 2012
 *
 * @generated c Wed Nov 14 22:52:26 2012
 */

#ifndef MAGMABLAS_C_H
#define MAGMABLAS_C_H

#define PRECISION_c

#ifdef __cplusplus
extern "C" {
#endif

  /*
   * Interface to clean
   */
float cpu_gpu_cdiff(
    magma_int_t m, magma_int_t n,
    const cuFloatComplex *hA, magma_int_t lda,
    const cuFloatComplex *dA, magma_int_t ldda );

// see also claset
void czero_32x32_block(
    cuFloatComplex *dA, magma_int_t ldda );

void czero_nbxnb_block(
    magma_int_t nb, cuFloatComplex *dA, magma_int_t ldda );

// see also claswp
void magmablas_cpermute_long2(
    magma_int_t n, cuFloatComplex *dAT, magma_int_t ldda,
    magma_int_t *ipiv, magma_int_t nb, magma_int_t ind );

void magmablas_cpermute_long3(
    /*magma_int_t n,*/ cuFloatComplex *dAT, magma_int_t ldda,
    const magma_int_t *ipiv, magma_int_t nb, magma_int_t ind );

  /*
   * Transpose functions
   */
void magmablas_cinplace_transpose(
    cuFloatComplex *dA, magma_int_t ldda, magma_int_t n );

void magmablas_ctranspose(
    cuFloatComplex       *odata, magma_int_t ldo,
    const cuFloatComplex *idata, magma_int_t ldi,
    magma_int_t m, magma_int_t n );

void magmablas_ctranspose2(
    cuFloatComplex       *odata, magma_int_t ldo,
    const cuFloatComplex *idata, magma_int_t ldi,
    magma_int_t m, magma_int_t n );

void magmablas_ctranspose2s(
    cuFloatComplex       *odata, magma_int_t ldo,
    const cuFloatComplex *idata, magma_int_t ldi,
    magma_int_t m, magma_int_t n,
    cudaStream_t *stream );

void magmablas_cgetmatrix_transpose(
    magma_int_t m, magma_int_t n,
    const cuFloatComplex *dAT,   magma_int_t ldda,
    cuFloatComplex       *hA,    magma_int_t lda,
    cuFloatComplex       *dwork, magma_int_t lddwork, magma_int_t nb );

void magmablas_csetmatrix_transpose(
    magma_int_t m, magma_int_t n,
    const cuFloatComplex *hA,    magma_int_t lda,
    cuFloatComplex       *dAT,   magma_int_t ldda,
    cuFloatComplex       *dwork, magma_int_t lddwork, magma_int_t nb );

  /*
   * Multi-GPU functions
   */
void magmablas_cgetmatrix_transpose_mgpu(
    magma_int_t ngpu, cudaStream_t stream[][2],
    cuFloatComplex **dAT, magma_int_t ldda,
    cuFloatComplex  *hA,  magma_int_t lda,
    cuFloatComplex **dB,  magma_int_t lddb,
    magma_int_t m, magma_int_t n, magma_int_t nb );

void magmablas_csetmatrix_transpose_mgpu(
    magma_int_t ngpu, cudaStream_t stream[][2],
    const cuFloatComplex  *hA,  magma_int_t lda,
    cuFloatComplex       **dAT, magma_int_t ldda,
    cuFloatComplex       **dB,  magma_int_t lddb,
    magma_int_t m, magma_int_t n, magma_int_t nb );

void magmablas_cgetmatrix_1D_bcyclic(
    magma_int_t m, magma_int_t n,
    cuFloatComplex *dA[], magma_int_t ldda,
    cuFloatComplex       *hA,   magma_int_t lda,
    magma_int_t ngpu, magma_int_t nb );

void magmablas_csetmatrix_1D_bcyclic(
    magma_int_t m, magma_int_t n,
    const cuFloatComplex *hA,   magma_int_t lda,
    cuFloatComplex       *dA[], magma_int_t ldda,
    magma_int_t ngpu, magma_int_t nb );

void magmablas_chemm_1gpu_old(
    char side, char uplo, magma_int_t m, magma_int_t n,
    cuFloatComplex alpha,
    cuFloatComplex *dA[], magma_int_t ldda,  magma_int_t offset,
    cuFloatComplex *dB[], magma_int_t lddb,
    cuFloatComplex beta,
    cuFloatComplex *dC[], magma_int_t lddc,
    cuFloatComplex *C,    magma_int_t ldc,
    magma_int_t ngpu, magma_int_t nb,
    cudaStream_t streams[][20], magma_int_t nstream );

void magmablas_chemm_1gpu(
    char side, char uplo, magma_int_t m, magma_int_t n,
    cuFloatComplex alpha,
    cuFloatComplex *dA[], magma_int_t ldda,  magma_int_t offset,
    cuFloatComplex *dB[], magma_int_t lddb,
    cuFloatComplex beta,
    cuFloatComplex *dC[], magma_int_t lddc,
    cuFloatComplex *C,    magma_int_t ldc,
    magma_int_t ngpu, magma_int_t nb,
    cudaStream_t streams[][20], magma_int_t nstream );


void magmablas_chemm_mgpu(
    char side, char uplo, magma_int_t m, magma_int_t n,
    cuFloatComplex alpha, cuFloatComplex *dA[], magma_int_t ldda,  magma_int_t offset,
                           cuFloatComplex *dB[], magma_int_t lddb,
    cuFloatComplex beta,  cuFloatComplex *dC[], magma_int_t lddc,
                           cuFloatComplex *dwork[],    magma_int_t lddwork,
                           cuFloatComplex *C,    magma_int_t ldc,
                           cuFloatComplex *work[], magma_int_t ldwork,
                           magma_int_t ngpu, magma_int_t nb, 
                           cudaStream_t streams[][20], magma_int_t nstream, 
                           cudaEvent_t redevents[][20],magma_int_t nbevents );

void magmablas_chemm_mgpu_com(
    char side, char uplo, magma_int_t m, magma_int_t n,
    cuFloatComplex alpha, cuFloatComplex *dA[], magma_int_t ldda,  magma_int_t offset,
                           cuFloatComplex *dB[], magma_int_t lddb,
    cuFloatComplex beta,  cuFloatComplex *dC[], magma_int_t lddc,
                           cuFloatComplex *dwork[],    magma_int_t lddwork,
                           cuFloatComplex *C,    magma_int_t ldc,
                           cuFloatComplex *work[], magma_int_t ldwork,
                           magma_int_t ngpu, magma_int_t nb, 
                           cudaStream_t streams[][20], magma_int_t nstream, 
                           cudaEvent_t redevents[][20],magma_int_t nbevents, 
                           magma_int_t gnode[MagmaMaxGPUs][MagmaMaxGPUs+2], magma_int_t nbcmplx );

void magmablas_chemm_mgpu_spec(
    char side, char uplo, magma_int_t m, magma_int_t n,
    cuFloatComplex alpha, cuFloatComplex *dA[], magma_int_t ldda,  magma_int_t offset,
                           cuFloatComplex *dB[], magma_int_t lddb,
    cuFloatComplex beta,  cuFloatComplex *dC[], magma_int_t lddc,
                           cuFloatComplex *dwork[],    magma_int_t lddwork,
                           cuFloatComplex *C,    magma_int_t ldc,
                           cuFloatComplex *work[], magma_int_t ldwork,
                           magma_int_t ngpu, magma_int_t nb, 
                           cudaStream_t streams[][20], magma_int_t nstream, 
                           cudaEvent_t redevents[][MagmaMaxGPUs*MagmaMaxGPUs],magma_int_t nbevents, 
                           magma_int_t gnode[MagmaMaxGPUs][MagmaMaxGPUs+2], magma_int_t nbcmplx );



void magmablas_cher2k_mgpu2(
    char uplo, char trans, magma_int_t n, magma_int_t k,
    cuFloatComplex alpha,
    cuFloatComplex *dA[], magma_int_t ldda, magma_int_t aoff,
    cuFloatComplex *dB[], magma_int_t lddb, magma_int_t boff,
    float beta,
    cuFloatComplex *dC[], magma_int_t lddc, magma_int_t offset,
    magma_int_t ngpu, magma_int_t nb,
    cudaStream_t streams[][20], magma_int_t nstream );

void magmablas_cher2k_mgpu_spec(
    char uplo, char trans, magma_int_t n, magma_int_t k,
    cuFloatComplex alpha, cuFloatComplex *dA[], magma_int_t lda, magma_int_t aoff,
                           cuFloatComplex *dB[], magma_int_t ldb, magma_int_t boff,
    float beta,           cuFloatComplex *dC[], magma_int_t ldc,  magma_int_t offset,
    magma_int_t ngpu, magma_int_t nb, cudaStream_t streams[][20], magma_int_t nstream );

  /*
   * LAPACK auxiliary functions
   */
void magmablas_cgeadd(
    magma_int_t m, magma_int_t n,
    cuFloatComplex alpha,
    const cuFloatComplex *dA, magma_int_t ldda,
    cuFloatComplex       *dB, magma_int_t lddb );

void magmablas_clacpy(
    char uplo,
    magma_int_t m, magma_int_t n,
    const cuFloatComplex *dA, magma_int_t ldda,
    cuFloatComplex       *dB, magma_int_t lddb );

float magmablas_clange(
    char norm,
    magma_int_t m, magma_int_t n,
    const cuFloatComplex *dA, magma_int_t ldda, float *dwork );

float magmablas_clanhe(
    char norm, char uplo,
    magma_int_t n,
    const cuFloatComplex *dA, magma_int_t ldda, float *dwork );

float magmablas_clansy(
    char norm, char uplo,
    magma_int_t n,
    const cuFloatComplex *dA, magma_int_t ldda, float *dwork );

void magmablas_clascl(
    char type, magma_int_t kl, magma_int_t ku,
    float cfrom, float cto,
    magma_int_t m, magma_int_t n,
    cuFloatComplex *dA, magma_int_t ldda, magma_int_t *info );

void magmablas_claset(
    char uplo, magma_int_t m, magma_int_t n,
    cuFloatComplex *dA, magma_int_t ldda );

void magmablas_claset_identity(
    magma_int_t m, magma_int_t n,
    cuFloatComplex *dA, magma_int_t ldda );

void magmablas_claswp(
    magma_int_t n,
    cuFloatComplex *dAT, magma_int_t ldda,
    magma_int_t i1,  magma_int_t i2,
    const magma_int_t *ipiv, magma_int_t inci );

void magmablas_claswpx(
    magma_int_t n,
    cuFloatComplex *dAT, magma_int_t ldx, magma_int_t ldy,
    magma_int_t i1, magma_int_t i2,
    const magma_int_t *ipiv, magma_int_t inci );

void magmablas_claswp2(
    magma_int_t n,
    cuFloatComplex* dAT, magma_int_t ldda,
    magma_int_t i1, magma_int_t i2,
    const magma_int_t *d_ipiv );

void magmablas_csymmetrize(
    char uplo, magma_int_t m, cuFloatComplex *dA, magma_int_t ldda );

void magmablas_csymmetrize_tiles(
    char uplo, magma_int_t m, cuFloatComplex *dA, magma_int_t ldda,
    magma_int_t ntile, magma_int_t mstride, magma_int_t nstride );

  /*
   * Level 1 BLAS
   */
void magmablas_cswap(
    magma_int_t n,
    cuFloatComplex *dA, magma_int_t ldda,
    cuFloatComplex *dB, magma_int_t lddb );

void magmablas_cswapblk(
    char storev,
    magma_int_t n,
    cuFloatComplex *dA, magma_int_t ldda,
    cuFloatComplex *dB, magma_int_t lddb,
    magma_int_t i1, magma_int_t i2,
    const magma_int_t *ipiv, magma_int_t inci,
    magma_int_t offset );

void magmablas_cswapdblk(
    magma_int_t n, magma_int_t nb,
    cuFloatComplex *dA, magma_int_t ldda, magma_int_t inca,
    cuFloatComplex *dB, magma_int_t lddb, magma_int_t incb );

  /*
   * Level 2 BLAS
   */
void magmablas_cgemv(
    char t, magma_int_t m, magma_int_t n,
    cuFloatComplex alpha,
    const cuFloatComplex *dA, magma_int_t ldda,
    const cuFloatComplex *dx, magma_int_t incx,
    cuFloatComplex beta,
    cuFloatComplex       *dy, magma_int_t incy );

#if defined(PRECISION_z ) || defined(PRECISION_c )
magma_int_t magmablas_chemv(
    char u, magma_int_t n,
    cuFloatComplex alpha,
    const cuFloatComplex *dA, magma_int_t ldda,
    const cuFloatComplex *dx, magma_int_t incx,
    cuFloatComplex beta,
    cuFloatComplex       *dy, magma_int_t incy );
#endif

magma_int_t magmablas_csymv(
    char u, magma_int_t n,
    cuFloatComplex alpha,
    const cuFloatComplex *dA, magma_int_t ldda,
    const cuFloatComplex *dx, magma_int_t incx,
    cuFloatComplex beta,
    cuFloatComplex       *dy, magma_int_t incy );

  /*
   * Level 3 BLAS
   */
void magmablas_cgemm(
    char tA, char tB,
    magma_int_t m, magma_int_t n, magma_int_t k,
    cuFloatComplex alpha,
    const cuFloatComplex *dA, magma_int_t ldda,
    const cuFloatComplex *dB, magma_int_t lddb,
    cuFloatComplex beta,
    cuFloatComplex       *dC, magma_int_t lddc );

void magmablas_cgemm_fermi80(
    char tA, char tB,
    magma_int_t m, magma_int_t n, magma_int_t k,
    cuFloatComplex alpha,
    const cuFloatComplex *dA, magma_int_t ldda,
    const cuFloatComplex *dB, magma_int_t lddb,
    cuFloatComplex beta,
    cuFloatComplex       *dC, magma_int_t lddc );

void magmablas_cgemm_fermi64(
    char tA, char tB,
    magma_int_t m, magma_int_t n, magma_int_t k,
    cuFloatComplex alpha,
    const cuFloatComplex *dA, magma_int_t ldda,
    const cuFloatComplex *dB, magma_int_t lddb,
    cuFloatComplex beta,
    cuFloatComplex       *dC, magma_int_t lddc );

void magmablas_chemm(
    char s, char u,
    magma_int_t m, magma_int_t n,
    cuFloatComplex alpha,
    const cuFloatComplex *dA, magma_int_t ldda,
    const cuFloatComplex *dB, magma_int_t lddb,
    cuFloatComplex beta,
    cuFloatComplex       *dC, magma_int_t lddc );

void magmablas_csymm(
    char s, char u,
    magma_int_t m, magma_int_t n,
    cuFloatComplex alpha,
    const cuFloatComplex *dA, magma_int_t ldda,
    const cuFloatComplex *dB, magma_int_t lddb,
    cuFloatComplex beta,
    cuFloatComplex       *dC, magma_int_t lddc );

void magmablas_csyrk(
    char u, char t,
    magma_int_t n, magma_int_t k,
    cuFloatComplex alpha,
    const cuFloatComplex *dA, magma_int_t ldda,
    cuFloatComplex beta,
    cuFloatComplex       *dC, magma_int_t lddc );

void magmablas_cherk(
    char u, char t,
    magma_int_t n, magma_int_t k,
    float  alpha,
    const cuFloatComplex *dA, magma_int_t ldda,
    float  beta,
    cuFloatComplex       *dC, magma_int_t lddc );

void magmablas_csyr2k(
    char u, char t,
    magma_int_t n, magma_int_t k,
    cuFloatComplex alpha,
    const cuFloatComplex *dA, magma_int_t ldda,
    const cuFloatComplex *dB, magma_int_t lddb,
    cuFloatComplex beta,
    cuFloatComplex       *dC, magma_int_t lddc );

void magmablas_cher2k(
    char u, char t,
    magma_int_t n, magma_int_t k,
    cuFloatComplex alpha,
    const cuFloatComplex *dA, magma_int_t ldda,
    const cuFloatComplex *dB, magma_int_t lddb,
    float  beta,
    cuFloatComplex       *dC, magma_int_t lddc );

void magmablas_ctrmm(
    char s, char u, char t,  char d,
    magma_int_t m, magma_int_t n,
    cuFloatComplex alpha,
    const cuFloatComplex *dA, magma_int_t ldda,
    cuFloatComplex       *dB, magma_int_t lddb );

void magmablas_ctrsm(
    char s, char u, char t, char d,
    magma_int_t m, magma_int_t n,
    cuFloatComplex alpha,
    const cuFloatComplex *dA, magma_int_t ldda,
    cuFloatComplex       *dB, magma_int_t lddb );


  /*
   * Wrappers for platform independence.
   * These wrap CUBLAS or AMD OpenCL BLAS functions.
   */

// ========================================
// copying vectors
// set copies host to device
// get copies device to host
// Add the function, file, and line for error-reporting purposes.

#define magma_csetvector(           n, hx_src, incx, dy_dst, incy ) \
        magma_csetvector_internal(  n, hx_src, incx, dy_dst, incy, __func__, __FILE__, __LINE__ )

#define magma_cgetvector(           n, dx_src, incx, hy_dst, incy ) \
        magma_cgetvector_internal(  n, dx_src, incx, hy_dst, incy, __func__, __FILE__, __LINE__ )

#define magma_csetvector_async(           n, hx_src, incx, dy_dst, incy, stream ) \
        magma_csetvector_async_internal(  n, hx_src, incx, dy_dst, incy, stream, __func__, __FILE__, __LINE__ )

#define magma_cgetvector_async(           n, dx_src, incx, hy_dst, incy, stream ) \
        magma_cgetvector_async_internal(  n, dx_src, incx, hy_dst, incy, stream, __func__, __FILE__, __LINE__ )

void magma_csetvector_internal(
    magma_int_t n,
    const cuFloatComplex *hx_src, magma_int_t incx,
    cuFloatComplex       *dy_dst, magma_int_t incy,
    const char* func, const char* file, int line );

void magma_cgetvector_internal(
    magma_int_t n,
    const cuFloatComplex *dx_src, magma_int_t incx,
    cuFloatComplex       *hy_dst, magma_int_t incy,
    const char* func, const char* file, int line );

void magma_csetvector_async_internal(
    magma_int_t n,
    const cuFloatComplex *hx_src, magma_int_t incx,
    cuFloatComplex       *dy_dst, magma_int_t incy,
    magma_stream_t stream,
    const char* func, const char* file, int line );

void magma_cgetvector_async_internal(
    magma_int_t n,
    const cuFloatComplex *dx_src, magma_int_t incx,
    cuFloatComplex       *hy_dst, magma_int_t incy,
    magma_stream_t stream,
    const char* func, const char* file, int line );


// ========================================
// copying sub-matrices (contiguous columns )
// set  copies host to device
// get  copies device to host
// copy copies device to device
// (with CUDA unified addressing, copy can be between same or different devices )
// Add the function, file, and line for error-reporting purposes.

#define magma_csetmatrix(           m, n, hA_src, lda, dB_dst, lddb ) \
        magma_csetmatrix_internal(  m, n, hA_src, lda, dB_dst, lddb, __func__, __FILE__, __LINE__ )

#define magma_cgetmatrix(           m, n, dA_src, ldda, hB_dst, ldb ) \
        magma_cgetmatrix_internal(  m, n, dA_src, ldda, hB_dst, ldb, __func__, __FILE__, __LINE__ )

#define magma_ccopymatrix(          m, n, dA_src, ldda, dB_dst, lddb ) \
        magma_ccopymatrix_internal( m, n, dA_src, ldda, dB_dst, lddb, __func__, __FILE__, __LINE__ )

#define magma_csetmatrix_async(           m, n, hA_src, lda, dB_dst, lddb, stream ) \
        magma_csetmatrix_async_internal(  m, n, hA_src, lda, dB_dst, lddb, stream, __func__, __FILE__, __LINE__ )

#define magma_cgetmatrix_async(           m, n, dA_src, ldda, hB_dst, ldb, stream ) \
        magma_cgetmatrix_async_internal(  m, n, dA_src, ldda, hB_dst, ldb, stream, __func__, __FILE__, __LINE__ )

#define magma_ccopymatrix_async(          m, n, dA_src, ldda, dB_dst, lddb, stream ) \
        magma_ccopymatrix_async_internal( m, n, dA_src, ldda, dB_dst, lddb, stream, __func__, __FILE__, __LINE__ )

void magma_csetmatrix_internal(
    magma_int_t m, magma_int_t n,
    const cuFloatComplex *hA_src, magma_int_t lda,
    cuFloatComplex       *dB_dst, magma_int_t lddb,
    const char* func, const char* file, int line );

void magma_cgetmatrix_internal(
    magma_int_t m, magma_int_t n,
    const cuFloatComplex *dA_src, magma_int_t ldda,
    cuFloatComplex       *hB_dst, magma_int_t ldb,
    const char* func, const char* file, int line );

void magma_ccopymatrix_internal(
    magma_int_t m, magma_int_t n,
    const cuFloatComplex *dA_src, magma_int_t ldda,
    cuFloatComplex       *dB_dst, magma_int_t lddb,
    const char* func, const char* file, int line );

void magma_csetmatrix_async_internal(
    magma_int_t m, magma_int_t n,
    const cuFloatComplex *hA_src, magma_int_t lda,
    cuFloatComplex       *dB_dst, magma_int_t lddb,
    magma_stream_t stream,
    const char* func, const char* file, int line );

void magma_cgetmatrix_async_internal(
    magma_int_t m, magma_int_t n,
    const cuFloatComplex *dA_src, magma_int_t ldda,
    cuFloatComplex       *hB_dst, magma_int_t ldb,
    magma_stream_t stream,
    const char* func, const char* file, int line );

void magma_ccopymatrix_async_internal(
    magma_int_t m, magma_int_t n,
    const cuFloatComplex *dA_src, magma_int_t ldda,
    cuFloatComplex       *dB_dst, magma_int_t lddb,
    magma_stream_t stream,
    const char* func, const char* file, int line );


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
    cuFloatComplex alpha,
    const cuFloatComplex *dA, magma_int_t ldda,
    const cuFloatComplex *dx, magma_int_t incx,
    cuFloatComplex beta,
    cuFloatComplex       *dy, magma_int_t incy );

void magma_chemv(
    magma_uplo_t uplo,
    magma_int_t n,
    cuFloatComplex alpha,
    const cuFloatComplex *dA, magma_int_t ldda,
    const cuFloatComplex *dx, magma_int_t incx,
    cuFloatComplex beta,
    cuFloatComplex       *dy, magma_int_t incy );

void magma_ctrsv(
    magma_uplo_t uplo, magma_trans_t trans, magma_diag_t diag,
    magma_int_t n,
    const cuFloatComplex *dA, magma_int_t ldda,
    cuFloatComplex       *dx, magma_int_t incx );

// ========================================
// Level 3 BLAS

void magma_cgemm(
    magma_trans_t transA, magma_trans_t transB,
    magma_int_t m, magma_int_t n, magma_int_t k,
    cuFloatComplex alpha,
    const cuFloatComplex *dA, magma_int_t ldda,
    const cuFloatComplex *dB, magma_int_t lddb,
    cuFloatComplex beta,
    cuFloatComplex       *dC, magma_int_t lddc );

void magma_chemm(
    magma_side_t side, magma_uplo_t uplo,
    magma_int_t m, magma_int_t n,
    cuFloatComplex alpha,
    const cuFloatComplex *dA, magma_int_t ldda,
    const cuFloatComplex *dB, magma_int_t lddb,
    cuFloatComplex beta,
    cuFloatComplex       *dC, magma_int_t lddc );

void magma_cherk(
    magma_uplo_t uplo, magma_trans_t trans,
    magma_int_t n, magma_int_t k,
    float alpha,
    const cuFloatComplex *dA, magma_int_t ldda,
    float beta,
    cuFloatComplex       *dC, magma_int_t lddc );

void magma_cher2k(
    magma_uplo_t uplo, magma_trans_t trans,
    magma_int_t n, magma_int_t k,
    cuFloatComplex alpha,
    const cuFloatComplex *dA, magma_int_t ldda,
    const cuFloatComplex *dB, magma_int_t lddb,
    float beta,
    cuFloatComplex       *dC, magma_int_t lddc );

void magma_ctrmm(
    magma_side_t side, magma_uplo_t uplo, magma_trans_t trans, magma_diag_t diag,
    magma_int_t m, magma_int_t n,
    cuFloatComplex alpha,
    const cuFloatComplex *dA, magma_int_t ldda,
    cuFloatComplex       *dB, magma_int_t lddb );

void magma_ctrsm(
    magma_side_t side, magma_uplo_t uplo, magma_trans_t trans, magma_diag_t diag,
    magma_int_t m, magma_int_t n,
    cuFloatComplex alpha,
    const cuFloatComplex *dA, magma_int_t ldda,
    cuFloatComplex       *dB, magma_int_t lddb );

#ifdef __cplusplus
}
#endif

#undef PRECISION_c

#endif  // MAGMABLAS_C_H
