/*
 *   -- MAGMA (version 1.3.0) --
 *      Univ. of Tennessee, Knoxville
 *      Univ. of California, Berkeley
 *      Univ. of Colorado, Denver
 *      November 2012
 *
 * @precisions normal z -> s d c
 */

#ifndef MAGMABLAS_Z_H
#define MAGMABLAS_Z_H

#define PRECISION_z

#ifdef __cplusplus
extern "C" {
#endif

  /*
   * Interface to clean
   */
double cpu_gpu_zdiff(
    magma_int_t m, magma_int_t n,
    const cuDoubleComplex *hA, magma_int_t lda,
    const cuDoubleComplex *dA, magma_int_t ldda );

// see also zlaset
void zzero_32x32_block(
    cuDoubleComplex *dA, magma_int_t ldda );

void zzero_nbxnb_block(
    magma_int_t nb, cuDoubleComplex *dA, magma_int_t ldda );

// see also zlaswp
void magmablas_zpermute_long2(
    magma_int_t n, cuDoubleComplex *dAT, magma_int_t ldda,
    magma_int_t *ipiv, magma_int_t nb, magma_int_t ind );

void magmablas_zpermute_long3(
    /*magma_int_t n,*/ cuDoubleComplex *dAT, magma_int_t ldda,
    const magma_int_t *ipiv, magma_int_t nb, magma_int_t ind );

  /*
   * Transpose functions
   */
void magmablas_zinplace_transpose(
    cuDoubleComplex *dA, magma_int_t ldda, magma_int_t n );

void magmablas_ztranspose(
    cuDoubleComplex       *odata, magma_int_t ldo,
    const cuDoubleComplex *idata, magma_int_t ldi,
    magma_int_t m, magma_int_t n );

void magmablas_ztranspose2(
    cuDoubleComplex       *odata, magma_int_t ldo,
    const cuDoubleComplex *idata, magma_int_t ldi,
    magma_int_t m, magma_int_t n );

void magmablas_ztranspose2s(
    cuDoubleComplex       *odata, magma_int_t ldo,
    const cuDoubleComplex *idata, magma_int_t ldi,
    magma_int_t m, magma_int_t n,
    cudaStream_t *stream );

void magmablas_zgetmatrix_transpose(
    magma_int_t m, magma_int_t n,
    const cuDoubleComplex *dAT,   magma_int_t ldda,
    cuDoubleComplex       *hA,    magma_int_t lda,
    cuDoubleComplex       *dwork, magma_int_t lddwork, magma_int_t nb );

void magmablas_zsetmatrix_transpose(
    magma_int_t m, magma_int_t n,
    const cuDoubleComplex *hA,    magma_int_t lda,
    cuDoubleComplex       *dAT,   magma_int_t ldda,
    cuDoubleComplex       *dwork, magma_int_t lddwork, magma_int_t nb );

  /*
   * Multi-GPU functions
   */
void magmablas_zgetmatrix_transpose_mgpu(
    magma_int_t ngpu, cudaStream_t stream[][2],
    cuDoubleComplex **dAT, magma_int_t ldda,
    cuDoubleComplex  *hA,  magma_int_t lda,
    cuDoubleComplex **dB,  magma_int_t lddb,
    magma_int_t m, magma_int_t n, magma_int_t nb );

void magmablas_zsetmatrix_transpose_mgpu(
    magma_int_t ngpu, cudaStream_t stream[][2],
    const cuDoubleComplex  *hA,  magma_int_t lda,
    cuDoubleComplex       **dAT, magma_int_t ldda,
    cuDoubleComplex       **dB,  magma_int_t lddb,
    magma_int_t m, magma_int_t n, magma_int_t nb );

void magmablas_zgetmatrix_1D_bcyclic(
    magma_int_t m, magma_int_t n,
    cuDoubleComplex *dA[], magma_int_t ldda,
    cuDoubleComplex       *hA,   magma_int_t lda,
    magma_int_t ngpu, magma_int_t nb );

void magmablas_zsetmatrix_1D_bcyclic(
    magma_int_t m, magma_int_t n,
    const cuDoubleComplex *hA,   magma_int_t lda,
    cuDoubleComplex       *dA[], magma_int_t ldda,
    magma_int_t ngpu, magma_int_t nb );

void magmablas_zhemm_1gpu_old(
    char side, char uplo, magma_int_t m, magma_int_t n,
    cuDoubleComplex alpha,
    cuDoubleComplex *dA[], magma_int_t ldda,  magma_int_t offset,
    cuDoubleComplex *dB[], magma_int_t lddb,
    cuDoubleComplex beta,
    cuDoubleComplex *dC[], magma_int_t lddc,
    cuDoubleComplex *C,    magma_int_t ldc,
    magma_int_t ngpu, magma_int_t nb,
    cudaStream_t streams[][20], magma_int_t nstream );

void magmablas_zhemm_1gpu(
    char side, char uplo, magma_int_t m, magma_int_t n,
    cuDoubleComplex alpha,
    cuDoubleComplex *dA[], magma_int_t ldda,  magma_int_t offset,
    cuDoubleComplex *dB[], magma_int_t lddb,
    cuDoubleComplex beta,
    cuDoubleComplex *dC[], magma_int_t lddc,
    cuDoubleComplex *C,    magma_int_t ldc,
    magma_int_t ngpu, magma_int_t nb,
    cudaStream_t streams[][20], magma_int_t nstream );


void magmablas_zhemm_mgpu(
    char side, char uplo, magma_int_t m, magma_int_t n,
    cuDoubleComplex alpha, cuDoubleComplex *dA[], magma_int_t ldda,  magma_int_t offset,
                           cuDoubleComplex *dB[], magma_int_t lddb,
    cuDoubleComplex beta,  cuDoubleComplex *dC[], magma_int_t lddc,
                           cuDoubleComplex *dwork[],    magma_int_t lddwork,
                           cuDoubleComplex *C,    magma_int_t ldc,
                           cuDoubleComplex *work[], magma_int_t ldwork,
                           magma_int_t ngpu, magma_int_t nb, 
                           cudaStream_t streams[][20], magma_int_t nstream, 
                           cudaEvent_t redevents[][20],magma_int_t nbevents );

void magmablas_zhemm_mgpu_com(
    char side, char uplo, magma_int_t m, magma_int_t n,
    cuDoubleComplex alpha, cuDoubleComplex *dA[], magma_int_t ldda,  magma_int_t offset,
                           cuDoubleComplex *dB[], magma_int_t lddb,
    cuDoubleComplex beta,  cuDoubleComplex *dC[], magma_int_t lddc,
                           cuDoubleComplex *dwork[],    magma_int_t lddwork,
                           cuDoubleComplex *C,    magma_int_t ldc,
                           cuDoubleComplex *work[], magma_int_t ldwork,
                           magma_int_t ngpu, magma_int_t nb, 
                           cudaStream_t streams[][20], magma_int_t nstream, 
                           cudaEvent_t redevents[][20],magma_int_t nbevents, 
                           magma_int_t gnode[MagmaMaxGPUs][MagmaMaxGPUs+2], magma_int_t nbcmplx );

void magmablas_zhemm_mgpu_spec(
    char side, char uplo, magma_int_t m, magma_int_t n,
    cuDoubleComplex alpha, cuDoubleComplex *dA[], magma_int_t ldda,  magma_int_t offset,
                           cuDoubleComplex *dB[], magma_int_t lddb,
    cuDoubleComplex beta,  cuDoubleComplex *dC[], magma_int_t lddc,
                           cuDoubleComplex *dwork[],    magma_int_t lddwork,
                           cuDoubleComplex *C,    magma_int_t ldc,
                           cuDoubleComplex *work[], magma_int_t ldwork,
                           magma_int_t ngpu, magma_int_t nb, 
                           cudaStream_t streams[][20], magma_int_t nstream, 
                           cudaEvent_t redevents[][MagmaMaxGPUs*MagmaMaxGPUs],magma_int_t nbevents, 
                           magma_int_t gnode[MagmaMaxGPUs][MagmaMaxGPUs+2], magma_int_t nbcmplx );



void magmablas_zher2k_mgpu2(
    char uplo, char trans, magma_int_t n, magma_int_t k,
    cuDoubleComplex alpha,
    cuDoubleComplex *dA[], magma_int_t ldda, magma_int_t aoff,
    cuDoubleComplex *dB[], magma_int_t lddb, magma_int_t boff,
    double beta,
    cuDoubleComplex *dC[], magma_int_t lddc, magma_int_t offset,
    magma_int_t ngpu, magma_int_t nb,
    cudaStream_t streams[][20], magma_int_t nstream );

void magmablas_zher2k_mgpu_spec(
    char uplo, char trans, magma_int_t n, magma_int_t k,
    cuDoubleComplex alpha, cuDoubleComplex *dA[], magma_int_t lda, magma_int_t aoff,
                           cuDoubleComplex *dB[], magma_int_t ldb, magma_int_t boff,
    double beta,           cuDoubleComplex *dC[], magma_int_t ldc,  magma_int_t offset,
    magma_int_t ngpu, magma_int_t nb, cudaStream_t streams[][20], magma_int_t nstream );

  /*
   * LAPACK auxiliary functions
   */
void magmablas_zgeadd(
    magma_int_t m, magma_int_t n,
    cuDoubleComplex alpha,
    const cuDoubleComplex *dA, magma_int_t ldda,
    cuDoubleComplex       *dB, magma_int_t lddb );

void magmablas_zlacpy(
    char uplo,
    magma_int_t m, magma_int_t n,
    const cuDoubleComplex *dA, magma_int_t ldda,
    cuDoubleComplex       *dB, magma_int_t lddb );

double magmablas_zlange(
    char norm,
    magma_int_t m, magma_int_t n,
    const cuDoubleComplex *dA, magma_int_t ldda, double *dwork );

double magmablas_zlanhe(
    char norm, char uplo,
    magma_int_t n,
    const cuDoubleComplex *dA, magma_int_t ldda, double *dwork );

double magmablas_zlansy(
    char norm, char uplo,
    magma_int_t n,
    const cuDoubleComplex *dA, magma_int_t ldda, double *dwork );

void magmablas_zlascl(
    char type, magma_int_t kl, magma_int_t ku,
    double cfrom, double cto,
    magma_int_t m, magma_int_t n,
    cuDoubleComplex *dA, magma_int_t ldda, magma_int_t *info );

void magmablas_zlaset(
    char uplo, magma_int_t m, magma_int_t n,
    cuDoubleComplex *dA, magma_int_t ldda );

void magmablas_zlaset_identity(
    magma_int_t m, magma_int_t n,
    cuDoubleComplex *dA, magma_int_t ldda );

void magmablas_zlaswp(
    magma_int_t n,
    cuDoubleComplex *dAT, magma_int_t ldda,
    magma_int_t i1,  magma_int_t i2,
    const magma_int_t *ipiv, magma_int_t inci );

void magmablas_zlaswpx(
    magma_int_t n,
    cuDoubleComplex *dAT, magma_int_t ldx, magma_int_t ldy,
    magma_int_t i1, magma_int_t i2,
    const magma_int_t *ipiv, magma_int_t inci );

void magmablas_zlaswp2(
    magma_int_t n,
    cuDoubleComplex* dAT, magma_int_t ldda,
    magma_int_t i1, magma_int_t i2,
    const magma_int_t *d_ipiv );

void magmablas_zsymmetrize(
    char uplo, magma_int_t m, cuDoubleComplex *dA, magma_int_t ldda );

void magmablas_zsymmetrize_tiles(
    char uplo, magma_int_t m, cuDoubleComplex *dA, magma_int_t ldda,
    magma_int_t ntile, magma_int_t mstride, magma_int_t nstride );

  /*
   * Level 1 BLAS
   */
void magmablas_zswap(
    magma_int_t n,
    cuDoubleComplex *dA, magma_int_t ldda,
    cuDoubleComplex *dB, magma_int_t lddb );

void magmablas_zswapblk(
    char storev,
    magma_int_t n,
    cuDoubleComplex *dA, magma_int_t ldda,
    cuDoubleComplex *dB, magma_int_t lddb,
    magma_int_t i1, magma_int_t i2,
    const magma_int_t *ipiv, magma_int_t inci,
    magma_int_t offset );

void magmablas_zswapdblk(
    magma_int_t n, magma_int_t nb,
    cuDoubleComplex *dA, magma_int_t ldda, magma_int_t inca,
    cuDoubleComplex *dB, magma_int_t lddb, magma_int_t incb );

  /*
   * Level 2 BLAS
   */
void magmablas_zgemv(
    char t, magma_int_t m, magma_int_t n,
    cuDoubleComplex alpha,
    const cuDoubleComplex *dA, magma_int_t ldda,
    const cuDoubleComplex *dx, magma_int_t incx,
    cuDoubleComplex beta,
    cuDoubleComplex       *dy, magma_int_t incy );

#if defined(PRECISION_z ) || defined(PRECISION_c )
magma_int_t magmablas_zhemv(
    char u, magma_int_t n,
    cuDoubleComplex alpha,
    const cuDoubleComplex *dA, magma_int_t ldda,
    const cuDoubleComplex *dx, magma_int_t incx,
    cuDoubleComplex beta,
    cuDoubleComplex       *dy, magma_int_t incy );
#endif

magma_int_t magmablas_zsymv(
    char u, magma_int_t n,
    cuDoubleComplex alpha,
    const cuDoubleComplex *dA, magma_int_t ldda,
    const cuDoubleComplex *dx, magma_int_t incx,
    cuDoubleComplex beta,
    cuDoubleComplex       *dy, magma_int_t incy );

  /*
   * Level 3 BLAS
   */
void magmablas_zgemm(
    char tA, char tB,
    magma_int_t m, magma_int_t n, magma_int_t k,
    cuDoubleComplex alpha,
    const cuDoubleComplex *dA, magma_int_t ldda,
    const cuDoubleComplex *dB, magma_int_t lddb,
    cuDoubleComplex beta,
    cuDoubleComplex       *dC, magma_int_t lddc );

void magmablas_zgemm_fermi80(
    char tA, char tB,
    magma_int_t m, magma_int_t n, magma_int_t k,
    cuDoubleComplex alpha,
    const cuDoubleComplex *dA, magma_int_t ldda,
    const cuDoubleComplex *dB, magma_int_t lddb,
    cuDoubleComplex beta,
    cuDoubleComplex       *dC, magma_int_t lddc );

void magmablas_zgemm_fermi64(
    char tA, char tB,
    magma_int_t m, magma_int_t n, magma_int_t k,
    cuDoubleComplex alpha,
    const cuDoubleComplex *dA, magma_int_t ldda,
    const cuDoubleComplex *dB, magma_int_t lddb,
    cuDoubleComplex beta,
    cuDoubleComplex       *dC, magma_int_t lddc );

void magmablas_zhemm(
    char s, char u,
    magma_int_t m, magma_int_t n,
    cuDoubleComplex alpha,
    const cuDoubleComplex *dA, magma_int_t ldda,
    const cuDoubleComplex *dB, magma_int_t lddb,
    cuDoubleComplex beta,
    cuDoubleComplex       *dC, magma_int_t lddc );

void magmablas_zsymm(
    char s, char u,
    magma_int_t m, magma_int_t n,
    cuDoubleComplex alpha,
    const cuDoubleComplex *dA, magma_int_t ldda,
    const cuDoubleComplex *dB, magma_int_t lddb,
    cuDoubleComplex beta,
    cuDoubleComplex       *dC, magma_int_t lddc );

void magmablas_zsyrk(
    char u, char t,
    magma_int_t n, magma_int_t k,
    cuDoubleComplex alpha,
    const cuDoubleComplex *dA, magma_int_t ldda,
    cuDoubleComplex beta,
    cuDoubleComplex       *dC, magma_int_t lddc );

void magmablas_zherk(
    char u, char t,
    magma_int_t n, magma_int_t k,
    double  alpha,
    const cuDoubleComplex *dA, magma_int_t ldda,
    double  beta,
    cuDoubleComplex       *dC, magma_int_t lddc );

void magmablas_zsyr2k(
    char u, char t,
    magma_int_t n, magma_int_t k,
    cuDoubleComplex alpha,
    const cuDoubleComplex *dA, magma_int_t ldda,
    const cuDoubleComplex *dB, magma_int_t lddb,
    cuDoubleComplex beta,
    cuDoubleComplex       *dC, magma_int_t lddc );

void magmablas_zher2k(
    char u, char t,
    magma_int_t n, magma_int_t k,
    cuDoubleComplex alpha,
    const cuDoubleComplex *dA, magma_int_t ldda,
    const cuDoubleComplex *dB, magma_int_t lddb,
    double  beta,
    cuDoubleComplex       *dC, magma_int_t lddc );

void magmablas_ztrmm(
    char s, char u, char t,  char d,
    magma_int_t m, magma_int_t n,
    cuDoubleComplex alpha,
    const cuDoubleComplex *dA, magma_int_t ldda,
    cuDoubleComplex       *dB, magma_int_t lddb );

void magmablas_ztrsm(
    char s, char u, char t, char d,
    magma_int_t m, magma_int_t n,
    cuDoubleComplex alpha,
    const cuDoubleComplex *dA, magma_int_t ldda,
    cuDoubleComplex       *dB, magma_int_t lddb );


  /*
   * Wrappers for platform independence.
   * These wrap CUBLAS or AMD OpenCL BLAS functions.
   */

// ========================================
// copying vectors
// set copies host to device
// get copies device to host
// Add the function, file, and line for error-reporting purposes.

#define magma_zsetvector(           n, hx_src, incx, dy_dst, incy ) \
        magma_zsetvector_internal(  n, hx_src, incx, dy_dst, incy, __func__, __FILE__, __LINE__ )

#define magma_zgetvector(           n, dx_src, incx, hy_dst, incy ) \
        magma_zgetvector_internal(  n, dx_src, incx, hy_dst, incy, __func__, __FILE__, __LINE__ )

#define magma_zsetvector_async(           n, hx_src, incx, dy_dst, incy, stream ) \
        magma_zsetvector_async_internal(  n, hx_src, incx, dy_dst, incy, stream, __func__, __FILE__, __LINE__ )

#define magma_zgetvector_async(           n, dx_src, incx, hy_dst, incy, stream ) \
        magma_zgetvector_async_internal(  n, dx_src, incx, hy_dst, incy, stream, __func__, __FILE__, __LINE__ )

void magma_zsetvector_internal(
    magma_int_t n,
    const cuDoubleComplex *hx_src, magma_int_t incx,
    cuDoubleComplex       *dy_dst, magma_int_t incy,
    const char* func, const char* file, int line );

void magma_zgetvector_internal(
    magma_int_t n,
    const cuDoubleComplex *dx_src, magma_int_t incx,
    cuDoubleComplex       *hy_dst, magma_int_t incy,
    const char* func, const char* file, int line );

void magma_zsetvector_async_internal(
    magma_int_t n,
    const cuDoubleComplex *hx_src, magma_int_t incx,
    cuDoubleComplex       *dy_dst, magma_int_t incy,
    magma_stream_t stream,
    const char* func, const char* file, int line );

void magma_zgetvector_async_internal(
    magma_int_t n,
    const cuDoubleComplex *dx_src, magma_int_t incx,
    cuDoubleComplex       *hy_dst, magma_int_t incy,
    magma_stream_t stream,
    const char* func, const char* file, int line );


// ========================================
// copying sub-matrices (contiguous columns )
// set  copies host to device
// get  copies device to host
// copy copies device to device
// (with CUDA unified addressing, copy can be between same or different devices )
// Add the function, file, and line for error-reporting purposes.

#define magma_zsetmatrix(           m, n, hA_src, lda, dB_dst, lddb ) \
        magma_zsetmatrix_internal(  m, n, hA_src, lda, dB_dst, lddb, __func__, __FILE__, __LINE__ )

#define magma_zgetmatrix(           m, n, dA_src, ldda, hB_dst, ldb ) \
        magma_zgetmatrix_internal(  m, n, dA_src, ldda, hB_dst, ldb, __func__, __FILE__, __LINE__ )

#define magma_zcopymatrix(          m, n, dA_src, ldda, dB_dst, lddb ) \
        magma_zcopymatrix_internal( m, n, dA_src, ldda, dB_dst, lddb, __func__, __FILE__, __LINE__ )

#define magma_zsetmatrix_async(           m, n, hA_src, lda, dB_dst, lddb, stream ) \
        magma_zsetmatrix_async_internal(  m, n, hA_src, lda, dB_dst, lddb, stream, __func__, __FILE__, __LINE__ )

#define magma_zgetmatrix_async(           m, n, dA_src, ldda, hB_dst, ldb, stream ) \
        magma_zgetmatrix_async_internal(  m, n, dA_src, ldda, hB_dst, ldb, stream, __func__, __FILE__, __LINE__ )

#define magma_zcopymatrix_async(          m, n, dA_src, ldda, dB_dst, lddb, stream ) \
        magma_zcopymatrix_async_internal( m, n, dA_src, ldda, dB_dst, lddb, stream, __func__, __FILE__, __LINE__ )

void magma_zsetmatrix_internal(
    magma_int_t m, magma_int_t n,
    const cuDoubleComplex *hA_src, magma_int_t lda,
    cuDoubleComplex       *dB_dst, magma_int_t lddb,
    const char* func, const char* file, int line );

void magma_zgetmatrix_internal(
    magma_int_t m, magma_int_t n,
    const cuDoubleComplex *dA_src, magma_int_t ldda,
    cuDoubleComplex       *hB_dst, magma_int_t ldb,
    const char* func, const char* file, int line );

void magma_zcopymatrix_internal(
    magma_int_t m, magma_int_t n,
    const cuDoubleComplex *dA_src, magma_int_t ldda,
    cuDoubleComplex       *dB_dst, magma_int_t lddb,
    const char* func, const char* file, int line );

void magma_zsetmatrix_async_internal(
    magma_int_t m, magma_int_t n,
    const cuDoubleComplex *hA_src, magma_int_t lda,
    cuDoubleComplex       *dB_dst, magma_int_t lddb,
    magma_stream_t stream,
    const char* func, const char* file, int line );

void magma_zgetmatrix_async_internal(
    magma_int_t m, magma_int_t n,
    const cuDoubleComplex *dA_src, magma_int_t ldda,
    cuDoubleComplex       *hB_dst, magma_int_t ldb,
    magma_stream_t stream,
    const char* func, const char* file, int line );

void magma_zcopymatrix_async_internal(
    magma_int_t m, magma_int_t n,
    const cuDoubleComplex *dA_src, magma_int_t ldda,
    cuDoubleComplex       *dB_dst, magma_int_t lddb,
    magma_stream_t stream,
    const char* func, const char* file, int line );


// ========================================
// Level 1 BLAS

void magma_zswap(
    magma_int_t n,
    cuDoubleComplex *dx, magma_int_t incx,
    cuDoubleComplex *dy, magma_int_t incy );

magma_int_t magma_izamax(
    magma_int_t n,
    cuDoubleComplex *dx, magma_int_t incx );

// ========================================
// Level 2 BLAS

void magma_zgemv(
    magma_trans_t transA,
    magma_int_t m, magma_int_t n,
    cuDoubleComplex alpha,
    const cuDoubleComplex *dA, magma_int_t ldda,
    const cuDoubleComplex *dx, magma_int_t incx,
    cuDoubleComplex beta,
    cuDoubleComplex       *dy, magma_int_t incy );

void magma_zhemv(
    magma_uplo_t uplo,
    magma_int_t n,
    cuDoubleComplex alpha,
    const cuDoubleComplex *dA, magma_int_t ldda,
    const cuDoubleComplex *dx, magma_int_t incx,
    cuDoubleComplex beta,
    cuDoubleComplex       *dy, magma_int_t incy );

void magma_ztrsv(
    magma_uplo_t uplo, magma_trans_t trans, magma_diag_t diag,
    magma_int_t n,
    const cuDoubleComplex *dA, magma_int_t ldda,
    cuDoubleComplex       *dx, magma_int_t incx );

// ========================================
// Level 3 BLAS

void magma_zgemm(
    magma_trans_t transA, magma_trans_t transB,
    magma_int_t m, magma_int_t n, magma_int_t k,
    cuDoubleComplex alpha,
    const cuDoubleComplex *dA, magma_int_t ldda,
    const cuDoubleComplex *dB, magma_int_t lddb,
    cuDoubleComplex beta,
    cuDoubleComplex       *dC, magma_int_t lddc );

void magma_zhemm(
    magma_side_t side, magma_uplo_t uplo,
    magma_int_t m, magma_int_t n,
    cuDoubleComplex alpha,
    const cuDoubleComplex *dA, magma_int_t ldda,
    const cuDoubleComplex *dB, magma_int_t lddb,
    cuDoubleComplex beta,
    cuDoubleComplex       *dC, magma_int_t lddc );

void magma_zherk(
    magma_uplo_t uplo, magma_trans_t trans,
    magma_int_t n, magma_int_t k,
    double alpha,
    const cuDoubleComplex *dA, magma_int_t ldda,
    double beta,
    cuDoubleComplex       *dC, magma_int_t lddc );

void magma_zher2k(
    magma_uplo_t uplo, magma_trans_t trans,
    magma_int_t n, magma_int_t k,
    cuDoubleComplex alpha,
    const cuDoubleComplex *dA, magma_int_t ldda,
    const cuDoubleComplex *dB, magma_int_t lddb,
    double beta,
    cuDoubleComplex       *dC, magma_int_t lddc );

void magma_ztrmm(
    magma_side_t side, magma_uplo_t uplo, magma_trans_t trans, magma_diag_t diag,
    magma_int_t m, magma_int_t n,
    cuDoubleComplex alpha,
    const cuDoubleComplex *dA, magma_int_t ldda,
    cuDoubleComplex       *dB, magma_int_t lddb );

void magma_ztrsm(
    magma_side_t side, magma_uplo_t uplo, magma_trans_t trans, magma_diag_t diag,
    magma_int_t m, magma_int_t n,
    cuDoubleComplex alpha,
    const cuDoubleComplex *dA, magma_int_t ldda,
    cuDoubleComplex       *dB, magma_int_t lddb );

#ifdef __cplusplus
}
#endif

#undef PRECISION_z

#endif  // MAGMABLAS_Z_H
