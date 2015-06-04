/*
 *   -- MAGMA (version 1.3.0) --
 *      Univ. of Tennessee, Knoxville
 *      Univ. of California, Berkeley
 *      Univ. of Colorado, Denver
 *      November 2012
 *
 * @generated d Wed Nov 14 22:52:26 2012
 */

#ifndef MAGMABLAS_D_H
#define MAGMABLAS_D_H

#define PRECISION_d

#ifdef __cplusplus
extern "C" {
#endif

  /*
   * Interface to clean
   */
double cpu_gpu_ddiff(
    magma_int_t m, magma_int_t n,
    const double *hA, magma_int_t lda,
    const double *dA, magma_int_t ldda );

// see also dlaset
void dzero_32x32_block(
    double *dA, magma_int_t ldda );

void dzero_nbxnb_block(
    magma_int_t nb, double *dA, magma_int_t ldda );

// see also dlaswp
void magmablas_dpermute_long2(
    magma_int_t n, double *dAT, magma_int_t ldda,
    magma_int_t *ipiv, magma_int_t nb, magma_int_t ind );

void magmablas_dpermute_long3(
    /*magma_int_t n,*/ double *dAT, magma_int_t ldda,
    const magma_int_t *ipiv, magma_int_t nb, magma_int_t ind );

  /*
   * Transpose functions
   */
void magmablas_dinplace_transpose(
    double *dA, magma_int_t ldda, magma_int_t n );

void magmablas_dtranspose(
    double       *odata, magma_int_t ldo,
    const double *idata, magma_int_t ldi,
    magma_int_t m, magma_int_t n );

void magmablas_dtranspose2(
    double       *odata, magma_int_t ldo,
    const double *idata, magma_int_t ldi,
    magma_int_t m, magma_int_t n );

void magmablas_dtranspose2s(
    double       *odata, magma_int_t ldo,
    const double *idata, magma_int_t ldi,
    magma_int_t m, magma_int_t n,
    cudaStream_t *stream );

void magmablas_dgetmatrix_transpose(
    magma_int_t m, magma_int_t n,
    const double *dAT,   magma_int_t ldda,
    double       *hA,    magma_int_t lda,
    double       *dwork, magma_int_t lddwork, magma_int_t nb );

void magmablas_dsetmatrix_transpose(
    magma_int_t m, magma_int_t n,
    const double *hA,    magma_int_t lda,
    double       *dAT,   magma_int_t ldda,
    double       *dwork, magma_int_t lddwork, magma_int_t nb );

  /*
   * Multi-GPU functions
   */
void magmablas_dgetmatrix_transpose_mgpu(
    magma_int_t ngpu, cudaStream_t stream[][2],
    double **dAT, magma_int_t ldda,
    double  *hA,  magma_int_t lda,
    double **dB,  magma_int_t lddb,
    magma_int_t m, magma_int_t n, magma_int_t nb );

void magmablas_dsetmatrix_transpose_mgpu(
    magma_int_t ngpu, cudaStream_t stream[][2],
    const double  *hA,  magma_int_t lda,
    double       **dAT, magma_int_t ldda,
    double       **dB,  magma_int_t lddb,
    magma_int_t m, magma_int_t n, magma_int_t nb );

void magmablas_dgetmatrix_1D_bcyclic(
    magma_int_t m, magma_int_t n,
    double *dA[], magma_int_t ldda,
    double       *hA,   magma_int_t lda,
    magma_int_t ngpu, magma_int_t nb );

void magmablas_dsetmatrix_1D_bcyclic(
    magma_int_t m, magma_int_t n,
    const double *hA,   magma_int_t lda,
    double       *dA[], magma_int_t ldda,
    magma_int_t ngpu, magma_int_t nb );

void magmablas_dsymm_1gpu_old(
    char side, char uplo, magma_int_t m, magma_int_t n,
    double alpha,
    double *dA[], magma_int_t ldda,  magma_int_t offset,
    double *dB[], magma_int_t lddb,
    double beta,
    double *dC[], magma_int_t lddc,
    double *C,    magma_int_t ldc,
    magma_int_t ngpu, magma_int_t nb,
    cudaStream_t streams[][20], magma_int_t nstream );

void magmablas_dsymm_1gpu(
    char side, char uplo, magma_int_t m, magma_int_t n,
    double alpha,
    double *dA[], magma_int_t ldda,  magma_int_t offset,
    double *dB[], magma_int_t lddb,
    double beta,
    double *dC[], magma_int_t lddc,
    double *C,    magma_int_t ldc,
    magma_int_t ngpu, magma_int_t nb,
    cudaStream_t streams[][20], magma_int_t nstream );


void magmablas_dsymm_mgpu(
    char side, char uplo, magma_int_t m, magma_int_t n,
    double alpha, double *dA[], magma_int_t ldda,  magma_int_t offset,
                           double *dB[], magma_int_t lddb,
    double beta,  double *dC[], magma_int_t lddc,
                           double *dwork[],    magma_int_t lddwork,
                           double *C,    magma_int_t ldc,
                           double *work[], magma_int_t ldwork,
                           magma_int_t ngpu, magma_int_t nb, 
                           cudaStream_t streams[][20], magma_int_t nstream, 
                           cudaEvent_t redevents[][20],magma_int_t nbevents );

void magmablas_dsymm_mgpu_com(
    char side, char uplo, magma_int_t m, magma_int_t n,
    double alpha, double *dA[], magma_int_t ldda,  magma_int_t offset,
                           double *dB[], magma_int_t lddb,
    double beta,  double *dC[], magma_int_t lddc,
                           double *dwork[],    magma_int_t lddwork,
                           double *C,    magma_int_t ldc,
                           double *work[], magma_int_t ldwork,
                           magma_int_t ngpu, magma_int_t nb, 
                           cudaStream_t streams[][20], magma_int_t nstream, 
                           cudaEvent_t redevents[][20],magma_int_t nbevents, 
                           magma_int_t gnode[MagmaMaxGPUs][MagmaMaxGPUs+2], magma_int_t nbcmplx );

void magmablas_dsymm_mgpu_spec(
    char side, char uplo, magma_int_t m, magma_int_t n,
    double alpha, double *dA[], magma_int_t ldda,  magma_int_t offset,
                           double *dB[], magma_int_t lddb,
    double beta,  double *dC[], magma_int_t lddc,
                           double *dwork[],    magma_int_t lddwork,
                           double *C,    magma_int_t ldc,
                           double *work[], magma_int_t ldwork,
                           magma_int_t ngpu, magma_int_t nb, 
                           cudaStream_t streams[][20], magma_int_t nstream, 
                           cudaEvent_t redevents[][MagmaMaxGPUs*MagmaMaxGPUs],magma_int_t nbevents, 
                           magma_int_t gnode[MagmaMaxGPUs][MagmaMaxGPUs+2], magma_int_t nbcmplx );



void magmablas_dsyr2k_mgpu2(
    char uplo, char trans, magma_int_t n, magma_int_t k,
    double alpha,
    double *dA[], magma_int_t ldda, magma_int_t aoff,
    double *dB[], magma_int_t lddb, magma_int_t boff,
    double beta,
    double *dC[], magma_int_t lddc, magma_int_t offset,
    magma_int_t ngpu, magma_int_t nb,
    cudaStream_t streams[][20], magma_int_t nstream );

void magmablas_dsyr2k_mgpu_spec(
    char uplo, char trans, magma_int_t n, magma_int_t k,
    double alpha, double *dA[], magma_int_t lda, magma_int_t aoff,
                           double *dB[], magma_int_t ldb, magma_int_t boff,
    double beta,           double *dC[], magma_int_t ldc,  magma_int_t offset,
    magma_int_t ngpu, magma_int_t nb, cudaStream_t streams[][20], magma_int_t nstream );

  /*
   * LAPACK auxiliary functions
   */
void magmablas_dgeadd(
    magma_int_t m, magma_int_t n,
    double alpha,
    const double *dA, magma_int_t ldda,
    double       *dB, magma_int_t lddb );

void magmablas_dlacpy(
    char uplo,
    magma_int_t m, magma_int_t n,
    const double *dA, magma_int_t ldda,
    double       *dB, magma_int_t lddb );

double magmablas_dlange(
    char norm,
    magma_int_t m, magma_int_t n,
    const double *dA, magma_int_t ldda, double *dwork );

double magmablas_dlansy(
    char norm, char uplo,
    magma_int_t n,
    const double *dA, magma_int_t ldda, double *dwork );

double magmablas_dlansy(
    char norm, char uplo,
    magma_int_t n,
    const double *dA, magma_int_t ldda, double *dwork );

void magmablas_dlascl(
    char type, magma_int_t kl, magma_int_t ku,
    double cfrom, double cto,
    magma_int_t m, magma_int_t n,
    double *dA, magma_int_t ldda, magma_int_t *info );

void magmablas_dlaset(
    char uplo, magma_int_t m, magma_int_t n,
    double *dA, magma_int_t ldda );

void magmablas_dlaset_identity(
    magma_int_t m, magma_int_t n,
    double *dA, magma_int_t ldda );

void magmablas_dlaswp(
    magma_int_t n,
    double *dAT, magma_int_t ldda,
    magma_int_t i1,  magma_int_t i2,
    const magma_int_t *ipiv, magma_int_t inci );

void magmablas_dlaswpx(
    magma_int_t n,
    double *dAT, magma_int_t ldx, magma_int_t ldy,
    magma_int_t i1, magma_int_t i2,
    const magma_int_t *ipiv, magma_int_t inci );

void magmablas_dlaswp2(
    magma_int_t n,
    double* dAT, magma_int_t ldda,
    magma_int_t i1, magma_int_t i2,
    const magma_int_t *d_ipiv );

void magmablas_dsymmetrize(
    char uplo, magma_int_t m, double *dA, magma_int_t ldda );

void magmablas_dsymmetrize_tiles(
    char uplo, magma_int_t m, double *dA, magma_int_t ldda,
    magma_int_t ntile, magma_int_t mstride, magma_int_t nstride );

  /*
   * Level 1 BLAS
   */
void magmablas_dswap(
    magma_int_t n,
    double *dA, magma_int_t ldda,
    double *dB, magma_int_t lddb );

void magmablas_dswapblk(
    char storev,
    magma_int_t n,
    double *dA, magma_int_t ldda,
    double *dB, magma_int_t lddb,
    magma_int_t i1, magma_int_t i2,
    const magma_int_t *ipiv, magma_int_t inci,
    magma_int_t offset );

void magmablas_dswapdblk(
    magma_int_t n, magma_int_t nb,
    double *dA, magma_int_t ldda, magma_int_t inca,
    double *dB, magma_int_t lddb, magma_int_t incb );

  /*
   * Level 2 BLAS
   */
void magmablas_dgemv(
    char t, magma_int_t m, magma_int_t n,
    double alpha,
    const double *dA, magma_int_t ldda,
    const double *dx, magma_int_t incx,
    double beta,
    double       *dy, magma_int_t incy );

#if defined(PRECISION_z ) || defined(PRECISION_c )
magma_int_t magmablas_dsymv(
    char u, magma_int_t n,
    double alpha,
    const double *dA, magma_int_t ldda,
    const double *dx, magma_int_t incx,
    double beta,
    double       *dy, magma_int_t incy );
#endif

magma_int_t magmablas_dsymv(
    char u, magma_int_t n,
    double alpha,
    const double *dA, magma_int_t ldda,
    const double *dx, magma_int_t incx,
    double beta,
    double       *dy, magma_int_t incy );

  /*
   * Level 3 BLAS
   */
void magmablas_dgemm(
    char tA, char tB,
    magma_int_t m, magma_int_t n, magma_int_t k,
    double alpha,
    const double *dA, magma_int_t ldda,
    const double *dB, magma_int_t lddb,
    double beta,
    double       *dC, magma_int_t lddc );

void magmablas_dgemm_fermi80(
    char tA, char tB,
    magma_int_t m, magma_int_t n, magma_int_t k,
    double alpha,
    const double *dA, magma_int_t ldda,
    const double *dB, magma_int_t lddb,
    double beta,
    double       *dC, magma_int_t lddc );

void magmablas_dgemm_fermi64(
    char tA, char tB,
    magma_int_t m, magma_int_t n, magma_int_t k,
    double alpha,
    const double *dA, magma_int_t ldda,
    const double *dB, magma_int_t lddb,
    double beta,
    double       *dC, magma_int_t lddc );

void magmablas_dsymm(
    char s, char u,
    magma_int_t m, magma_int_t n,
    double alpha,
    const double *dA, magma_int_t ldda,
    const double *dB, magma_int_t lddb,
    double beta,
    double       *dC, magma_int_t lddc );

void magmablas_dsymm(
    char s, char u,
    magma_int_t m, magma_int_t n,
    double alpha,
    const double *dA, magma_int_t ldda,
    const double *dB, magma_int_t lddb,
    double beta,
    double       *dC, magma_int_t lddc );

void magmablas_dsyrk(
    char u, char t,
    magma_int_t n, magma_int_t k,
    double alpha,
    const double *dA, magma_int_t ldda,
    double beta,
    double       *dC, magma_int_t lddc );

void magmablas_dsyrk(
    char u, char t,
    magma_int_t n, magma_int_t k,
    double  alpha,
    const double *dA, magma_int_t ldda,
    double  beta,
    double       *dC, magma_int_t lddc );

void magmablas_dsyr2k(
    char u, char t,
    magma_int_t n, magma_int_t k,
    double alpha,
    const double *dA, magma_int_t ldda,
    const double *dB, magma_int_t lddb,
    double beta,
    double       *dC, magma_int_t lddc );

void magmablas_dsyr2k(
    char u, char t,
    magma_int_t n, magma_int_t k,
    double alpha,
    const double *dA, magma_int_t ldda,
    const double *dB, magma_int_t lddb,
    double  beta,
    double       *dC, magma_int_t lddc );

void magmablas_dtrmm(
    char s, char u, char t,  char d,
    magma_int_t m, magma_int_t n,
    double alpha,
    const double *dA, magma_int_t ldda,
    double       *dB, magma_int_t lddb );

void magmablas_dtrsm(
    char s, char u, char t, char d,
    magma_int_t m, magma_int_t n,
    double alpha,
    const double *dA, magma_int_t ldda,
    double       *dB, magma_int_t lddb );


  /*
   * Wrappers for platform independence.
   * These wrap CUBLAS or AMD OpenCL BLAS functions.
   */

// ========================================
// copying vectors
// set copies host to device
// get copies device to host
// Add the function, file, and line for error-reporting purposes.

#define magma_dsetvector(           n, hx_src, incx, dy_dst, incy ) \
        magma_dsetvector_internal(  n, hx_src, incx, dy_dst, incy, __func__, __FILE__, __LINE__ )

#define magma_dgetvector(           n, dx_src, incx, hy_dst, incy ) \
        magma_dgetvector_internal(  n, dx_src, incx, hy_dst, incy, __func__, __FILE__, __LINE__ )

#define magma_dsetvector_async(           n, hx_src, incx, dy_dst, incy, stream ) \
        magma_dsetvector_async_internal(  n, hx_src, incx, dy_dst, incy, stream, __func__, __FILE__, __LINE__ )

#define magma_dgetvector_async(           n, dx_src, incx, hy_dst, incy, stream ) \
        magma_dgetvector_async_internal(  n, dx_src, incx, hy_dst, incy, stream, __func__, __FILE__, __LINE__ )

void magma_dsetvector_internal(
    magma_int_t n,
    const double *hx_src, magma_int_t incx,
    double       *dy_dst, magma_int_t incy,
    const char* func, const char* file, int line );

void magma_dgetvector_internal(
    magma_int_t n,
    const double *dx_src, magma_int_t incx,
    double       *hy_dst, magma_int_t incy,
    const char* func, const char* file, int line );

void magma_dsetvector_async_internal(
    magma_int_t n,
    const double *hx_src, magma_int_t incx,
    double       *dy_dst, magma_int_t incy,
    magma_stream_t stream,
    const char* func, const char* file, int line );

void magma_dgetvector_async_internal(
    magma_int_t n,
    const double *dx_src, magma_int_t incx,
    double       *hy_dst, magma_int_t incy,
    magma_stream_t stream,
    const char* func, const char* file, int line );


// ========================================
// copying sub-matrices (contiguous columns )
// set  copies host to device
// get  copies device to host
// copy copies device to device
// (with CUDA unified addressing, copy can be between same or different devices )
// Add the function, file, and line for error-reporting purposes.

#define magma_dsetmatrix(           m, n, hA_src, lda, dB_dst, lddb ) \
        magma_dsetmatrix_internal(  m, n, hA_src, lda, dB_dst, lddb, __func__, __FILE__, __LINE__ )

#define magma_dgetmatrix(           m, n, dA_src, ldda, hB_dst, ldb ) \
        magma_dgetmatrix_internal(  m, n, dA_src, ldda, hB_dst, ldb, __func__, __FILE__, __LINE__ )

#define magma_dcopymatrix(          m, n, dA_src, ldda, dB_dst, lddb ) \
        magma_dcopymatrix_internal( m, n, dA_src, ldda, dB_dst, lddb, __func__, __FILE__, __LINE__ )

#define magma_dsetmatrix_async(           m, n, hA_src, lda, dB_dst, lddb, stream ) \
        magma_dsetmatrix_async_internal(  m, n, hA_src, lda, dB_dst, lddb, stream, __func__, __FILE__, __LINE__ )

#define magma_dgetmatrix_async(           m, n, dA_src, ldda, hB_dst, ldb, stream ) \
        magma_dgetmatrix_async_internal(  m, n, dA_src, ldda, hB_dst, ldb, stream, __func__, __FILE__, __LINE__ )

#define magma_dcopymatrix_async(          m, n, dA_src, ldda, dB_dst, lddb, stream ) \
        magma_dcopymatrix_async_internal( m, n, dA_src, ldda, dB_dst, lddb, stream, __func__, __FILE__, __LINE__ )

void magma_dsetmatrix_internal(
    magma_int_t m, magma_int_t n,
    const double *hA_src, magma_int_t lda,
    double       *dB_dst, magma_int_t lddb,
    const char* func, const char* file, int line );

void magma_dgetmatrix_internal(
    magma_int_t m, magma_int_t n,
    const double *dA_src, magma_int_t ldda,
    double       *hB_dst, magma_int_t ldb,
    const char* func, const char* file, int line );

void magma_dcopymatrix_internal(
    magma_int_t m, magma_int_t n,
    const double *dA_src, magma_int_t ldda,
    double       *dB_dst, magma_int_t lddb,
    const char* func, const char* file, int line );

void magma_dsetmatrix_async_internal(
    magma_int_t m, magma_int_t n,
    const double *hA_src, magma_int_t lda,
    double       *dB_dst, magma_int_t lddb,
    magma_stream_t stream,
    const char* func, const char* file, int line );

void magma_dgetmatrix_async_internal(
    magma_int_t m, magma_int_t n,
    const double *dA_src, magma_int_t ldda,
    double       *hB_dst, magma_int_t ldb,
    magma_stream_t stream,
    const char* func, const char* file, int line );

void magma_dcopymatrix_async_internal(
    magma_int_t m, magma_int_t n,
    const double *dA_src, magma_int_t ldda,
    double       *dB_dst, magma_int_t lddb,
    magma_stream_t stream,
    const char* func, const char* file, int line );


// ========================================
// Level 1 BLAS

void magma_dswap(
    magma_int_t n,
    double *dx, magma_int_t incx,
    double *dy, magma_int_t incy );

magma_int_t magma_idamax(
    magma_int_t n,
    double *dx, magma_int_t incx );

// ========================================
// Level 2 BLAS

void magma_dgemv(
    magma_trans_t transA,
    magma_int_t m, magma_int_t n,
    double alpha,
    const double *dA, magma_int_t ldda,
    const double *dx, magma_int_t incx,
    double beta,
    double       *dy, magma_int_t incy );

void magma_dsymv(
    magma_uplo_t uplo,
    magma_int_t n,
    double alpha,
    const double *dA, magma_int_t ldda,
    const double *dx, magma_int_t incx,
    double beta,
    double       *dy, magma_int_t incy );

void magma_dtrsv(
    magma_uplo_t uplo, magma_trans_t trans, magma_diag_t diag,
    magma_int_t n,
    const double *dA, magma_int_t ldda,
    double       *dx, magma_int_t incx );

// ========================================
// Level 3 BLAS

void magma_dgemm(
    magma_trans_t transA, magma_trans_t transB,
    magma_int_t m, magma_int_t n, magma_int_t k,
    double alpha,
    const double *dA, magma_int_t ldda,
    const double *dB, magma_int_t lddb,
    double beta,
    double       *dC, magma_int_t lddc );

void magma_dsymm(
    magma_side_t side, magma_uplo_t uplo,
    magma_int_t m, magma_int_t n,
    double alpha,
    const double *dA, magma_int_t ldda,
    const double *dB, magma_int_t lddb,
    double beta,
    double       *dC, magma_int_t lddc );

void magma_dsyrk(
    magma_uplo_t uplo, magma_trans_t trans,
    magma_int_t n, magma_int_t k,
    double alpha,
    const double *dA, magma_int_t ldda,
    double beta,
    double       *dC, magma_int_t lddc );

void magma_dsyr2k(
    magma_uplo_t uplo, magma_trans_t trans,
    magma_int_t n, magma_int_t k,
    double alpha,
    const double *dA, magma_int_t ldda,
    const double *dB, magma_int_t lddb,
    double beta,
    double       *dC, magma_int_t lddc );

void magma_dtrmm(
    magma_side_t side, magma_uplo_t uplo, magma_trans_t trans, magma_diag_t diag,
    magma_int_t m, magma_int_t n,
    double alpha,
    const double *dA, magma_int_t ldda,
    double       *dB, magma_int_t lddb );

void magma_dtrsm(
    magma_side_t side, magma_uplo_t uplo, magma_trans_t trans, magma_diag_t diag,
    magma_int_t m, magma_int_t n,
    double alpha,
    const double *dA, magma_int_t ldda,
    double       *dB, magma_int_t lddb );

#ifdef __cplusplus
}
#endif

#undef PRECISION_d

#endif  // MAGMABLAS_D_H
