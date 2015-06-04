/*
 *   -- MAGMA (version 1.2.1) --
 *      Univ. of Tennessee, Knoxville
 *      Univ. of California, Berkeley
 *      Univ. of Colorado, Denver
 *      June 2012
 *
 * @generated d Thu Jun 28 12:30:02 2012
 */

#ifndef _MAGMABLAS_D_H_
#define _MAGMABLAS_D_H_

#define PRECISION_d

#ifdef __cplusplus
extern "C" {
#endif

  /*
   * Interface to clean
   */
double cpu_gpu_ddiff(             magma_int_t M, magma_int_t N, 
                  double * a, magma_int_t lda, 
                  double *da, magma_int_t ldda);

void dzero_32x32_block(           double *, magma_int_t);

void dzero_nbxnb_block(           magma_int_t, double *, magma_int_t);

void magmablas_dpermute_long(     double *, magma_int_t, 
                  magma_int_t *, magma_int_t, magma_int_t);

void magmablas_dpermute_long2(magma_int_t n, double *, magma_int_t, 
                  magma_int_t *, magma_int_t, magma_int_t);

void magmablas_dpermute_long3( double *dAT, magma_int_t lda, 
                               magma_int_t *ipiv, magma_int_t nb, magma_int_t ind );

  /*
   * Transpose functions
   */
void magmablas_dinplace_transpose(double *, magma_int_t, magma_int_t);

void magmablas_dtranspose(        double *, magma_int_t, 
                  double *, magma_int_t, 
                  magma_int_t, magma_int_t);

void magmablas_dtranspose2(       double *, magma_int_t, 
                  double *, magma_int_t, 
                  magma_int_t, magma_int_t);

void magmablas_dtranspose2s(double *odata, magma_int_t ldo,
                       double *idata, magma_int_t ldi,
                       magma_int_t m, magma_int_t n, cudaStream_t *stream );

void magmablas_dgetmatrix_transpose(  magma_int_t m, magma_int_t n,
                                      double *dat, magma_int_t ldda,
                                      double  *ha, magma_int_t lda,
                                      double  *dB, magma_int_t lddb, magma_int_t nb );
void magmablas_dsetmatrix_transpose(  magma_int_t m, magma_int_t n,
                                      double  *ha, magma_int_t lda, 
                                      double *dat, magma_int_t ldda,
                                      double  *dB, magma_int_t lddb, magma_int_t nb );

  /*
   * Multi-GPU functions
   */
void magmablas_dgetmatrix_transpose_mgpu(
                  magma_int_t num_gpus, cudaStream_t **stream0,
                  double **dat, magma_int_t ldda,
                  double   *ha, magma_int_t lda,
                  double  **dB, magma_int_t lddb,
                  magma_int_t m, magma_int_t n, magma_int_t nb);

void magmablas_dsetmatrix_transpose_mgpu(
                  magma_int_t num_gpus, cudaStream_t **stream0,
                  double  *ha,  magma_int_t lda, 
                  double **dat, magma_int_t ldda, magma_int_t starti,
                  double **dB,  magma_int_t lddb,
                  magma_int_t m, magma_int_t n, magma_int_t nb);

void magmablas_dgetmatrix_1D_bcyclic( magma_int_t m, magma_int_t n,
                                      double  *da[], magma_int_t ldda,
                                      double  *ha, magma_int_t lda,
                                      magma_int_t num_gpus, magma_int_t nb );

void magmablas_dsetmatrix_1D_bcyclic( magma_int_t m, magma_int_t n,
                                      double  *ha, magma_int_t lda,
                                      double  *da[], magma_int_t ldda,
                                      magma_int_t num_gpus, magma_int_t nb );

void magmablas_dsymm_mgpu(
    char side, char uplo, magma_int_t m, magma_int_t n,
    double alpha, double *dA[], magma_int_t ldda,  magma_int_t offset,
                           double *dB[], magma_int_t lddb,
    double beta,  double *dC[], magma_int_t lddc,
                           double *C,    magma_int_t ldc,
    magma_int_t ngpu, magma_int_t nb, cudaStream_t streams[][20], magma_int_t nstream );

void magmablas_dsyr2k_mgpu2(
    char uplo, char trans, magma_int_t n, magma_int_t k,
    double alpha, double *dA[], magma_int_t lda,
                           double *dB[], magma_int_t ldb,
    double beta,           double *dC[], magma_int_t ldc,  magma_int_t offset,
    magma_int_t ngpu, magma_int_t nb, cudaStream_t streams[][10], magma_int_t nstream );

  /*
   * LAPACK auxiliary functions
   */
void   magmablas_dlacpy( char uplo, 
             magma_int_t m, magma_int_t n, 
             double *A, magma_int_t lda, 
             double *B, magma_int_t ldb);

double magmablas_dlange( char norm, 
             magma_int_t m, magma_int_t n, 
             double *A, magma_int_t lda, double *WORK);

double magmablas_dlansy( char norm, char uplo, 
             magma_int_t n,
             double *A, magma_int_t lda, double *WORK);

double magmablas_dlansy( char norm, char uplo,
             magma_int_t n, 
             double *A, magma_int_t lda, double *WORK);

void   magmablas_dlascl( char type, magma_int_t kl, magma_int_t ku,
             double cfrom, double cto,
             magma_int_t m, magma_int_t n,
             double *A, magma_int_t lda, magma_int_t *info );

void   magmablas_dlaset( char uplo, magma_int_t m, magma_int_t n,
             double *A, magma_int_t lda);

void   magmablas_dlaset_identity(
             magma_int_t m, magma_int_t n,
             double *A, magma_int_t lda);

void   magmablas_dlaswp( magma_int_t N, 
             double *dAT, magma_int_t lda, 
             magma_int_t i1,  magma_int_t i2, 
             magma_int_t *ipiv, magma_int_t inci );

void   magmablas_dlaswpx(magma_int_t N, 
             double *dAT, magma_int_t ldx, magma_int_t ldy, 
             magma_int_t i1, magma_int_t i2,
             magma_int_t *ipiv, magma_int_t inci );

  /*
   * Level 1 BLAS
   */
void   magmablas_dswap(   magma_int_t N, 
              double *dA1, magma_int_t lda1, 
              double *dA2, magma_int_t lda2 );

void   magmablas_dswapblk(char storev, 
              magma_int_t N, 
              double *dA1, magma_int_t lda1, 
              double *dA2, magma_int_t lda2,
              magma_int_t i1, magma_int_t i2, 
              magma_int_t *ipiv, magma_int_t inci, 
              magma_int_t offset);

void magmablas_dswapdblk(magma_int_t n, magma_int_t nb,
             double *dA1, magma_int_t ldda1, magma_int_t inca1,
             double *dA2, magma_int_t ldda2, magma_int_t inca2 );

  /*
   * Level 2 BLAS
   */
void magmablas_dgemv(char t, magma_int_t M, magma_int_t N, 
             double alpha,
             double *A, magma_int_t lda, 
             double * X, magma_int_t incX, 
             double beta, 
             double *Y, magma_int_t incY);

#if defined(PRECISION_z) || defined(PRECISION_c)
magma_int_t magmablas_dsymv(char u, magma_int_t N, 
                            double alpha, 
                            double *A, magma_int_t lda, 
                            double *X, magma_int_t incX, 
                            double beta, 
                            double *Y, magma_int_t incY);
#endif

magma_int_t magmablas_dsymv(char u, magma_int_t N, 
                            double alpha, 
                            double *A, magma_int_t lda, 
                            double *X, magma_int_t incX, 
                            double beta, 
                            double *Y, magma_int_t incY);

  /*
   * Level 3 BLAS
   */
void magmablas_dgemm(char tA, char tB,
             magma_int_t m, magma_int_t n, magma_int_t k, 
             double alpha,
             const double *A, magma_int_t lda, 
             const double *B, magma_int_t ldb, 
             double beta,
             double *C, magma_int_t ldc);

void magmablas_dgemm_fermi80(char tA, char tB, 
                 magma_int_t m, magma_int_t n, magma_int_t k,
                 double alpha, 
                 const double *A, magma_int_t lda, 
                 const double *B, magma_int_t ldb,
                 double beta, 
                 double *C, magma_int_t ldc);

void magmablas_dgemm_fermi64(char tA, char tB, 
                 magma_int_t m, magma_int_t n, magma_int_t k,
                 double alpha, 
                 const double *A, magma_int_t lda, 
                 const double *B, magma_int_t ldb, 
                 double beta, 
                 double *C, magma_int_t ldc);

void magmablas_dsymm(char s, char u,          
             magma_int_t m, magma_int_t n,
             double alpha, 
             const double *A, magma_int_t lda,
             const double *B, magma_int_t ldb,
             double beta, 
             double *C, magma_int_t ldc);

void magmablas_dsymm(char s, char u,
             magma_int_t m, magma_int_t n,
             double alpha, 
             const double *A, magma_int_t lda, 
             const double *B, magma_int_t ldb,
             double beta,
             double *C, magma_int_t ldc);

void magmablas_dsyrk(char u, char t,
             magma_int_t n, magma_int_t k, 
             double alpha, 
             const double *A, magma_int_t lda,
             double beta,
             double *C, magma_int_t ldc);

void magmablas_dsyrk(char u, char t,
             magma_int_t n, magma_int_t k, 
             double  alpha, 
             const double *A, magma_int_t lda,
             double  beta, 
             double *C, magma_int_t ldc);

void magmablas_dsyr2k(char u, char t,
              magma_int_t n, magma_int_t k,
              double alpha, 
              const double *A, magma_int_t lda,
              const double *B, magma_int_t ldb, 
              double beta, 
              double *C, magma_int_t ldc);

void magmablas_dsyr2k(char u, char t,
              magma_int_t n, magma_int_t k, 
              double alpha, 
              const double *A, magma_int_t lda, 
              const double *B, magma_int_t ldb,
              double  beta,
              double *C, magma_int_t ldc);

void magmablas_dtrmm(char s, char u, char t,  char d, 
             magma_int_t m, magma_int_t n,
             double alpha,
             const double *A, magma_int_t lda,
             double *B, magma_int_t ldb);

void magmablas_dtrsm(char s, char u, char t, char d,
             magma_int_t m, magma_int_t n,
             double alpha,
             /*const*/ double *A, magma_int_t lda,
             double *B, magma_int_t ldb);


  /*
   * Wrappers for platform independence.
   * These wrap CUBLAS or AMD OpenCL BLAS functions.
   */

// ========================================
// copying vectors
// set copies host to device
// get copies device to host

void magma_dsetvector(
    magma_int_t n,
    double const *hx_src, magma_int_t incx,
    double       *dy_dst, magma_int_t incy );

void magma_dgetvector(
    magma_int_t n,
    double const *dx_src, magma_int_t incx,
    double       *hy_dst, magma_int_t incy );

void magma_dsetvector_async(
    magma_int_t n,
    double const *hx_src, magma_int_t incx,
    double       *dy_dst, magma_int_t incy,
    magma_stream_t stream );

void magma_dgetvector_async(
    magma_int_t n,
    double const *dx_src, magma_int_t incx,
    double       *hy_dst, magma_int_t incy,
    magma_stream_t stream );


// ========================================
// copying sub-matrices (contiguous columns)
// set copies host to device
// get copies device to host
// cpy copies device to device (with CUDA unified addressing, can be same or different devices)

void magma_dsetmatrix(
    magma_int_t m, magma_int_t n,
    double const *hA_src, magma_int_t lda,
    double       *dB_dst, magma_int_t ldb );

void magma_dgetmatrix(
    magma_int_t m, magma_int_t n,
    double const *dA_src, magma_int_t lda,
    double       *hB_dst, magma_int_t ldb );

void magma_dsetmatrix_async(
    magma_int_t m, magma_int_t n,
    double const *hA_src, magma_int_t lda,
    double       *dB_dst, magma_int_t ldb,
    magma_stream_t stream );

void magma_dgetmatrix_async(
    magma_int_t m, magma_int_t n,
    double const *dA_src, magma_int_t lda,
    double       *hB_dst, magma_int_t ldb,
    magma_stream_t stream );

void magma_dcopymatrix(
    magma_int_t m, magma_int_t n,
    double const *dA_src, magma_int_t lda,
    double       *dB_dst, magma_int_t ldb );

void magma_dcopymatrix_async(
    magma_int_t m, magma_int_t n,
    double const *dA_src, magma_int_t lda,
    double       *dB_dst, magma_int_t ldb,
    magma_stream_t stream );


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
    double alpha, double const *dA, magma_int_t lda,
                           double const *dx, magma_int_t incx,
    double beta,  double       *dy, magma_int_t incy );

void magma_dsymv(
    magma_uplo_t uplo,
    magma_int_t n,
    double alpha, double const *dA, magma_int_t lda,
                           double const *dx, magma_int_t incx,
    double beta,  double       *dy, magma_int_t incy );

void magma_dtrsv(
    magma_uplo_t uplo, magma_trans_t trans, magma_diag_t diag, 
    magma_int_t n, 
    double const *dA, magma_int_t lda, 
    double       *dx, magma_int_t incx );

// ========================================
// Level 3 BLAS

void magma_dgemm(
    magma_trans_t transA, magma_trans_t transB,
    magma_int_t m, magma_int_t n, magma_int_t k,
    double alpha, double const *dA, magma_int_t lda,
                           double const *dB, magma_int_t ldb,
    double beta,  double       *dC, magma_int_t ldc );

void magma_dsymm(
    magma_side_t side, magma_uplo_t uplo,
    magma_int_t m, magma_int_t n,
    double alpha, double const *dA, magma_int_t lda,
                           double const *dB, magma_int_t ldb,
    double beta,  double       *dC, magma_int_t ldc );

void magma_dsyrk(
    magma_uplo_t uplo, magma_trans_t trans,
    magma_int_t n, magma_int_t k,
    double alpha, double const *dA, magma_int_t lda,
    double beta,  double       *dC, magma_int_t ldc );

void magma_dsyr2k(
    magma_uplo_t uplo, magma_trans_t trans,
    magma_int_t n, magma_int_t k,
    double alpha, double const *dA, magma_int_t lda,
                           double const *dB, magma_int_t ldb,
    double beta,           double       *dC, magma_int_t ldc );

void magma_dtrmm(
    magma_side_t side, magma_uplo_t uplo, magma_trans_t trans, magma_diag_t diag,
    magma_int_t m, magma_int_t n,
    double alpha, double const *dA, magma_int_t lda,
                           double       *dB, magma_int_t ldb );

void magma_dtrsm(
    magma_side_t side, magma_uplo_t uplo, magma_trans_t trans, magma_diag_t diag,
    magma_int_t m, magma_int_t n,
    double alpha, double const *dA, magma_int_t lda,
                           double       *dB, magma_int_t ldb );

#ifdef __cplusplus
}
#endif

#undef PRECISION_d
#endif
