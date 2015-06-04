/*
 *   -- MAGMA (version 1.2.0) --
 *      Univ. of Tennessee, Knoxville
 *      Univ. of California, Berkeley
 *      Univ. of Colorado, Denver
 *      May 2012
 *
 * @generated d Tue May 15 18:17:05 2012
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
double cpu_gpu_ddiff(             int M, int N, 
                  double * a, int lda, 
                  double *da, int ldda);
void dzero_32x32_block(           double *, magma_int_t);
void dzero_nbxnb_block(           magma_int_t, double *, magma_int_t);
void magmablas_dinplace_transpose(double *, magma_int_t, magma_int_t);
void magmablas_dpermute_long(     double *, magma_int_t, 
                  magma_int_t *, magma_int_t, magma_int_t);
void magmablas_dpermute_long2(    double *, magma_int_t, 
                  magma_int_t *, magma_int_t, magma_int_t);
void magmablas_dpermute_long3( double *dAT, int lda, 
                               int *ipiv, int nb, int ind );
void magmablas_dtranspose(        double *, magma_int_t, 
                  double *, magma_int_t, 
                  magma_int_t, magma_int_t);
void magmablas_dtranspose2(       double *, magma_int_t, 
                  double *, magma_int_t, 
                  magma_int_t, magma_int_t);
void magmablas_dtranspose2s(double *odata, int ldo,
                       double *idata, int ldi,
                       int m, int n, cudaStream_t *stream );

void magmablas_dgetmatrix_transpose(  int m, int n,
                                      double *dat, int ldda,
                                      double  *ha, int lda,
                                      double  *dB, int lddb, int nb );
void magmablas_dgetmatrix_transpose2( int m, int n,
                                      double **dat, int *ldda,
                                      double  *ha,  int  lda,
                                      double **dB,  int  lddb, int nb,
                                      int num_gpus, cudaStream_t stream[][2] );
void magmablas_dsetmatrix_transpose(  int m, int n,
                                      double  *ha, int lda, 
                                      double *dat, int ldda,
                                      double  *dB, int lddb, int nb );
void magmablas_dsetmatrix_transpose2( int m, int n,
                                      double  *ha,  int  lda, 
                                      double **dat, int *ldda,
                                      double **dB,  int  lddb, int nb,
                                      int num_gpus, cudaStream_t stream[][2] );
void magmablas_dgetmatrix_1D_bcyclic( int m, int n,
                                      double  *da[], int ldda,
                                      double  *ha, int lda,
                                      int num_gpus, int nb );
void magmablas_dsetmatrix_1D_bcyclic( int m, int n,
                                      double  *ha, int lda,
                                      double  *da[], int ldda,
                                      int num_gpus, int nb );

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
void   magmablas_dlascl( char type, int kl, int ku,
             double cfrom, double cto,
             int m, int n,
             double *A, int lda, int *info );
void   magmablas_dlaset( char uplo, magma_int_t m, magma_int_t n,
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
