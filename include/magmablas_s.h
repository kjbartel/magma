/*
 *   -- MAGMA (version 1.1) --
 *      Univ. of Tennessee, Knoxville
 *      Univ. of California, Berkeley
 *      Univ. of Colorado, Denver
 *      November 2011
 *
 * @generated s Sun Nov 13 20:47:58 2011
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
float cpu_gpu_sdiff(             int M, int N, 
                  float * a, int lda, 
                  float *da, int ldda);
void szero_32x32_block(           float *, magma_int_t);
void szero_nbxnb_block(           magma_int_t, float *, magma_int_t);
void magmablas_sinplace_transpose(float *, magma_int_t, magma_int_t);
void magmablas_spermute_long(     float *, magma_int_t, 
                  magma_int_t *, magma_int_t, magma_int_t);
void magmablas_spermute_long2(    float *, magma_int_t, 
                  magma_int_t *, magma_int_t, magma_int_t);
void magmablas_spermute_long3( float *dAT, int lda, 
                               int *ipiv, int nb, int ind );
void magmablas_stranspose(        float *, magma_int_t, 
                  float *, magma_int_t, 
                  magma_int_t, magma_int_t);
void magmablas_stranspose2(       float *, magma_int_t, 
                  float *, magma_int_t, 
                  magma_int_t, magma_int_t);
void magmablas_stranspose2s(float *odata, int ldo,
                       float *idata, int ldi,
                       int m, int n, cudaStream_t *stream );

void magmablas_sgetmatrix_transpose(  int m, int n,
                                      float *dat, int ldda,
                                      float  *ha, int lda,
                                      float  *dB, int lddb, int nb );
void magmablas_sgetmatrix_transpose2( int m, int n,
                                      float **dat, int *ldda,
                                      float  *ha,  int  lda,
                                      float **dB,  int  lddb, int nb,
                                      int num_gpus, cudaStream_t stream[][2] );
void magmablas_ssetmatrix_transpose(  int m, int n,
                                      float  *ha, int lda, 
                                      float *dat, int ldda,
                                      float  *dB, int lddb, int nb );
void magmablas_ssetmatrix_transpose2( int m, int n,
                                      float  *ha,  int  lda, 
                                      float **dat, int *ldda,
                                      float **dB,  int  lddb, int nb,
                                      int num_gpus, cudaStream_t stream[][2] );
void magmablas_sgetmatrix_1D_bcyclic( int m, int n,
                                      float  *da[], int ldda,
                                      float  *ha, int lda,
                                      int num_gpus, int nb );
void magmablas_ssetmatrix_1D_bcyclic( int m, int n,
                                      float  *ha, int lda,
                                      float  *da[], int ldda,
                                      int num_gpus, int nb );

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
void   magmablas_slascl( char type, int kl, int ku,
             float cfrom, float cto,
             int m, int n,
             float *A, int lda, int *info );
void   magmablas_slaset( char uplo, magma_int_t m, magma_int_t n,
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
   * Workspace interface (alphabetical order)
   */
magma_int_t magmablasw_ssymv(char u, magma_int_t N, 
                 float alpha, 
                 float *A, magma_int_t lda, 
                 float *X, magma_int_t incX, 
                 float beta, 
                 float *Y, magma_int_t incY,
                 float *dWork);

#ifdef __cplusplus
}
#endif

#undef PRECISION_s
#endif
