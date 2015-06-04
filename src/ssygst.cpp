/*
    -- MAGMA (version 1.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2011

       @author Raffaele Solca

       @generated s Sun Nov 13 20:48:30 2011
*/

#include "common_magma.h"

#define A(i, j) (a+(j)*lda + (i))
#define B(i, j) (b+(j)*ldb + (i))

#define dA(i, j) (dw+(j)*ldda + (i))
#define dB(i, j) (dw+n*ldda+(j)*lddb + (i))

extern "C" magma_int_t
magma_ssygst(magma_int_t itype, char uplo, magma_int_t n,
             float *a, magma_int_t lda,
             float *b, magma_int_t ldb, magma_int_t *info)
{
/*
  -- MAGMA (version 1.1) --
     Univ. of Tennessee, Knoxville
     Univ. of California, Berkeley
     Univ. of Colorado, Denver
     November 2011

 
   Purpose
   =======
   
   SSYGST reduces a real symmetric-definite generalized
   eigenproblem to standard form.
   
   If ITYPE = 1, the problem is A*x = lambda*B*x,
   and A is overwritten by inv(U**T)*A*inv(U) or inv(L)*A*inv(L**T)
   
   If ITYPE = 2 or 3, the problem is A*B*x = lambda*x or
   B*A*x = lambda*x, and A is overwritten by U*A*U**T or L**T*A*L.
   
   B must have been previously factorized as U**T*U or L*L**T by SPOTRF.
   
   Arguments
   =========
   
   ITYPE   (input) INTEGER
           = 1: compute inv(U**T)*A*inv(U) or inv(L)*A*inv(L**T);
           = 2 or 3: compute U*A*U**T or L**T*A*L.
   
   UPLO    (input) CHARACTER*1
           = 'U':  Upper triangle of A is stored and B is factored as
                   U**T*U;
           = 'L':  Lower triangle of A is stored and B is factored as
                   L*L**T.
   
   N       (input) INTEGER
           The order of the matrices A and B.  N >= 0.
   
   A       (input/output) COMPLEX*16 array, dimension (LDA,N)
           On entry, the symmetric matrix A.  If UPLO = 'U', the leading
           N-by-N upper triangular part of A contains the upper
           triangular part of the matrix A, and the strictly lower
           triangular part of A is not referenced.  If UPLO = 'L', the
           leading N-by-N lower triangular part of A contains the lower
           triangular part of the matrix A, and the strictly upper
           triangular part of A is not referenced.
   
           On exit, if INFO = 0, the transformed matrix, stored in the
           same format as A.
   
   LDA     (input) INTEGER
           The leading dimension of the array A.  LDA >= max(1,N).
   
   B       (input) COMPLEX*16 array, dimension (LDB,N)
           The triangular factor from the Cholesky factorization of B,
           as returned by SPOTRF.
   
   LDB     (input) INTEGER
           The leading dimension of the array B.  LDB >= max(1,N).
   
   INFO    (output) INTEGER
           = 0:  successful exit
           < 0:  if INFO = -i, the i-th argument had an illegal value
   
   =====================================================================*/
  
  char uplo_[2] = {uplo, 0};
  magma_int_t        nb;
  magma_int_t        k, kb, kb2;
  float    zone  = MAGMA_S_ONE;
  float    mzone  = MAGMA_S_NEG_ONE;
  float    zhalf  = MAGMA_S_HALF;
  float    mzhalf  = MAGMA_S_NEG_HALF;
  float   *dw;
  magma_int_t        ldda = n;
  magma_int_t        lddb = n;
  float             done  = (float) 1.0;
  long int           upper = lapackf77_lsame(uplo_, "U");
  
  /* Test the input parameters. */
  *info = 0;
  if (itype<1 || itype>3){
    *info = -1;
  }else if ((! upper) && (! lapackf77_lsame(uplo_, "L"))) {
    *info = -2;
  } else if (n < 0) {
    *info = -3;
  } else if (ldda < max(1,n)) {
    *info = -5;
  }else if (lddb < max(1,n)) {
    *info = -7;
  }
    if (*info != 0) {
        magma_xerbla( __func__, -(*info) );
        return MAGMA_ERR_ILLEGAL_VALUE;
    }
  
  /* Quick return */
  if ( n == 0 )
    return MAGMA_SUCCESS;
  
  if (cudaSuccess != cudaMalloc( (void**)&dw, 2*n*n*sizeof(float) ) ) {
    *info = -6;
    return MAGMA_ERR_CUBLASALLOC;
  }
  
  nb = magma_get_ssygst_nb(n);
  
  static cudaStream_t stream[2];
  cudaStreamCreate(&stream[0]);
  cudaStreamCreate(&stream[1]);

  cublasSetMatrix(n, n, sizeof(float), A(0, 0), lda, dA(0, 0), ldda);
  cublasSetMatrix(n, n, sizeof(float), B(0, 0), ldb, dB(0, 0), lddb);
  
  /* Use hybrid blocked code */
    
    if (itype==1) {
      if (upper) {
        
        /* Compute inv(U')*A*inv(U) */
        
        for(k = 0; k<n; k+=nb){
          kb = min(n-k,nb);
          kb2= min(n-k-nb,nb);
          
          /* Update the upper triangle of A(k:n,k:n) */
          
          lapackf77_shegs2( &itype, uplo_, &kb, A(k,k), &lda, B(k,k), &ldb, info);
          
          cudaMemcpy2DAsync(dA(k, k), ldda * sizeof(float),
                             A(k, k), lda  * sizeof(float),
                            sizeof(float)*kb, kb,
                            cudaMemcpyHostToDevice, stream[0]);
          
          if(k+kb<n){
            
            cublasStrsm(MagmaLeft, MagmaUpper, MagmaTrans, MagmaNonUnit,
                        kb, n-k-kb,
                        zone, dB(k,k), lddb, 
                        dA(k,k+kb), ldda); 
            
            cudaStreamSynchronize(stream[0]);
            
            cublasSsymm(MagmaLeft, MagmaUpper,
                        kb, n-k-kb,
                        mzhalf, dA(k,k), ldda,
                        dB(k,k+kb), lddb,
                        zone, dA(k, k+kb), ldda);
            
            cublasSsyr2k(MagmaUpper, MagmaTrans,
                         n-k-kb, kb,
                         mzone, dA(k,k+kb), ldda,
                         dB(k,k+kb), lddb,
                         done, dA(k+kb,k+kb), ldda);
            
            cudaMemcpy2DAsync(  A(k+kb, k+kb), lda*sizeof(float),
                              dA(k+kb, k+kb), ldda*sizeof(float),
                              sizeof(float)*kb2, kb2,
                              cudaMemcpyDeviceToHost, stream[1]);
            
            cublasSsymm(MagmaLeft, MagmaUpper,
                        kb, n-k-kb,
                        mzhalf, dA(k,k), ldda,
                        dB(k,k+kb), lddb,
                        zone, dA(k, k+kb), ldda);
            
            cublasStrsm(MagmaRight, MagmaUpper, MagmaNoTrans, MagmaNonUnit,
                        kb, n-k-kb,
                        zone ,dB(k+kb,k+kb), lddb,
                        dA(k,k+kb), ldda);
          
            cudaStreamSynchronize(stream[1]);
            
          }
        
        }
        
        cudaStreamSynchronize(stream[0]);
        
      } else {
        
        /* Compute inv(L)*A*inv(L') */
        
        for(k = 0; k<n; k+=nb){
          kb= min(n-k,nb);
          kb2= min(n-k-nb,nb);
          
          /* Update the lower triangle of A(k:n,k:n) */
          
          lapackf77_shegs2( &itype, uplo_, &kb, A(k,k), &lda, B(k,k), &ldb, info);
          
          cudaMemcpy2DAsync(dA(k, k), ldda * sizeof(float),
                            A(k, k), lda  * sizeof(float),
                            sizeof(float)*kb, kb,
                            cudaMemcpyHostToDevice, stream[0]);
          
          if(k+kb<n){
            
            cublasStrsm(MagmaRight, MagmaLower, MagmaTrans, MagmaNonUnit,
                        n-k-kb, kb,
                        zone, dB(k,k), lddb, 
                        dA(k+kb,k), ldda);
            
            cudaStreamSynchronize(stream[0]);
            
            cublasSsymm(MagmaRight, MagmaLower,
                        n-k-kb, kb,
                        mzhalf, dA(k,k), ldda,
                        dB(k+kb,k), lddb,
                        zone, dA(k+kb, k), ldda);
            
            cublasSsyr2k(MagmaLower, MagmaNoTrans,
                         n-k-kb, kb,
                         mzone, dA(k+kb,k), ldda,
                         dB(k+kb,k), lddb,
                         done, dA(k+kb,k+kb), ldda);
            
            cudaMemcpy2DAsync( A(k+kb, k+kb), lda *sizeof(float),
                              dA(k+kb, k+kb), ldda*sizeof(float),
                              sizeof(float)*kb2, kb2,
                              cudaMemcpyDeviceToHost, stream[1]);
            
            cublasSsymm(MagmaRight, MagmaLower,
                        n-k-kb, kb,
                        mzhalf, dA(k,k), ldda,
                        dB(k+kb,k), lddb,
                        zone, dA(k+kb, k), ldda);
            
            cublasStrsm(MagmaLeft, MagmaLower, MagmaNoTrans, MagmaNonUnit,
                        n-k-kb, kb,
                        zone, dB(k+kb,k+kb), lddb, 
                        dA(k+kb,k), ldda);            
          }

          cudaStreamSynchronize(stream[1]);
          
        }
        
      }
      
      cudaStreamSynchronize(stream[0]);
      
    } else {
      
      if (upper) {
        
        /* Compute U*A*U' */
        
        for(k = 0; k<n; k+=nb){
          kb= min(n-k,nb);
          
          cudaMemcpy2DAsync( A(k, k), lda *sizeof(float),
                            dA(k, k), ldda*sizeof(float),
                            sizeof(float)*kb, kb,
                            cudaMemcpyDeviceToHost, stream[0]);
          
          /* Update the upper triangle of A(1:k+kb-1,1:k+kb-1) */
          if(k>0){
            
            cublasStrmm(MagmaLeft, MagmaUpper, MagmaNoTrans, MagmaNonUnit,
                        k, kb,
                        zone ,dB(0,0), lddb,
                        dA(0,k), ldda);
            
            cublasSsymm(MagmaRight, MagmaUpper,
                        k, kb,
                        zhalf, dA(k,k), ldda,
                        dB(0,k), lddb,
                        zone, dA(0, k), ldda);
            
            cudaStreamSynchronize(stream[1]);
            
            cublasSsyr2k(MagmaUpper, MagmaNoTrans,
                         k, kb,
                         zone, dA(0,k), ldda,
                         dB(0,k), lddb,
                         done, dA(0,0), ldda);
            
            cublasSsymm(MagmaRight, MagmaUpper,
                        k, kb,
                        zhalf, dA(k,k), ldda,
                        dB(0,k), lddb,
                        zone, dA(0, k), ldda);
            
            cublasStrmm(MagmaRight, MagmaUpper, MagmaTrans, MagmaNonUnit,
                        k, kb,
                        zone, dB(k,k), lddb, 
                        dA(0,k), ldda);
            
          }
          
          cudaStreamSynchronize(stream[0]);
          
          lapackf77_shegs2( &itype, uplo_, &kb, A(k, k), &lda, B(k, k), &ldb, info);
          
          cudaMemcpy2DAsync(dA(k, k), ldda * sizeof(float),
                             A(k, k), lda  * sizeof(float),
                            sizeof(float)*kb, kb,
                            cudaMemcpyHostToDevice, stream[1]);
          
        }
        
        cudaStreamSynchronize(stream[1]);
        
      } else {
        
        /* Compute L'*A*L */
        
        for(k = 0; k<n; k+=nb){
          kb= min(n-k,nb);
          
          cudaMemcpy2DAsync( A(k, k), lda *sizeof(float),
                            dA(k, k), ldda*sizeof(float),
                            sizeof(float)*kb, kb,
                            cudaMemcpyDeviceToHost, stream[0]);
          
          /* Update the lower triangle of A(1:k+kb-1,1:k+kb-1) */
          if(k>0){ 
            
            cublasStrmm(MagmaRight, MagmaLower, MagmaNoTrans, MagmaNonUnit,
                        kb, k,
                        zone ,dB(0,0), lddb,
                        dA(k,0), ldda);
            
            cublasSsymm(MagmaLeft, MagmaLower,
                        kb, k,
                        zhalf, dA(k,k), ldda,
                        dB(k,0), lddb,
                        zone, dA(k, 0), ldda);
            
            cudaStreamSynchronize(stream[1]);
            
            cublasSsyr2k(MagmaLower, MagmaTrans,
                         k, kb,
                         zone, dA(k,0), ldda,
                         dB(k,0), lddb,
                         done, dA(0,0), ldda);
            
            cublasSsymm(MagmaLeft, MagmaLower,
                        kb, k,
                        zhalf, dA(k,k), ldda,
                        dB(k,0), lddb,
                        zone, dA(k, 0), ldda);
            
            cublasStrmm(MagmaLeft, MagmaLower, MagmaTrans, MagmaNonUnit,
                        kb, k,
                        zone, dB(k,k), lddb, 
                        dA(k,0), ldda);
          }
          
          cudaStreamSynchronize(stream[0]);
          
          lapackf77_shegs2( &itype, uplo_, &kb, A(k,k), &lda, B(k,k), &ldb, info);
          
          cudaMemcpy2DAsync(dA(k, k), ldda * sizeof(float),
                             A(k, k), lda  * sizeof(float),
                            sizeof(float)*kb, kb,
                            cudaMemcpyHostToDevice, stream[1]);
        }
        
        cudaStreamSynchronize(stream[1]);
        
      }
  }
  
  cublasGetMatrix(n, n, sizeof(float), dA(0, 0), ldda, A(0, 0), lda);

  cudaStreamDestroy(stream[0]);
  cudaStreamDestroy(stream[1]); 
  
  cublasFree(dw);
  
  return MAGMA_SUCCESS;
} /* magma_ssygst_gpu */
