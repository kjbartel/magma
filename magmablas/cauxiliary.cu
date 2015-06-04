/*
    -- MAGMA (version 1.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2011

       @generated c Sun Nov 13 20:48:36 2011

*/
#include "common_magma.h"

/* ////////////////////////////////////////////////////////////////////////////
   -- This is an auxiliary routine called from cgehrd.  The routine is called
      in 16 blocks, 32 thread per block and initializes to zero the 1st 
      32x32 block of A.
*/

__global__ void cset_to_zero(cuFloatComplex *A, int lda){
    int ind = blockIdx.x*lda + threadIdx.x;
    
    A += ind;
    A[0] = MAGMA_C_ZERO;
//   A[16*lda] = 0.;
}

__global__ void cset_nbxnb_to_zero(int nb, cuFloatComplex *A, int lda){
   int ind = blockIdx.x*lda + threadIdx.x, i, j;

   A += ind;
   for(i=0; i<nb; i+=32){
     for(j=0; j<nb; j+=32)
         A[j] = MAGMA_C_ZERO;
     A += 32*lda;
   }
}

void czero_32x32_block(cuFloatComplex *A, int lda)
{
  // cset_to_zero<<< 16, 32, 0, magma_stream >>>(A, lda);
  cset_to_zero<<< 32, 32, 0, magma_stream >>>(A, lda);
}

void czero_nbxnb_block(int nb, cuFloatComplex *A, int lda)
{
  cset_nbxnb_to_zero<<< 32, 32, 0, magma_stream >>>(nb, A, lda);
}

/* ////////////////////////////////////////////////////////////////////////////
   -- GPU kernel for initializing a matrix by 0
*/
#define claset_threads 64

__global__ void claset(int m, int n, cuFloatComplex *A, int lda){
   int ibx = blockIdx.x * claset_threads;
   int iby = blockIdx.y * 32;

   int ind = ibx + threadIdx.x;

   A += ind + __mul24(iby, lda);

   #pragma unroll
   for(int i=0; i<32; i++)
     if (iby+i < n && ind < m)
        A[i*lda] = MAGMA_C_ZERO;
}

__global__ void clasetlower(int m, int n, cuFloatComplex *A, int lda){
   int ibx = blockIdx.x * claset_threads;
   int iby = blockIdx.y * 32;

   int ind = ibx + threadIdx.x;

   A += ind + __mul24(iby, lda);

   #pragma unroll
   for(int i=0; i<32; i++)
     if (iby+i < n && ind < m && ind > i+iby)
        A[i*lda] = MAGMA_C_ZERO;
}

__global__ void clasetupper(int m, int n, cuFloatComplex *A, int lda){
   int ibx = blockIdx.x * claset_threads;
   int iby = blockIdx.y * 32;

   int ind = ibx + threadIdx.x;

   A += ind + __mul24(iby, lda);

   #pragma unroll
   for(int i=0; i<32; i++)
     if (iby+i < n && ind < m && ind < i+iby)
        A[i*lda] = MAGMA_C_ZERO;
}


/* ////////////////////////////////////////////////////////////////////////////
   -- Set the m x n matrix pointed by A to 0 on the GPU.
*/
extern "C" void
magmablas_claset(char uplo, magma_int_t m, magma_int_t n, 
                 cuFloatComplex *A, magma_int_t lda)
{
   dim3 threads(claset_threads, 1, 1);
   dim3 grid(m/claset_threads+(m % claset_threads != 0), n/32+(n%32!=0));

   if (m!=0 && n !=0)
     if (uplo == MagmaLower)
        clasetlower<<< grid, threads, 0, magma_stream >>> (m, n, A, lda);        
     else if (uplo == MagmaUpper)
        clasetupper<<< grid, threads, 0, magma_stream >>> (m, n, A, lda);
     else
        claset<<< grid, threads, 0, magma_stream >>> (m, n, A, lda);
}

/* ////////////////////////////////////////////////////////////////////////////
   -- Given two matrices, 'a' on the CPU and 'da' on the GPU, this function
      returns the Frobenious norm of the difference of the two matrices.
      The function is used for debugging.
*/
float cpu_gpu_cdiff(int M, int N, cuFloatComplex * a, int lda, cuFloatComplex *da, int ldda)
{
  int one = 1, j;
  cuFloatComplex mone = MAGMA_C_NEG_ONE;
  float  work[1];
  cuFloatComplex *ha = (cuFloatComplex*)malloc( M * N * sizeof(cuFloatComplex));
  float res;

  cublasGetMatrix(M, N, sizeof(cuFloatComplex), da, ldda, ha, M);
  for(j=0; j<N; j++)
    blasf77_caxpy(&M, &mone, a+j*lda, &one, ha+j*M, &one);
  res = lapackf77_clange("f", &M, &N, ha, &M, work);

  free(ha);
  return res;
}

/* ////////////////////////////////////////////////////////////////////////////
 -- GPU kernel for setting 0 in the nb-1 upper subdiagonals and 1 in the diagonal
    @author Raffaele Solca
 */
__global__ void csetdiag1subdiag0_L(int k, cuFloatComplex *A, int lda){
  
  int nb = blockDim.x;
  int ibx = blockIdx.x * nb;
  
  int ind = ibx + threadIdx.x + 1;
  
  A += ind - nb + __mul24((ibx), lda);
  
  cuFloatComplex tmp = MAGMA_C_ZERO;
  if(threadIdx.x == nb-1)
    tmp = MAGMA_C_ONE;
  
#pragma unroll
  for(int i=0; i<nb; i++)
    if (ibx+i < k && ind + i  >= nb){
      A[i*(lda+1)] = tmp;
    }
  
}

/* ////////////////////////////////////////////////////////////////////////////
 -- GPU kernel for setting 0 in the nb-1 lower subdiagonals and 1 in the diagonal
    @author Raffaele Solca
 */

__global__ void csetdiag1subdiag0_U(int k, cuFloatComplex *A, int lda){

  int nb = blockDim.x;
  int ibx = blockIdx.x * nb;

  int ind = ibx + threadIdx.x;

  A += ind + __mul24((ibx), lda);

  cuFloatComplex tmp = MAGMA_C_ZERO;
  if(threadIdx.x == 0)
    tmp = MAGMA_C_ONE;

#pragma unroll
  for(int i=0; i<nb; i++)
    if (ibx+i < k && ind + i < k){
      A[i*(lda+1)] = tmp;
    }

}

/* ////////////////////////////////////////////////////////////////////////////
 -- Set 1s in the diagonal and 0s in the nb-1 lower (UPLO='U') or 
    upper (UPLO='L') subdiagonals 
    @author Raffaele Solca
 */
extern "C" void
magmablas_csetdiag1subdiag0(char uplo, magma_int_t k, magma_int_t nb,
                 cuFloatComplex *A, magma_int_t lda)
{
  
  dim3 threads(nb, 1, 1);
  dim3 grid((k-1)/nb+1);
  if(k>lda)  
    fprintf(stderr,"wrong second argument of csetdiag1subdiag0");
  if(uplo == MagmaLower)
    csetdiag1subdiag0_L<<< grid, threads, 0, magma_stream >>> (k, A, lda);
  else if(uplo == MagmaUpper){
    csetdiag1subdiag0_U<<< grid, threads, 0, magma_stream >>> (k, A, lda);
  }
  else 
    fprintf(stderr,"wrong first argument of csetdiag1subdiag0");

  return;
}
