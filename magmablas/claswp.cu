/*
    -- MAGMA (version 1.3.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2012

       @generated c Wed Nov 14 22:53:48 2012
       
       @author Mathieu Faverge
       @author Ichitaro Yamazaki
       @author Mark Gates
*/
#include "common_magma.h"

// MAX_PIVOTS is maximum number of pivots to apply in each kernel launch
// NTHREADS is number of threads in a block
#define MAX_PIVOTS 32
#define NTHREADS   64

typedef struct {
    cuFloatComplex *dAT;
    int n, lda, j0, npivots;
    int ipiv[MAX_PIVOTS];
} claswp_params_t;


// Matrix A is stored row-wise in dAT.
// Divide matrix A into block-columns of NTHREADS columns each.
// Each GPU block processes one block-column of A.
// Each thread goes down a column of A,
// swapping rows according to pivots stored in params.
__global__ void claswp_kernel( claswp_params_t params )
{
    unsigned int tid = threadIdx.x + blockDim.x*blockIdx.x;
    if( tid < params.n ) {
        int lda = params.lda;
        cuFloatComplex *dAT = params.dAT + tid + params.j0*lda;
        cuFloatComplex *A1  = dAT;
        
        for( int i1 = 0; i1 < params.npivots; ++i1 ) {
            int i2 = params.ipiv[i1];
            cuFloatComplex *A2 = dAT + i2*lda;
            cuFloatComplex temp = *A1;
            *A1 = *A2;
            *A2 = temp;
            A1 += lda;  // A1 = dA + i1*ldx
        }
    }
}


// Launch claswp kernel with ceil( n / NTHREADS ) blocks of NTHREADS threads each.
extern "C" void claswp( claswp_params_t &params )
{
    int blocks = (params.n + NTHREADS - 1) / NTHREADS;
    claswp_kernel<<< blocks, NTHREADS, 0, magma_stream >>>( params );
}


// Swap rows of A, stored row-wise.
// This version updates each entry of ipiv by adding ind.
// It is used in cgetrf, cgetrf_gpu, cgetrf_mgpu, cgetrf_ooc.
extern "C" void
magmablas_cpermute_long2( magma_int_t n, cuFloatComplex *dAT, magma_int_t lda,
                          magma_int_t *ipiv, magma_int_t nb, magma_int_t ind )
{
    for( int k = 0; k < nb; k += MAX_PIVOTS ) {
        int npivots = min( MAX_PIVOTS, nb-k );
        // fields are:             dAT  n  lda  j0       npivots
        claswp_params_t params = { dAT, n, lda, ind + k, npivots };
        for( int j = 0; j < npivots; ++j ) {
            params.ipiv[j] = ipiv[ind + k + j] - k - 1;
            ipiv[ind + k + j] += ind;
        }
        claswp( params );
    }
}


// Swap rows of A, stored row-wise.
// This version assumes ind has already been added to ipiv.
// It is used in cgetrf_mgpu, cgetrf_ooc.
extern "C" void
magmablas_cpermute_long3( cuFloatComplex *dAT, magma_int_t lda,
                          const magma_int_t *ipiv, magma_int_t nb, magma_int_t ind )
{
    for( int k = 0; k < nb; k += MAX_PIVOTS ) {
        int npivots = min( MAX_PIVOTS, nb-k );
        // fields are:             dAT  n    lda  j0       npivots
        claswp_params_t params = { dAT, lda, lda, ind + k, npivots };
        for( int j = 0; j < MAX_PIVOTS; ++j ) {
            params.ipiv[j] = ipiv[ind + k + j] - k - 1 - ind;
        }
        claswp( params );
    }
}


// Swap rows of A, stored row-wise.
// This interface is identical to LAPACK's laswp interface.
// It is used in cgessm, cgetrf_incpiv.
extern "C" void
magmablas_claswp( magma_int_t n, cuFloatComplex *dAT, magma_int_t lda,
                  magma_int_t i1, magma_int_t i2,
                  const magma_int_t *ipiv, magma_int_t inci )
{
    for( int k = i1-1; k < i2; k += MAX_PIVOTS ) {
        int npivots = min( MAX_PIVOTS, i2-k );
        // fields are:             dAT        n  lda  j0 npivots
        claswp_params_t params = { dAT+k*lda, n, lda, 0, npivots };
        for( int j = 0; j < npivots; ++j ) {
            params.ipiv[j] = ipiv[(k+j)*inci] - k - 1;
        }
        claswp( params );
    }
}


// ------------------------------------------------------------
// Extended version has stride in both directions (ldx, ldy)
// to handle both row-wise and column-wise storage.

typedef struct {
    cuFloatComplex *dA;
    int n, ldx, ldy, j0, npivots;
    int ipiv[MAX_PIVOTS];
} claswpx_params_t;


// Matrix A is stored row-wise in dA.
// Divide matrix A into block-columns of NTHREADS columns each.
// Each GPU block processes one block-column of A.
// Each thread goes down a column of A,
// swapping rows according to pivots stored in params.
__global__ void claswpx_kernel( claswpx_params_t params )
{
    unsigned int tid = threadIdx.x + blockDim.x*blockIdx.x;
    if( tid < params.n ) {
        int ldx = params.ldx;
        cuFloatComplex *dA = params.dA + tid*params.ldy + params.j0*ldx;
        cuFloatComplex *A1  = dA;
        
        for( int i1 = 0; i1 < params.npivots; ++i1 ) {
            int i2 = params.ipiv[i1];
            cuFloatComplex *A2 = dA + i2*ldx;
            cuFloatComplex temp = *A1;
            *A1 = *A2;
            *A2 = temp;
            A1 += ldx;  // A1 = dA + i1*ldx
        }
    }
}


// Launch claswpx kernel with ceil( n / NTHREADS ) blocks of NTHREADS threads each.
extern "C" void claswpx( claswpx_params_t &params )
{
    int blocks = (params.n + NTHREADS - 1) / NTHREADS;
    claswpx_kernel<<< blocks, NTHREADS, 0, magma_stream >>>( params );
}


// Swap rows of A.
// For A stored row-wise,    set ldx=lda and ldy=1.
// For A stored column-wise, set ldx=1   and ldy=lda.
// Otherwise, this interface is identical to LAPACK's laswp interface.
extern "C" void
magmablas_claswpx( magma_int_t n, cuFloatComplex *dA, magma_int_t ldx, magma_int_t ldy,
                   magma_int_t i1, magma_int_t i2,
                   const magma_int_t *ipiv, magma_int_t inci )
{
    for( int k = i1-1; k < i2; k += MAX_PIVOTS ) {
        int npivots = min( MAX_PIVOTS, i2-k );
        // fields are:              dA        n  ldx  ldy  j0 npivots
        claswpx_params_t params = { dA+k*ldx, n, ldx, ldy, 0, npivots };
        for( int j = 0; j < npivots; ++j ) {
            params.ipiv[j] = ipiv[(k+j)*inci] - k - 1;
        }
        claswpx( params );
    }
}


// ------------------------------------------------------------
// This version takes d_ipiv on the GPU. Thus it does not pass pivots
// as an argument using a structure, avoiding all the argument size
// limitations of CUDA and OpenCL. It also needs just one kernel launch
// with all the pivots, instead of multiple kernel launches with small
// batches of pivots. On Fermi, it is faster than magmablas_claswp
// (including copying pivots to the GPU).

__global__ void claswp2_kernel( int n, cuFloatComplex *dAT, int lda, int npivots, const magma_int_t* d_ipiv )
{
    unsigned int tid = threadIdx.x + blockDim.x*blockIdx.x;
    if( tid < n ) {
        dAT += tid;
        cuFloatComplex *A1  = dAT;
        
        for( int i1 = 0; i1 < npivots; ++i1 ) {
            int i2 = d_ipiv[i1] - 1;  // Fortran index
            cuFloatComplex *A2 = dAT + i2*lda;
            cuFloatComplex temp = *A1;
            *A1 = *A2;
            *A2 = temp;
            A1 += lda;  // A1 = dA + i1*ldx
        }
    }
}

// Swap rows of A, stored row-wise.
// d_ipiv is vector of pivots stored on the GPU,
// unlike magmablas_claswp where ipiv is stored on the CPU.
// This interface is identical to LAPACK's laswp interface.
extern "C" void
magmablas_claswp2( magma_int_t n, cuFloatComplex* dAT, magma_int_t lda,
                   magma_int_t i1, magma_int_t i2,
                   const magma_int_t *d_ipiv )
{
    int blocks = (n + NTHREADS - 1) / NTHREADS;
    claswp2_kernel<<< blocks, NTHREADS, 0, magma_stream >>>(
        n, dAT + (i1-1)*lda, lda, i2-(i1-1), d_ipiv );
}
