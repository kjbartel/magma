/*
    -- MAGMA (version 1.2.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       June 2012

       @generated s Thu Jun 28 12:31:18 2012

*/
#include "common_magma.h"

#define BLOCK_SIZE 64

typedef struct {
        float *A;
        int n, lda, j0;
        short ipiv[BLOCK_SIZE];
} slaswp_params_t;

typedef struct {
        float *A;
        int n, lda, j0, npivots;
        short ipiv[BLOCK_SIZE];
} slaswp_params_t2;

/*********************************************************
 *
 * LAPACK Swap: permute a set of lines following ipiv
 *
 ********************************************************/
typedef struct {
    float *A;
    int n, ldx, ldy, j0, npivots;
    short ipiv[BLOCK_SIZE];
} slaswpx_params_t;

__global__ void myslaswpx( slaswpx_params_t params )
{
    unsigned int y = threadIdx.x + __mul24(blockDim.x, blockIdx.x);
    unsigned int offset1 = __mul24( y, params.ldy);
    if( y < params.n )
    {
        int ldx = params.ldx;
        float *A = params.A + offset1 + ldx * params.j0;
        float *Ai = A;
        
        for( int i = 0; i < params.npivots; i++ )
        {
            int j = params.ipiv[i];
            float *p2 = A + j*ldx;
            float temp = *Ai;
            *Ai = *p2;
            *p2 = temp;
            Ai += ldx;
        }
    }
}

extern "C" void slaswpx( slaswpx_params_t &params )
{
         int blocksize = 64;
        dim3 blocks = (params.n+blocksize-1) / blocksize;
        myslaswpx<<< blocks, blocksize, 0, magma_stream >>>( params );
}

/*
 * Old version
 */
__global__ void myslaswp2( slaswp_params_t2 params )
{
        unsigned int tid = threadIdx.x + __mul24(blockDim.x, blockIdx.x);
        if( tid < params.n )
        {
                int lda = params.lda;
                float *A = params.A + tid + lda * params.j0;

                for( int i = 0; i < params.npivots; i++ )
                {
                         int j = params.ipiv[i];
                        float *p1 = A + i*lda;
                        float *p2 = A + j*lda;
                        float temp = *p1;
                        *p1 = *p2;
                        *p2 = temp;
                }
        }
}

extern "C" void slaswp2( slaswp_params_t &params );

extern "C" void slaswp3( slaswp_params_t2 &params )
{
         int blocksize = 64;
        dim3 blocks = (params.n+blocksize-1) / blocksize;
        myslaswp2<<< blocks, blocksize, 0, magma_stream >>>( params );
}


extern "C" void 
magmablas_spermute_long2( magma_int_t n, float *dAT, magma_int_t lda,
                          magma_int_t *ipiv, magma_int_t nb, magma_int_t ind )
{
        int k;

        for( k = 0; k < nb-BLOCK_SIZE; k += BLOCK_SIZE )
        {
                //slaswp_params_t params = { dAT, lda, lda, ind + k };
                slaswp_params_t2 params = { dAT, n, lda, ind + k, BLOCK_SIZE };
                for( int j = 0; j < BLOCK_SIZE; j++ )
                {
                        params.ipiv[j] = ipiv[ind + k + j] - k - 1;
                        ipiv[ind + k + j] += ind;
                }
                //slaswp2( params );
                slaswp3( params );
        }

        int num_pivots = nb - k;

        slaswp_params_t2 params = { dAT, n, lda, ind + k, num_pivots};
        for( int j = 0; j < num_pivots; j++ )
        {
            params.ipiv[j] = ipiv[ind + k + j] - k - 1;
            ipiv[ind + k + j] += ind;
        }
        slaswp3( params );
}

extern "C" void 
magmablas_slaswp( magma_int_t n, float *dAT, magma_int_t lda, 
                  magma_int_t i1, magma_int_t i2, magma_int_t *ipiv, magma_int_t inci )
{
  int k;
  
  for( k=(i1-1); k<i2; k+=BLOCK_SIZE )
    {
      int sb = min(BLOCK_SIZE, i2-k);
      //slaswp_params_t params = { dAT, lda, lda, ind + k };
      slaswp_params_t2 params = { dAT+k*lda, n, lda, 0, sb };
      for( int j = 0; j < sb; j++ )
        {
          params.ipiv[j] = ipiv[(k+j)*inci] - k - 1;
        }
      slaswp3( params );
    }
}

extern "C" void 
magmablas_slaswpx( magma_int_t n, float *dAT, magma_int_t ldx, magma_int_t ldy, 
                   magma_int_t i1, magma_int_t i2, magma_int_t *ipiv, magma_int_t inci )
{
  int k;
  
  for( k=(i1-1); k<i2; k+=BLOCK_SIZE )
    {
      int sb = min(BLOCK_SIZE, i2-k);
      //slaswp_params_t params = { dAT, lda, lda, ind + k };
      slaswpx_params_t params = { dAT+k*ldx, n, ldx, ldy, 0, sb };
      for( int j = 0; j < sb; j++ )
        {
          params.ipiv[j] = ipiv[(k+j)*inci] - k - 1;
        }
      slaswpx( params );
    }
}

#undef BLOCK_SIZE
