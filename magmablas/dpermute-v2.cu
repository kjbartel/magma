/*
    -- MAGMA (version 1.2.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       May 2012

       @generated d Tue May 15 18:18:01 2012

*/
#include "common_magma.h"

#define BLOCK_SIZE 64

typedef struct {
        double *A;
        int n, lda, j0;
        short ipiv[BLOCK_SIZE];
} dlaswp_params_t;

typedef struct {
        double *A;
        int n, lda, j0, npivots;
        short ipiv[BLOCK_SIZE];
} dlaswp_params_t2;

/*********************************************************
 *
 * LAPACK Swap: permute a set of lines following ipiv
 *
 ********************************************************/
typedef struct {
    double *A;
    int n, ldx, ldy, j0, npivots;
    short ipiv[BLOCK_SIZE];
} dlaswpx_params_t;

__global__ void mydlaswpx( dlaswpx_params_t params )
{
    unsigned int y = threadIdx.x + __mul24(blockDim.x, blockIdx.x);
    unsigned int offset1 = __mul24( y, params.ldy);
    if( y < params.n )
    {
        int ldx = params.ldx;
        double *A = params.A + offset1 + ldx * params.j0;
        double *Ai = A;
        
        for( int i = 0; i < params.npivots; i++ )
        {
            int j = params.ipiv[i];
            double *p2 = A + j*ldx;
            double temp = *Ai;
            *Ai = *p2;
            *p2 = temp;
            Ai += ldx;
        }
    }
}

extern "C" void dlaswpx( dlaswpx_params_t &params )
{
         int blocksize = 64;
        dim3 blocks = (params.n+blocksize-1) / blocksize;
        mydlaswpx<<< blocks, blocksize, 0, magma_stream >>>( params );
}

/*
 * Old version
 */
__global__ void mydlaswp2( dlaswp_params_t2 params )
{
        unsigned int tid = threadIdx.x + __mul24(blockDim.x, blockIdx.x);
        if( tid < params.n )
        {
                int lda = params.lda;
                double *A = params.A + tid + lda * params.j0;

                for( int i = 0; i < params.npivots; i++ )
                {
                         int j = params.ipiv[i];
                        double *p1 = A + i*lda;
                        double *p2 = A + j*lda;
                        double temp = *p1;
                        *p1 = *p2;
                        *p2 = temp;
                }
        }
}

extern "C" void dlaswp2( dlaswp_params_t &params );

extern "C" void dlaswp3( dlaswp_params_t2 &params )
{
         int blocksize = 64;
        dim3 blocks = (params.n+blocksize-1) / blocksize;
        mydlaswp2<<< blocks, blocksize, 0, magma_stream >>>( params );
}


extern "C" void 
magmablas_dpermute_long2( double *dAT, int lda, int *ipiv, int nb, int ind )
{
        int k;

        for( k = 0; k < nb-BLOCK_SIZE; k += BLOCK_SIZE )
        {
                //dlaswp_params_t params = { dAT, lda, lda, ind + k };
                dlaswp_params_t2 params = { dAT, lda, lda, ind + k, BLOCK_SIZE };
                for( int j = 0; j < BLOCK_SIZE; j++ )
                {
                        params.ipiv[j] = ipiv[ind + k + j] - k - 1;
                        ipiv[ind + k + j] += ind;
                }
                //dlaswp2( params );
                dlaswp3( params );
        }

        int num_pivots = nb - k;

        dlaswp_params_t2 params = { dAT, lda, lda, ind + k, num_pivots};
        for( int j = 0; j < num_pivots; j++ )
        {
            params.ipiv[j] = ipiv[ind + k + j] - k - 1;
            ipiv[ind + k + j] += ind;
        }
        dlaswp3( params );
}

extern "C" void 
magmablas_dlaswp( int n, double *dAT, int lda, 
                  int i1, int i2, int *ipiv, int inci )
{
  int k;
  
  for( k=(i1-1); k<i2; k+=BLOCK_SIZE )
    {
      int sb = min(BLOCK_SIZE, i2-k);
      //dlaswp_params_t params = { dAT, lda, lda, ind + k };
      dlaswp_params_t2 params = { dAT+k*lda, n, lda, 0, sb };
      for( int j = 0; j < sb; j++ )
        {
          params.ipiv[j] = ipiv[(k+j)*inci] - k - 1;
        }
      dlaswp3( params );
    }
}

extern "C" void 
magmablas_dlaswpx( int n, double *dAT, int ldx, int ldy, 
                   int i1, int i2, int *ipiv, int inci )
{
  int k;
  
  for( k=(i1-1); k<i2; k+=BLOCK_SIZE )
    {
      int sb = min(BLOCK_SIZE, i2-k);
      //dlaswp_params_t params = { dAT, lda, lda, ind + k };
      dlaswpx_params_t params = { dAT+k*ldx, n, ldx, ldy, 0, sb };
      for( int j = 0; j < sb; j++ )
        {
          params.ipiv[j] = ipiv[(k+j)*inci] - k - 1;
        }
      dlaswpx( params );
    }
}

#undef BLOCK_SIZE
