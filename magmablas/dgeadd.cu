/*
    -- MAGMA (version 1.5.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date September 2014

       @generated from zgeadd.cu normal z -> d, Wed Sep 17 15:08:23 2014
       @author Mark Gates
*/
#include "common_magma.h"

#define BLK_X 64
#define BLK_Y 32

/*
    Divides matrix into ceil( m/BLK_X ) x ceil( n/BLK_Y ) blocks.
    Each block has BLK_X threads.
    Each thread loops across one row, updating BLK_Y entries.

    Code similar to dlaset.
*/
__global__
void dgeadd_full(
    int m, int n,
    double alpha,
    const double *dA, int ldda,
    double       *dB, int lddb )
{
    int ind = blockIdx.x*BLK_X + threadIdx.x;
    int iby = blockIdx.y*BLK_Y;
    /* check if full block-column */
    bool full = (iby + BLK_Y <= n);
    /* do only rows inside matrix */
    if ( ind < m ) {
        dA += ind + iby*ldda;
        dB += ind + iby*lddb;
        if ( full ) {
            // full block-column
            #pragma unroll
            for( int j=0; j < BLK_Y; ++j ) {
                dB[j*lddb] = alpha*dA[j*ldda] + dB[j*lddb];
            }
        }
        else {
            // partial block-column
            for( int j=0; j < BLK_Y && iby+j < n; ++j ) {
                dB[j*lddb] = alpha*dA[j*ldda] + dB[j*lddb];
            }
        }
    }
}


/**
    Purpose
    -------
    ZGEADD adds two matrices, dB = alpha*dA + dB.
    
    Arguments
    ---------
    
    @param[in]
    m       INTEGER
            The number of rows of the matrix dA.  M >= 0.
    
    @param[in]
    n       INTEGER
            The number of columns of the matrix dA.  N >= 0.
    
    @param[in]
    alpha   DOUBLE_PRECISION
            The scalar alpha.
            
    @param[in]
    dA      DOUBLE_PRECISION array, dimension (LDDA,N)
            The m by n matrix dA.
    
    @param[in]
    ldda    INTEGER
            The leading dimension of the array dA.  LDDA >= max(1,M).
            
    @param[in,out]
    dB      DOUBLE_PRECISION array, dimension (LDDB,N)
            The m by n matrix dB.
    
    @param[in]
    lddb    INTEGER
            The leading dimension of the array dB.  LDDB >= max(1,M).
    
    @param[in]
    queue   magma_queue_t
            Queue to execute in.

    @ingroup magma_daux2
    ********************************************************************/
extern "C" void
magmablas_dgeadd_q(
    magma_int_t m, magma_int_t n,
    double alpha,
    const double *dA, magma_int_t ldda,
    double       *dB, magma_int_t lddb,
    magma_queue_t queue )
{
    magma_int_t info = 0;
    if ( m < 0 )
        info = -1;
    else if ( n < 0 )
        info = -2;
    else if ( ldda < max(1,m))
        info = -5;
    else if ( lddb < max(1,m))
        info = -7;
    
    if ( info != 0 ) {
        magma_xerbla( __func__, -(info) );
        return;
    }
    
    if ( m == 0 || n == 0 )
        return;
    
    dim3 threads( BLK_X );
    dim3 grid( (m + BLK_X - 1)/BLK_X, (n + BLK_Y - 1)/BLK_Y );
    
    dgeadd_full<<< grid, threads, 0, queue >>>
        ( m, n, alpha, dA, ldda, dB, lddb );
}


/**
    @see magmablas_dgeadd_q
    @ingroup magma_daux2
    ********************************************************************/
extern "C" void
magmablas_dgeadd(
    magma_int_t m, magma_int_t n,
    double alpha,
    const double *dA, magma_int_t ldda,
    double       *dB, magma_int_t lddb )
{
    magmablas_dgeadd_q( m, n, alpha, dA, ldda, dB, lddb, magma_stream );
}
