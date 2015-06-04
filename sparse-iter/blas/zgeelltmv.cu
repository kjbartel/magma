/*
    -- MAGMA (version 1.5.0-beta2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2014

       @precisions normal z -> c d s

*/

#include "common_magma.h"

#if (GPUSHMEM < 200)
   #define BLOCK_SIZE 128
#else
   #define BLOCK_SIZE 512
#endif


// ELL SpMV kernel
//Michael Garland
__global__ void 
zgeelltmv_kernel( int num_rows, 
                 int num_cols,
                 int num_cols_per_row,
                 magmaDoubleComplex alpha, 
                 magmaDoubleComplex *d_val, 
                 magma_index_t *d_colind,
                 magmaDoubleComplex *d_x,
                 magmaDoubleComplex beta, 
                 magmaDoubleComplex *d_y)
{
    int row = blockDim.x * blockIdx.x + threadIdx.x ;
    if(row < num_rows ){
        magmaDoubleComplex dot = MAGMA_Z_MAKE(0.0, 0.0);
        for ( int n = 0; n < num_cols_per_row ; n ++){
            int col = d_colind [ num_rows * n + row ];
            magmaDoubleComplex val = d_val [ num_rows * n + row ];
            if( val != 0)
                dot += val * d_x[col ];
        }
        d_y[ row ] = dot * alpha + beta * d_y [ row ];
    }
}

// shifted ELL SpMV kernel
//Michael Garland
__global__ void 
zgeelltmv_kernel_shift( int num_rows, 
                        int num_cols,
                        int num_cols_per_row,
                        magmaDoubleComplex alpha, 
                        magmaDoubleComplex lambda, 
                        magmaDoubleComplex *d_val, 
                        magma_index_t *d_colind,
                        magmaDoubleComplex *d_x,
                        magmaDoubleComplex beta, 
                        int offset,
                        int blocksize,
                        magma_index_t *add_rows,
                        magmaDoubleComplex *d_y){

    int row = blockDim.x * blockIdx.x + threadIdx.x ;
    if(row < num_rows ){
        magmaDoubleComplex dot = MAGMA_Z_MAKE(0.0, 0.0);
        for ( int n = 0; n < num_cols_per_row ; n ++){
            int col = d_colind [ num_rows * n + row ];
            magmaDoubleComplex val = d_val [ num_rows * n + row ];
            if( val != 0)
                dot += val * d_x[col ];
        }
        if( row<blocksize )
            d_y[ row ] = dot * alpha - lambda 
                    * d_x[ offset+row ] + beta * d_y [ row ];
        else
            d_y[ row ] = dot * alpha - lambda 
                    * d_x[ add_rows[row-blocksize] ] + beta * d_y [ row ];            
    }
}




/*  -- MAGMA (version 1.5.0-beta2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2014

    Purpose
    =======
    
    This routine computes y = alpha *  A^t *  x + beta * y on the GPU.
    Input format is ELL.
    
    Arguments
    =========

    magma_int_t m                   number of rows in A
    magma_int_t n                   number of columns in A 
    magmaDoubleComplex alpha        scalar multiplier
    magmaDoubleComplex *d_val       array containing values of A in ELL
    magma_int_t *d_colind           columnindices of A in ELL
    magmaDoubleComplex *d_x         input vector x
    magmaDoubleComplex beta         scalar multiplier
    magmaDoubleComplex *d_y         input/output vector y

    ======================================================================    */

extern "C" magma_int_t
magma_zgeelltmv(   magma_trans_t transA,
                   magma_int_t m, magma_int_t n,
                   magma_int_t nnz_per_row,
                   magmaDoubleComplex alpha,
                   magmaDoubleComplex *d_val,
                   magma_index_t *d_colind,
                   magmaDoubleComplex *d_x,
                   magmaDoubleComplex beta,
                   magmaDoubleComplex *d_y ){



   dim3 grid( (m+BLOCK_SIZE-1)/BLOCK_SIZE, 1, 1);

   zgeelltmv_kernel<<< grid, BLOCK_SIZE, 0, magma_stream >>>
                  ( m, n, nnz_per_row, alpha, d_val, d_colind, d_x, beta, d_y );


   return MAGMA_SUCCESS;
}


/*  -- MAGMA (version 1.5.0-beta2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2014

    Purpose
    =======
    
    This routine computes y = alpha *( A - lambda I ) * x + beta * y on the GPU.
    Input format is ELL.
    
    Arguments
    =========

    magma_int_t m                   number of rows in A
    magma_int_t n                   number of columns in A 
    magmaDoubleComplex alpha        scalar multiplier
    magmaDoubleComplex lambda       scalar multiplier
    magmaDoubleComplex *d_val       array containing values of A in ELL
    magma_int_t *d_colind           columnindices of A in ELL
    magmaDoubleComplex *d_x         input vector x
    magmaDoubleComplex beta         scalar multiplier
    magmaDoubleComplex *d_y         input/output vector y

    ======================================================================    */

extern "C" magma_int_t
magma_zgeelltmv_shift( magma_trans_t transA,
                       magma_int_t m, magma_int_t n,
                       magma_int_t nnz_per_row,
                       magmaDoubleComplex alpha,
                       magmaDoubleComplex lambda,
                       magmaDoubleComplex *d_val,
                       magma_index_t *d_colind,
                       magmaDoubleComplex *d_x,
                       magmaDoubleComplex beta,
                       int offset,
                       int blocksize,
                       magma_index_t *add_rows,
                       magmaDoubleComplex *d_y ){



   dim3 grid( (m+BLOCK_SIZE-1)/BLOCK_SIZE, 1, 1);
   magmaDoubleComplex tmp_shift;
   //magma_zsetvector(1,&lambda,1,&tmp_shift,1); 
   tmp_shift = lambda;
   zgeelltmv_kernel_shift<<< grid, BLOCK_SIZE, 0, magma_stream >>>
                  ( m, n, nnz_per_row, alpha, tmp_shift, d_val, d_colind, d_x, 
                            beta, offset, blocksize, add_rows, d_y );


   return MAGMA_SUCCESS;
}


