/*
    -- MAGMA (version 1.5.0-beta2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2014

       @generated from magma_z_blaswrapper.cpp normal z -> d, Fri May 30 10:41:36 2014
       @author Hartwig Anzt

*/

//#include <cusparse_v2.h>

#include "common_magma.h"
#include "../../include/magmablas.h"
#include "../include/magmasparse_types.h"
#include "../include/magmasparse.h"




/*  -- MAGMA (version 1.5.0-beta2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2014

    Purpose
    =======

    For a given input matrix A and vectors x, y and scalars alpha, beta
    the wrapper determines the suitable SpMV computing
              y = alpha * A * x + beta * y.  
    Arguments
    =========

    double alpha      scalar alpha
    magma_d_sparse_matrix A       sparse matrix A    
    magma_d_vector x              input vector x  
    magmaDoublComplex beta        scalar beta
    magma_d_vector y              input vector y      

    ========================================================================  */

magma_int_t
magma_d_spmv(     double alpha, magma_d_sparse_matrix A, 
                magma_d_vector x, double beta, magma_d_vector y )
{
    if( A.memory_location != x.memory_location || 
                            x.memory_location != y.memory_location ){
    printf("error: linear algebra objects are not located in same memory!\n");
    printf("memory locations are: %d   %d   %d\n", 
                    A.memory_location, x.memory_location, y.memory_location );
    return MAGMA_ERR_INVALID_PTR;
    }

    // DEV case
    if( A.memory_location == Magma_DEV ){
        if( A.num_cols == x.num_rows ){
             if( A.storage_type == Magma_CSR 
                            || A.storage_type == Magma_CSRL 
                            || A.storage_type == Magma_CSRU ){
                 //printf("using CSR kernel for SpMV: ");
                 magma_dgecsrmv( MagmaNoTrans, A.num_rows, A.num_cols, alpha, 
                                 A.val, A.row, A.col, x.val, beta, y.val );
                 //printf("done.\n");
                 return MAGMA_SUCCESS;
             }
             else if( A.storage_type == Magma_ELLPACK ){
                 //printf("using ELLPACK kernel for SpMV: ");
                 magma_dgeellmv( MagmaNoTrans, A.num_rows, A.num_cols, 
                    A.max_nnz_row, alpha, A.val, A.col, x.val, beta, y.val );
                 //printf("done.\n");
                 return MAGMA_SUCCESS;
             }
             else if( A.storage_type == Magma_ELL ){
                 //printf("using ELL kernel for SpMV: ");
                 magma_dgeelltmv( MagmaNoTrans, A.num_rows, A.num_cols, 
                    A.max_nnz_row, alpha, A.val, A.col, x.val, beta, y.val );
                 //printf("done.\n");
                 return MAGMA_SUCCESS;
             }
             else if( A.storage_type == Magma_ELLRT ){
                 //printf("using ELLRT kernel for SpMV: ");
                 magma_dgeellrtmv( MagmaNoTrans, A.num_rows, A.num_cols, 
                            A.max_nnz_row, alpha, A.val, A.col, A.row, x.val, 
                         beta, y.val, A.alignment, A.blocksize );
                 //printf("done.\n");
                 return MAGMA_SUCCESS;
             }
             else if( A.storage_type == Magma_SELLC ){
                 //printf("using SELLC kernel for SpMV: ");
                 magma_dgesellcmv( MagmaNoTrans, A.num_rows, A.num_cols, 
                    A.blocksize, A.numblocks, A.alignment, 
                    alpha, A.val, A.col, A.row, x.val, beta, y.val );

                 //printf("done.\n");
                 return MAGMA_SUCCESS;
             }
             else if( A.storage_type == Magma_SELLP ){
                 //printf("using SELLP kernel for SpMV: ");
                 magma_dgesellpmv( MagmaNoTrans, A.num_rows, A.num_cols, 
                    A.blocksize, A.numblocks, A.alignment, 
                    alpha, A.val, A.col, A.row, x.val, beta, y.val );

                 //printf("done.\n");
                 return MAGMA_SUCCESS;
             }
             else if( A.storage_type == Magma_DENSE ){
                 //printf("using DENSE kernel for SpMV: ");
                 magmablas_dgemv( MagmaNoTrans, A.num_rows, A.num_cols, alpha, 
                                 A.val, A.num_rows, x.val, 1, beta,  y.val, 1 );
                 //printf("done.\n");
                 return MAGMA_SUCCESS;
             }
/*             else if( A.storage_type == Magma_BCSR ){
                 //printf("using CUSPARSE BCSR kernel for SpMV: ");
                // CUSPARSE context //
                cusparseHandle_t cusparseHandle = 0;
                cusparseStatus_t cusparseStatus;
                cusparseStatus = cusparseCreate(&cusparseHandle);
                cusparseMatDescr_t descr = 0;
                cusparseStatus = cusparseCreateMatDescr(&descr);
                // end CUSPARSE context //
                cusparseDirection_t dirA = CUSPARSE_DIRECTION_ROW;
                int mb = (A.num_rows + A.blocksize-1)/A.blocksize;
                int nb = (A.num_cols + A.blocksize-1)/A.blocksize;
                cusparseDbsrmv( cusparseHandle, dirA, 
                    CUSPARSE_OPERATION_NON_TRANSPOSE, mb, nb, A.numblocks, 
                    &alpha, descr, A.val, A.row, A.col, A.blocksize, x.val, 
                    &beta, y.val );
                 //printf("done.\n");
                 return MAGMA_SUCCESS;
             }*/
             else {
                 printf("error: format not supported.\n");
                 return MAGMA_ERR_NOT_SUPPORTED;
             }
        }
        else if( A.num_cols < x.num_rows ){
            magma_int_t num_vecs = x.num_rows / A.num_cols;
            if( A.storage_type == Magma_CSR ){
                 //printf("using CSR kernel for SpMV: ");
                 magma_dmgecsrmv( MagmaNoTrans, A.num_rows, A.num_cols, 
                    num_vecs, alpha, A.val, A.row, A.col, x.val, beta, y.val );
                 //printf("done.\n");
                 return MAGMA_SUCCESS;
             }
             else if( A.storage_type == Magma_ELLPACK ){
                 //printf("using ELLPACK kernel for SpMV: ");
                 magma_dmgeellmv( MagmaNoTrans, A.num_rows, A.num_cols, 
            num_vecs, A.max_nnz_row, alpha, A.val, A.col, x.val, beta, y.val );
                 //printf("done.\n");
                 return MAGMA_SUCCESS;
             }
             else if( A.storage_type == Magma_ELL ){
                 //printf("using ELL kernel for SpMV: ");
                 magma_dmgeelltmv( MagmaNoTrans, A.num_rows, A.num_cols, 
                        num_vecs, A.max_nnz_row, alpha, A.val, 
                        A.col, x.val, beta, y.val );
                 //printf("done.\n");
                 return MAGMA_SUCCESS;
             }else if( A.storage_type == Magma_SELLP ){
                 //printf("using SELLP kernel for SpMV: ");
                 magma_dmgesellpmv( MagmaNoTrans, A.num_rows, A.num_cols, 
                    num_vecs, A.blocksize, A.numblocks, A.alignment, 
                    alpha, A.val, A.col, A.row, x.val, beta, y.val );

                 //printf("done.\n");
                 return MAGMA_SUCCESS;
             }/*
             if( A.storage_type == Magma_DENSE ){
                 //printf("using DENSE kernel for SpMV: ");
                 magmablas_dmgemv( MagmaNoTrans, A.num_rows, A.num_cols, 
                            num_vecs, alpha, A.val, A.num_rows, x.val, 1, 
                            beta,  y.val, 1 );
                 //printf("done.\n");
                 return MAGMA_SUCCESS;
             }*/
             else {
                 printf("error: format not supported.\n");
                 return MAGMA_ERR_NOT_SUPPORTED;
             }
        }
         
         
    }
    // CPU case missing!     
    else{
        printf("error: CPU not yet supported.\n");
        return MAGMA_ERR_NOT_SUPPORTED;
    }
    return MAGMA_SUCCESS;
}



/*  -- MAGMA (version 1.5.0-beta2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2014

    Purpose
    =======

    For a given input matrix A and vectors x, y and scalars alpha, beta
    the wrapper determines the suitable SpMV computing
              y = alpha * ( A - lambda I ) * x + beta * y.  
    Arguments
    =========

    double alpha      scalar alpha
    magma_d_sparse_matrix A       sparse matrix A   
    double alpha      scalar lambda 
    magma_d_vector x              input vector x  
    magmaDoublComplex beta        scalar beta
    magma_d_vector y              input vector y      

    ========================================================================  */


magma_int_t
magma_d_spmv_shift( double alpha, 
                    magma_d_sparse_matrix A, 
                    double lambda,
                    magma_d_vector x, 
                    double beta, 
                    magma_int_t offset, 
                    magma_int_t blocksize, 
                    magma_index_t *add_rows, 
                    magma_d_vector y ){

    if( A.memory_location != x.memory_location 
                || x.memory_location != y.memory_location ){
    printf("error: linear algebra objects are not located in same memory!\n");
    printf("memory locations are: %d   %d   %d\n", 
                    A.memory_location, x.memory_location, y.memory_location );
    return MAGMA_ERR_INVALID_PTR;
    }
    // DEV case
    if( A.memory_location == Magma_DEV ){
         if( A.storage_type == Magma_CSR ){
             //printf("using CSR kernel for SpMV: ");
             magma_dgecsrmv_shift( MagmaNoTrans, A.num_rows, A.num_cols, 
                alpha, lambda, A.val, A.row, A.col, x.val, beta, offset, 
                blocksize, add_rows, y.val );
             //printf("done.\n");
             return MAGMA_SUCCESS;
         }
         else if( A.storage_type == Magma_ELLPACK ){
             //printf("using ELLPACK kernel for SpMV: ");
             magma_dgeellmv_shift( MagmaNoTrans, A.num_rows, A.num_cols, 
                A.max_nnz_row, alpha, lambda, A.val, A.col, x.val, beta, offset, 
                blocksize, add_rows, y.val );
             //printf("done.\n");
             return MAGMA_SUCCESS;
         }
         else if( A.storage_type == Magma_ELL ){
             //printf("using ELL kernel for SpMV: ");
             magma_dgeelltmv_shift( MagmaNoTrans, A.num_rows, A.num_cols, 
                A.max_nnz_row, alpha, lambda, A.val, A.col, x.val, beta, offset, 
                blocksize, add_rows, y.val );
             //printf("done.\n");
             return MAGMA_SUCCESS;
         }
         else {
             printf("error: format not supported.\n");
             return MAGMA_ERR_NOT_SUPPORTED;
         }
    }
    // CPU case missing!     
    else{
        printf("error: CPU not yet supported.\n");
        return MAGMA_ERR_NOT_SUPPORTED;
    }

}