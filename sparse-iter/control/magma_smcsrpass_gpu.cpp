/*
    -- MAGMA (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from magma_zmcsrpass_gpu.cpp normal z -> s, Fri Jan 30 19:00:32 2015
       @author Hartwig Anzt
*/

//  in this file, many routines are taken from 
//  the IO functions provided by MatrixMarket

#include <fstream>
#include <stdlib.h>
#include <string>
#include <sstream>
#include <iostream>
#include <ostream>
#include <assert.h>
#include <stdio.h>

#include "magmasparse_s.h"
#include "magma.h"
#include "mmio.h"


using namespace std;




/**
    Purpose
    -------

    Passes a CSR matrix to MAGMA (located on DEV).

    Arguments
    ---------

    @param[in]
    m           magma_int_t 
                number of rows

    @param[in]
    n           magma_int_t 
                number of columns

    @param[in]
    row         magmaIndex_ptr 
                row pointer

    @param[in]
    col         magmaIndex_ptr 
                column indices

    @param[in]
    val         magmaFloat_ptr 
                array containing matrix entries

    @param[out]
    A           magma_s_sparse_matrix*
                matrix in magma sparse matrix format
    @param[in]
    queue       magma_queue_t
                Queue to execute in.

    @ingroup magmasparse_saux
    ********************************************************************/

extern "C"
magma_int_t
magma_scsrset_gpu(
    magma_int_t m, 
    magma_int_t n, 
    magmaIndex_ptr row, 
    magmaIndex_ptr col, 
    magmaFloat_ptr val,
    magma_s_sparse_matrix *A,
    magma_queue_t queue )
{
    A->num_rows = m;
    A->num_cols = n;
    magma_index_t nnz;
    magma_index_getvector( 1, row+m, 1, &nnz, 1 );
    A->nnz = (magma_int_t) nnz;
    A->storage_type = Magma_CSR;
    A->memory_location = Magma_DEV;
    A->dval = val;
    A->dcol = col;
    A->drow = row;

    return MAGMA_SUCCESS;
}


/**
    Purpose
    -------

    Passes a MAGMA matrix to CSR structure (located on DEV).

    Arguments
    ---------

    @param[in]
    A           magma_s_sparse_matrix
                magma sparse matrix in CSR format

    @param[out]
    m           magma_int_t 
                number of rows

    @param[out]
    n           magma_int_t 
                number of columns

    @param[out]
    row         magmaIndex_ptr 
                row pointer

    @param[out]
    col         magmaIndex_ptr 
                column indices

    @param[out]
    val         magmaFloat_ptr 
                array containing matrix entries

    @param[in]
    queue       magma_queue_t
                Queue to execute in.

    @ingroup magmasparse_saux
    ********************************************************************/

extern "C"
magma_int_t
magma_scsrget_gpu(
    magma_s_sparse_matrix A,
    magma_int_t *m, 
    magma_int_t *n, 
    magmaIndex_ptr *row, 
    magmaIndex_ptr *col, 
    magmaFloat_ptr *val,
    magma_queue_t queue )
{
    if ( A.memory_location == Magma_DEV && A.storage_type == Magma_CSR ) {

        *m = A.num_rows;
        *n = A.num_cols;
        *val = A.dval;
        *col = A.dcol;
        *row = A.drow;
    } else {
        magma_s_sparse_matrix A_DEV, A_CSR;
        magma_s_mconvert( A, &A_CSR, A.storage_type, Magma_CSR, queue ); 
        magma_s_mtransfer( A_CSR, &A_DEV, A.memory_location, Magma_DEV, queue ); 
        magma_scsrget_gpu( A_DEV, m, n, row, col, val, queue );
        magma_s_mfree( &A_CSR, queue );
        magma_s_mfree( &A_DEV, queue );
    }
    return MAGMA_SUCCESS;
}


