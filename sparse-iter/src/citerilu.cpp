/*
    -- MAGMA (version 1.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date

       @author Hartwig Anzt

       @generated from ziterilu.cpp normal z -> c, Tue May  5 14:27:23 2015
*/
#include "common_magmasparse.h"

#define PRECISION_c


/**
    Purpose
    -------

    Prepares the ILU preconditioner via the iterative ILU iteration.

    Arguments
    ---------

    @param[in]
    A           magma_c_matrix
                input matrix A

    @param[in][out]
    precond     magma_c_preconditioner*
                preconditioner parameters
                
    @param[in]
    queue       magma_queue_t
                Queue to execute in.

    @ingroup magmasparse_cgepr
    ********************************************************************/
extern "C"
magma_int_t
magma_citerilusetup(
    magma_c_matrix A,
    magma_c_matrix b,
    magma_c_preconditioner *precond,
    magma_queue_t queue )
{
    magma_int_t info = 0;

    cusparseHandle_t cusparseHandle=NULL;
    cusparseMatDescr_t descrL=NULL;
    cusparseMatDescr_t descrU=NULL;

    magma_c_matrix hAh={Magma_CSR}, hA={Magma_CSR}, hL={Magma_CSR}, hU={Magma_CSR},
    hAcopy={Magma_CSR}, hAL={Magma_CSR}, hAU={Magma_CSR}, hAUt={Magma_CSR},
    hUT={Magma_CSR}, hAtmp={Magma_CSR}, hACSRCOO={Magma_CSR}, dAinitguess={Magma_CSR},
    dL={Magma_CSR}, dU={Magma_CSR}, DL={Magma_CSR}, RL={Magma_CSR}, DU={Magma_CSR}, RU={Magma_CSR};

    // copy original matrix as CSRCOO to device
    CHECK( magma_cmtransfer(A, &hAh, A.memory_location, Magma_CPU, queue ));
    CHECK( magma_cmconvert( hAh, &hA, hAh.storage_type, Magma_CSR , queue ));
    magma_cmfree(&hAh, queue );

    CHECK( magma_cmtransfer( hA, &hAcopy, Magma_CPU, Magma_CPU , queue ));

    // in case using fill-in
    CHECK( magma_csymbilu( &hAcopy, precond->levels, &hAL, &hAUt,  queue ));
    // add a unit diagonal to L for the algorithm
    CHECK( magma_cmLdiagadd( &hAL , queue ));
    // transpose U for the algorithm
    CHECK( magma_c_cucsrtranspose(  hAUt, &hAU , queue ));
    magma_cmfree( &hAUt , queue );

    // ---------------- initial guess ------------------- //
    CHECK( magma_cmconvert( hAcopy, &hACSRCOO, Magma_CSR, Magma_CSRCOO , queue ));
    CHECK( magma_cmtransfer( hACSRCOO, &dAinitguess, Magma_CPU, Magma_DEV , queue ));
    magma_cmfree(&hACSRCOO, queue );
    magma_cmfree(&hAcopy, queue );

    // transfer the factor L and U
    CHECK( magma_cmtransfer( hAL, &dL, Magma_CPU, Magma_DEV , queue ));
    CHECK( magma_cmtransfer( hAU, &dU, Magma_CPU, Magma_DEV , queue ));
    magma_cmfree(&hAL, queue );
    magma_cmfree(&hAU, queue );

    for(int i=0; i<precond->sweeps; i++){
        CHECK( magma_citerilu_csr( dAinitguess, dL, dU , queue ));
    }

    CHECK( magma_cmtransfer( dL, &hL, Magma_DEV, Magma_CPU , queue ));
    CHECK( magma_cmtransfer( dU, &hU, Magma_DEV, Magma_CPU , queue ));
    CHECK( magma_c_cucsrtranspose(  hU, &hUT , queue ));

    magma_cmfree(&dL, queue );
    magma_cmfree(&dU, queue );
    magma_cmfree(&hU, queue );
    CHECK( magma_cmlumerge( hL, hUT, &hAtmp, queue ));

    magma_cmfree(&hL, queue );
    magma_cmfree(&hUT, queue );

    CHECK( magma_cmtransfer( hAtmp, &precond->M, Magma_CPU, Magma_DEV , queue ));

    hAL.diagorder_type = Magma_UNITY;
    CHECK( magma_cmconvert(hAtmp, &hAL, Magma_CSR, Magma_CSRL, queue ));
    hAL.storage_type = Magma_CSR;
    CHECK( magma_cmconvert(hAtmp, &hAU, Magma_CSR, Magma_CSRU, queue ));
    hAU.storage_type = Magma_CSR;

    magma_cmfree(&hAtmp, queue );

    CHECK( magma_ccsrsplit( 256, hAL, &DL, &RL , queue ));
    CHECK( magma_ccsrsplit( 256, hAU, &DU, &RU , queue ));

    CHECK( magma_cmtransfer( DL, &precond->LD, Magma_CPU, Magma_DEV , queue ));
    CHECK( magma_cmtransfer( DU, &precond->UD, Magma_CPU, Magma_DEV , queue ));

    // for cusparse uncomment this
    CHECK( magma_cmtransfer( hAL, &precond->L, Magma_CPU, Magma_DEV , queue ));
    CHECK( magma_cmtransfer( hAU, &precond->U, Magma_CPU, Magma_DEV , queue ));
    
/*

    //-- for ba-solve uncomment this

    if( RL.nnz != 0 )
        CHECK( magma_cmtransfer( RL, &precond->L, Magma_CPU, Magma_DEV , queue ));
    else{
        precond->L.nnz = 0;
        precond->L.val = NULL;
        precond->L.col = NULL;
        precond->L.row = NULL;
        precond->L.blockinfo = NULL;
    }

    if( RU.nnz != 0 )
        CHECK( magma_cmtransfer( RU, &precond->U, Magma_CPU, Magma_DEV , queue ));
    else{
        precond->U.nnz = 0;
        precond->L.val = NULL;
        precond->L.col = NULL;
        precond->L.row = NULL;
        precond->L.blockinfo = NULL;
    }

    //-- for ba-solve uncomment this
*/

        // extract the diagonal of L into precond->d
    CHECK( magma_cjacobisetup_diagscal( precond->L, &precond->d, queue ));
    CHECK( magma_cvinit( &precond->work1, Magma_DEV, hA.num_rows, 1, MAGMA_C_ZERO, queue ));
    
    // extract the diagonal of U into precond->d2
    CHECK( magma_cjacobisetup_diagscal( precond->U, &precond->d2, queue ));
    CHECK( magma_cvinit( &precond->work2, Magma_DEV, hA.num_rows, 1, MAGMA_C_ZERO, queue ));

    magma_cmfree(&hAL, queue );
    magma_cmfree(&hAU, queue );
    magma_cmfree(&DL, queue );
    magma_cmfree(&RL, queue );
    magma_cmfree(&DU, queue );
    magma_cmfree(&RU, queue );

    // CUSPARSE context //
    CHECK_CUSPARSE( cusparseCreate( &cusparseHandle ));
    CHECK_CUSPARSE( cusparseCreateMatDescr( &descrL ));
    CHECK_CUSPARSE( cusparseSetMatType( descrL, CUSPARSE_MATRIX_TYPE_TRIANGULAR ));
    CHECK_CUSPARSE( cusparseSetMatDiagType( descrL, CUSPARSE_DIAG_TYPE_UNIT ));
    CHECK_CUSPARSE( cusparseSetMatIndexBase( descrL, CUSPARSE_INDEX_BASE_ZERO ));
    CHECK_CUSPARSE( cusparseSetMatFillMode( descrL, CUSPARSE_FILL_MODE_LOWER ));
    CHECK_CUSPARSE( cusparseCreateSolveAnalysisInfo( &precond->cuinfoL ));
    CHECK_CUSPARSE( cusparseCcsrsv_analysis( cusparseHandle,
        CUSPARSE_OPERATION_NON_TRANSPOSE, precond->L.num_rows,
        precond->L.nnz, descrL,
        precond->L.val, precond->L.row, precond->L.col, precond->cuinfoL ));
    CHECK_CUSPARSE( cusparseCreateMatDescr( &descrU ));
    CHECK_CUSPARSE( cusparseSetMatType( descrU, CUSPARSE_MATRIX_TYPE_TRIANGULAR ));
    CHECK_CUSPARSE( cusparseSetMatDiagType( descrU, CUSPARSE_DIAG_TYPE_NON_UNIT ));
    CHECK_CUSPARSE( cusparseSetMatIndexBase( descrU, CUSPARSE_INDEX_BASE_ZERO ));
    CHECK_CUSPARSE( cusparseSetMatFillMode( descrU, CUSPARSE_FILL_MODE_UPPER ));
    CHECK_CUSPARSE( cusparseCreateSolveAnalysisInfo( &precond->cuinfoU ));
    CHECK_CUSPARSE( cusparseCcsrsv_analysis( cusparseHandle,
        CUSPARSE_OPERATION_NON_TRANSPOSE, precond->U.num_rows,
        precond->U.nnz, descrU,
        precond->U.val, precond->U.row, precond->U.col, precond->cuinfoU ));


cleanup:
    cusparseDestroy( cusparseHandle );
    cusparseDestroyMatDescr( descrL );
    cusparseDestroyMatDescr( descrU );
    cusparseHandle=NULL;
    descrL=NULL;
    descrU=NULL;
    magma_cmfree( &hAh, queue );
    magma_cmfree( &hA, queue );
    magma_cmfree( &hL, queue );
    magma_cmfree( &hU, queue );
    magma_cmfree( &hAcopy, queue );
    magma_cmfree( &hAL, queue );
    magma_cmfree( &hAU, queue );
    magma_cmfree( &hAUt, queue );
    magma_cmfree( &hUT, queue );
    magma_cmfree( &hAtmp, queue );
    magma_cmfree( &hACSRCOO, queue );
    magma_cmfree( &dAinitguess, queue );
    magma_cmfree( &dL, queue );
    magma_cmfree( &dU, queue );
    magma_cmfree( &DL, queue );
    magma_cmfree( &DU, queue );
    magma_cmfree( &RL, queue );
    magma_cmfree( &RU, queue );

    return info;
}






/**
    Purpose
    -------

    Prepares the IC preconditioner via the iterative IC iteration.

    Arguments
    ---------

    @param[in]
    A           magma_c_matrix
                input matrix A

    @param[in][out]
    precond     magma_c_preconditioner*
                preconditioner parameters
                
    @param[in]
    queue       magma_queue_t
                Queue to execute in.

    @ingroup magmasparse_chepr
    ********************************************************************/
extern "C"
magma_int_t
magma_citericsetup(
    magma_c_matrix A,
    magma_c_matrix b,
    magma_c_preconditioner *precond,
    magma_queue_t queue )
{
    magma_int_t info = 0;
    
    cusparseHandle_t cusparseHandle=NULL;
    cusparseMatDescr_t descrL=NULL;
    cusparseMatDescr_t descrU=NULL;

    magma_c_matrix hAh={Magma_CSR}, hA={Magma_CSR}, hAtmp={Magma_CSR},
    hAL={Magma_CSR}, hAUt={Magma_CSR}, hALt={Magma_CSR}, hM={Magma_CSR},
    hACSRCOO={Magma_CSR}, dAinitguess={Magma_CSR}, dL={Magma_CSR};
    magma_c_matrix d_h={Magma_CSR};


    // copy original matrix as CSRCOO to device
    CHECK( magma_cmtransfer(A, &hAh, A.memory_location, Magma_CPU, queue ));
    CHECK( magma_cmconvert( hAh, &hA, hAh.storage_type, Magma_CSR , queue ));
    magma_cmfree(&hAh, queue );

    // in case using fill-in
    CHECK( magma_csymbilu( &hA, precond->levels, &hAL, &hAUt , queue ));

    // need only lower triangular
    magma_cmfree(&hAUt, queue );
    magma_cmfree(&hAL, queue );
    CHECK( magma_cmconvert( hA, &hAtmp, Magma_CSR, Magma_CSRL , queue ));
    magma_cmfree(&hA, queue );

    // ---------------- initial guess ------------------- //
    CHECK( magma_cmconvert( hAtmp, &hACSRCOO, Magma_CSR, Magma_CSRCOO , queue ));
    //int blocksize = 1;
    //magma_cmreorder( hACSRCOO, n, blocksize, blocksize, blocksize, &hAinitguess , queue );
    CHECK( magma_cmtransfer( hACSRCOO, &dAinitguess, Magma_CPU, Magma_DEV , queue ));
    magma_cmfree(&hACSRCOO, queue );
    CHECK( magma_cmtransfer( hAtmp, &dL, Magma_CPU, Magma_DEV , queue ));
    magma_cmfree(&hAtmp, queue );

    for(int i=0; i<precond->sweeps; i++){
        CHECK( magma_citeric_csr( dAinitguess, dL , queue ));
    }
    CHECK( magma_cmtransfer( dL, &hAL, Magma_DEV, Magma_CPU , queue ));
    magma_cmfree(&dL, queue );
    magma_cmfree(&dAinitguess, queue );


    // for CUSPARSE
    CHECK( magma_cmtransfer( hAL, &precond->M, Magma_CPU, Magma_DEV , queue ));

    // Jacobi setup
    CHECK( magma_cjacobisetup_matrix( precond->M, &precond->L, &precond->d , queue ));

    // for Jacobi, we also need U
    CHECK( magma_c_cucsrtranspose(   hAL, &hALt , queue ));
    CHECK( magma_cjacobisetup_matrix( hALt, &hM, &d_h , queue ));

    CHECK( magma_cmtransfer( hM, &precond->U, Magma_CPU, Magma_DEV , queue ));

    magma_cmfree(&hM, queue );

    magma_cmfree(&d_h, queue );


        // copy the matrix to precond->L and (transposed) to precond->U
    CHECK( magma_cmtransfer(precond->M, &(precond->L), Magma_DEV, Magma_DEV, queue ));
    CHECK( magma_cmtranspose( precond->L, &(precond->U), queue ));

    // extract the diagonal of L into precond->d
    CHECK( magma_cjacobisetup_diagscal( precond->L, &precond->d, queue ));
    CHECK( magma_cvinit( &precond->work1, Magma_DEV, hAL.num_rows, 1, MAGMA_C_ZERO, queue ));

    // extract the diagonal of U into precond->d2
    CHECK( magma_cjacobisetup_diagscal( precond->U, &precond->d2, queue ));
    CHECK( magma_cvinit( &precond->work2, Magma_DEV, hAL.num_rows, 1, MAGMA_C_ZERO, queue ));


    magma_cmfree(&hAL, queue );
    magma_cmfree(&hALt, queue );


    // CUSPARSE context //
    CHECK_CUSPARSE( cusparseCreate( &cusparseHandle ));
    CHECK_CUSPARSE( cusparseCreateMatDescr( &descrL ));
    CHECK_CUSPARSE( cusparseSetMatType( descrL, CUSPARSE_MATRIX_TYPE_TRIANGULAR ));
    CHECK_CUSPARSE( cusparseSetMatDiagType( descrL, CUSPARSE_DIAG_TYPE_NON_UNIT ));
    CHECK_CUSPARSE( cusparseSetMatIndexBase( descrL, CUSPARSE_INDEX_BASE_ZERO ));
    CHECK_CUSPARSE( cusparseSetMatFillMode( descrL, CUSPARSE_FILL_MODE_LOWER ));
    CHECK_CUSPARSE( cusparseCreateSolveAnalysisInfo( &precond->cuinfoL ));
    CHECK_CUSPARSE( cusparseCcsrsv_analysis( cusparseHandle,
        CUSPARSE_OPERATION_NON_TRANSPOSE, precond->M.num_rows,
        precond->M.nnz, descrL,
        precond->M.val, precond->M.row, precond->M.col, precond->cuinfoL ));
    CHECK_CUSPARSE( cusparseCreateMatDescr( &descrU ));
    CHECK_CUSPARSE( cusparseSetMatType( descrU, CUSPARSE_MATRIX_TYPE_TRIANGULAR ));
    CHECK_CUSPARSE( cusparseSetMatDiagType( descrU, CUSPARSE_DIAG_TYPE_NON_UNIT ));
    CHECK_CUSPARSE( cusparseSetMatIndexBase( descrU, CUSPARSE_INDEX_BASE_ZERO ));
    CHECK_CUSPARSE( cusparseSetMatFillMode( descrU, CUSPARSE_FILL_MODE_LOWER ));
    CHECK_CUSPARSE( cusparseCreateSolveAnalysisInfo( &precond->cuinfoU ));
    CHECK_CUSPARSE( cusparseCcsrsv_analysis( cusparseHandle,
        CUSPARSE_OPERATION_TRANSPOSE, precond->M.num_rows,
        precond->M.nnz, descrU,
        precond->M.val, precond->M.row, precond->M.col, precond->cuinfoU ));

    
    cleanup:
    cusparseDestroy( cusparseHandle );
    cusparseDestroyMatDescr( descrL );
    cusparseDestroyMatDescr( descrU );
    cusparseHandle=NULL;
    descrL=NULL;
    descrU=NULL;    
    magma_cmfree( &hAh, queue );
    magma_cmfree( &hA, queue );
    magma_cmfree( &hAtmp, queue );
    magma_cmfree( &hAL, queue );
    magma_cmfree( &hAUt, queue );
    magma_cmfree( &hALt, queue );
    magma_cmfree( &hM, queue );
    magma_cmfree( &hACSRCOO, queue );
    magma_cmfree( &dAinitguess, queue );
    magma_cmfree( &dL, queue );
    magma_cmfree( &d_h, queue );
    
    return info;
}


/**
    Purpose
    -------

    Updates an existing preconditioner via additional iterative IC sweeps for
    previous factorization initial guess (PFIG).
    See  Anzt et al., Parallel Computing, 2015.

    Arguments
    ---------

    @param[in]
    precond         magma_c_preconditioner*
                    preconditioner parameters

    @param[in]
    magma_int_t     number of updates
    
    @param[in]
    queue       magma_queue_t
                Queue to execute in.
                
    @ingroup magmasparse_chepr
    ********************************************************************/
extern "C"
magma_int_t
magma_citericupdate(
    magma_c_matrix A,
    magma_c_preconditioner *precond,
    magma_int_t updates,
    magma_queue_t queue )
{
    magma_int_t info = 0;

    magma_c_matrix hALt={Magma_CSR};
    magma_c_matrix d_h={Magma_CSR};


    // copy original matrix as CSRCOO to device

    for(int i=0; i<updates; i++){
        CHECK( magma_citeric_csr( A, precond->M , queue ));
    }
    //magma_cmtransfer( precond->M, &precond->M, Magma_DEV, Magma_DEV , queue );
    magma_cmfree(&precond->L, queue );
    magma_cmfree(&precond->U, queue );
    magma_cmfree( &precond->d , queue );


    // Jacobi setup
    CHECK( magma_cjacobisetup_matrix( precond->M, &precond->L, &precond->d , queue ));

    // for Jacobi, we also need U
    CHECK( magma_c_cucsrtranspose(   precond->M, &hALt , queue ));
    CHECK( magma_cjacobisetup_matrix( hALt, &precond->U, &d_h , queue ));


cleanup:
    magma_cmfree(&d_h, queue );
    magma_cmfree(&hALt, queue );
    
    return info;
}


/**
    Purpose
    -------

    Performs the left triangular solves using the IC preconditioner via Jacobi.

    Arguments
    ---------

    @param[in]
    b           magma_c_matrix
                RHS

    @param[out]
    x           magma_c_matrix*
                vector to precondition

    @param[in]
    precond     magma_c_preconditioner*
                preconditioner parameters
                
    @param[in]
    queue       magma_queue_t
                Queue to execute in.

    @ingroup magmasparse_cgepr
    ********************************************************************/
extern "C"
magma_int_t
magma_capplyiteric_l(
    magma_c_matrix b,
    magma_c_matrix *x,
    magma_c_preconditioner *precond,
    magma_queue_t queue )
{
    magma_int_t info = 0;
    
    magma_int_t dofs = precond->L.num_rows;
    magma_c_solver_par jacobiiter_par;
    jacobiiter_par.maxiter = precond->maxiter;

    // compute c = D^{-1}b and copy c as initial guess to x
    CHECK( magma_cjacobisetup_vector_gpu( dofs, b, precond->d,
                                                precond->work1, x, queue ));
    // Jacobi iterator
    CHECK( magma_cjacobiiter_precond( precond->L, x, &jacobiiter_par, precond , queue ));

cleanup:
    return info;
}


/**
    Purpose
    -------

    Performs the right triangular solves using the IC preconditioner via Jacobi.

    Arguments
    ---------

    @param[in]
    b           magma_c_matrix
                RHS

    @param[out]
    x           magma_c_matrix*
                vector to precondition

    @param[in]
    precond     magma_c_preconditioner*
                preconditioner parameters
                
    @param[in]
    queue       magma_queue_t
                Queue to execute in.

    @ingroup magmasparse_cgepr
    ********************************************************************/
extern "C"
magma_int_t
magma_capplyiteric_r(
    magma_c_matrix b,
    magma_c_matrix *x,
    magma_c_preconditioner *precond,
    magma_queue_t queue )
{
    magma_int_t info = 0;

    magma_int_t dofs = precond->U.num_rows;
    magma_c_solver_par jacobiiter_par;
    jacobiiter_par.maxiter = precond->maxiter;

    // compute c = D^{-1}b and copy c as initial guess to x
    CHECK( magma_cjacobisetup_vector_gpu( dofs, b, precond->d,
                                                precond->work1, x, queue ));

    // Jacobi iterator
    CHECK( magma_cjacobiiter_precond( precond->U, x, &jacobiiter_par, precond , queue ));
    
cleanup:
    return info;
}
