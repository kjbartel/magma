/*
    -- MAGMA (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @generated from magma_z_solver_wrapper.cpp normal z -> d, Sun May  3 11:22:59 2015
       @author Hartwig Anzt

*/
#include "common_magmasparse.h"




/**
    Purpose
    -------

    ALlows the user to choose a solver.

    Arguments
    ---------

    @param[in]
    A           magma_d_matrix
                sparse matrix A

    @param[in]
    b           magma_d_matrix
                input vector b

    @param[in]
    x           magma_d_matrix*
                output vector x

    @param[in]
    zopts     magma_dopts
              options for solver and preconditioner
    @param[in]
    queue       magma_queue_t
                Queue to execute in.

    @ingroup magmasparse_daux
    ********************************************************************/

extern "C" magma_int_t
magma_d_solver(
    magma_d_matrix A, magma_d_matrix b,
    magma_d_matrix *x, magma_dopts *zopts,
    magma_queue_t queue )
{
    magma_int_t info = 0;
    
    // make sure RHS is a dense matrix
    if ( b.storage_type != Magma_DENSE ) {
        printf( "error: sparse RHS not yet supported.\n" );
        return MAGMA_ERR_NOT_SUPPORTED;
    }
    if( b.num_cols == 1 ){
    // preconditioner
        if ( zopts->solver_par.solver != Magma_ITERREF ) {
            int stat = magma_d_precondsetup( A, b, &zopts->precond_par, queue );
            if (  stat != MAGMA_SUCCESS ){
                printf("error: bad preconditioner.\n");
                return MAGMA_ERR_BADPRECOND; 
            }
        }
        switch( zopts->solver_par.solver ) {
            case  Magma_CG:
                    CHECK( magma_dcg_res( A, b, x, &zopts->solver_par, queue )); break;
            case  Magma_CGMERGE:
                    CHECK( magma_dcg_merge( A, b, x, &zopts->solver_par, queue )); break;
            case  Magma_PCG:
                    CHECK( magma_dpcg( A, b, x, &zopts->solver_par, &zopts->precond_par, queue )); break;
            case  Magma_BICGSTAB:
                    CHECK( magma_dbicgstab( A, b, x, &zopts->solver_par, queue )); break;
            case  Magma_BICGSTABMERGE: 
                    CHECK( magma_dbicgstab_merge( A, b, x, &zopts->solver_par, queue )); break;
            case  Magma_PBICGSTAB:
                    CHECK( magma_dpbicgstab( A, b, x, &zopts->solver_par, &zopts->precond_par, queue )); break;
            case  Magma_GMRES:
                    CHECK( magma_dfgmres( A, b, x, &zopts->solver_par, &zopts->precond_par, queue )); break;
            case  Magma_PGMRES:
                    CHECK( magma_dfgmres( A, b, x, &zopts->solver_par, &zopts->precond_par, queue )); break;
            case  Magma_LOBPCG:
                    CHECK( magma_dlobpcg( A, &zopts->solver_par, &zopts->precond_par, queue )); break;
            case  Magma_ITERREF:
                    CHECK( magma_diterref( A, b, x, &zopts->solver_par, &zopts->precond_par, queue )); break;
            case  Magma_JACOBI:
                    CHECK( magma_djacobi( A, b, x, &zopts->solver_par, queue )); break;
            case  Magma_BAITER:
                    CHECK( magma_dbaiter( A, b, x, &zopts->solver_par, queue ) ); break;
            default:
                    printf("error: solver class not supported.\n"); break;
        }
    }
    else{
  // preconditioner
        if ( zopts->solver_par.solver != Magma_ITERREF ) {
            int stat = magma_d_precondsetup( A, b, &zopts->precond_par, queue );
            if (  stat != MAGMA_SUCCESS ){
                printf("error: bad preconditioner.\n");
                return MAGMA_ERR_BADPRECOND; 
            }
        }
        switch( zopts->solver_par.solver ) {
            case  Magma_CG:
                    CHECK( magma_dbpcg( A, b, x, &zopts->solver_par, &zopts->precond_par, queue )); break;
            case  Magma_PCG:
                    CHECK( magma_dbpcg( A, b, x, &zopts->solver_par, &zopts->precond_par, queue )); break;
            case  Magma_LOBPCG:
                    CHECK( magma_dlobpcg( A, &zopts->solver_par, &zopts->precond_par, queue )); break;
            default:
                    printf("error: only 1 RHS supported for this solver class.\n"); break;
        }
    }
cleanup:
    return info; 
}


