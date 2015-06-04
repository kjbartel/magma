/*
    -- MAGMA (version 1.6.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date November 2014

       @author Hartwig Anzt 

       @generated from zcg_merge.cpp normal z -> c, Sat Nov 15 19:54:22 2014
*/

#include "common_magma.h"
#include "magmasparse.h"

#include <assert.h>

#define RTOLERANCE     lapackf77_slamch( "E" )
#define ATOLERANCE     lapackf77_slamch( "E" )


/**
    Purpose
    -------

    Solves a system of linear equations
       A * X = B
    where A is a complex Hermitian N-by-N positive definite matrix A.
    This is a GPU implementation of the Conjugate Gradient method in variant,
    where multiple operations are merged into one compute kernel.    

    Arguments
    ---------

    @param[in]
    A           magma_c_sparse_matrix
                input matrix A

    @param[in]
    b           magma_c_vector
                RHS b

    @param[in,out]
    x           magma_c_vector*
                solution approximation

    @param[in,out]
    solver_par  magma_c_solver_par*
                solver parameters
    @param[in]
    queue       magma_queue_t
                Queue to execute in.

    @ingroup magmasparse_chesv
    ********************************************************************/

extern "C" magma_int_t
magma_ccg_merge(
    magma_c_sparse_matrix A, magma_c_vector b, magma_c_vector *x,  
    magma_c_solver_par *solver_par,
    magma_queue_t queue )
{
    // set queue for old dense routines
    magma_queue_t orig_queue;
    magmablasGetKernelStream( &orig_queue );

    // prepare solver feedback
    solver_par->solver = Magma_CGMERGE;
    solver_par->numiter = 0;
    solver_par->info = MAGMA_SUCCESS; 

    // some useful variables
    magmaFloatComplex c_zero = MAGMA_C_ZERO, c_one = MAGMA_C_ONE;
    magma_int_t dofs = A.num_rows;

    // GPU stream
    magma_queue_t stream[2];
    magma_event_t event[1];
    magma_queue_create( &stream[0] );
    magma_queue_create( &stream[1] );
    magma_event_create( &event[0] );

    // GPU workspace
    magma_c_vector r, d, z;
    magma_c_vinit( &r, Magma_DEV, dofs, c_zero, queue );
    magma_c_vinit( &d, Magma_DEV, dofs, c_zero, queue );
    magma_c_vinit( &z, Magma_DEV, dofs, c_zero, queue );
    
    magmaFloatComplex *d1, *d2, *skp;
    d1 = NULL;
    d2 = NULL;
    skp = NULL;
    magma_int_t stat_dev = 0, stat_cpu = 0;
    stat_dev += magma_cmalloc( &d1, dofs*(1) );
    stat_dev += magma_cmalloc( &d2, dofs*(1) );
    // array for the parameters
    stat_dev += magma_cmalloc( &skp, 6 );       
    // skp = [alpha|beta|gamma|rho|tmp1|tmp2]
    if( stat_dev != 0 ){
        magma_free( d1 );
        magma_free( d2 );
        magma_free( skp );
        printf("error: memory allocation.\n");
        return MAGMA_ERR_DEVICE_ALLOC;
    }

    // solver variables
    magmaFloatComplex alpha, beta, gamma, rho, tmp1, *skp_h;
    float nom, nom0, r0, betanom, den;

    // solver setup
    magma_cscal( dofs, c_zero, x->dval, 1) ;                     // x = 0
    magma_ccopy( dofs, b.dval, 1, r.dval, 1 );                    // r = b
    magma_ccopy( dofs, b.dval, 1, d.dval, 1 );                    // d = b
    nom0 = betanom = magma_scnrm2( dofs, r.dval, 1 );               
    nom = nom0 * nom0;                                           // nom = r' * r
    magma_c_spmv( c_one, A, d, c_zero, z, queue );                      // z = A d
    den = MAGMA_C_REAL( magma_cdotc(dofs, d.dval, 1, z.dval, 1) ); // den = d'* z
    solver_par->init_res = nom0;
    
    // array on host for the parameters
    stat_cpu += magma_cmalloc_cpu( &skp_h, 6 );
    if( stat_cpu != 0 ){
        magma_free( d1 );
        magma_free( d2 );
        magma_free( skp );
        magma_free_cpu( skp_h );
        printf("error: memory allocation.\n");
        return MAGMA_ERR_HOST_ALLOC;
    }
    
    alpha = rho = gamma = tmp1 = c_one; 
    beta =  magma_cdotc(dofs, r.dval, 1, r.dval, 1);
    skp_h[0]=alpha; 
    skp_h[1]=beta; 
    skp_h[2]=gamma; 
    skp_h[3]=rho; 
    skp_h[4]=tmp1; 
    skp_h[5]=MAGMA_C_MAKE(nom, 0.0);

    magma_csetvector( 6, skp_h, 1, skp, 1 );
    
    if ( (r0 = nom * solver_par->epsilon) < ATOLERANCE ) 
        r0 = ATOLERANCE;
    if ( nom < r0 ) {
        magmablasSetKernelStream( orig_queue );
        return MAGMA_SUCCESS;
    }
    // check positive definite
    if (den <= 0.0) {
        printf("Operator A is not postive definite. (Ar,r) = %f\n", den);
        magmablasSetKernelStream( orig_queue );
        return MAGMA_NONSPD;
        solver_par->info = MAGMA_NONSPD;;
    }
    
    //Chronometry
    real_Double_t tempo1, tempo2;
    tempo1 = magma_sync_wtime( queue );
    if ( solver_par->verbose > 0 ) {
        solver_par->res_vec[0] = (real_Double_t) nom0;
        solver_par->timing[0] = 0.0;
    }
    
    // start iteration
    for( solver_par->numiter= 1; solver_par->numiter<solver_par->maxiter; 
                                                    solver_par->numiter++ ) {

        magmablasSetKernelStream(stream[0]);
        
        // computes SpMV and dot product
        magma_ccgmerge_spmv1(  A, d1, d2, d.dval, z.dval, skp, queue ); 
            
        // updates x, r, computes scalars and updates d
        magma_ccgmerge_xrbeta( dofs, d1, d2, x->dval, r.dval, d.dval, z.dval, skp, queue ); 

        // check stopping criterion (asynchronous copy)
        magma_cgetvector_async( 1 , skp+1, 1, 
                                                    skp_h+1, 1, stream[1] );
        betanom = sqrt(MAGMA_C_REAL(skp_h[1]));

        if ( solver_par->verbose > 0 ) {
            tempo2 = magma_sync_wtime( queue );
            if ( (solver_par->numiter)%solver_par->verbose==0 ) {
                solver_par->res_vec[(solver_par->numiter)/solver_par->verbose] 
                        = (real_Double_t) betanom;
                solver_par->timing[(solver_par->numiter)/solver_par->verbose] 
                        = (real_Double_t) tempo2-tempo1;
            }
        }

        if (  betanom  < r0 ) {
            break;
        }

    } 
    tempo2 = magma_sync_wtime( queue );
    solver_par->runtime = (real_Double_t) tempo2-tempo1;
    float residual;
    magma_cresidual( A, b, *x, &residual, queue );
    solver_par->iter_res = betanom;
    solver_par->final_res = residual;

    if ( solver_par->numiter < solver_par->maxiter) {
        solver_par->info = MAGMA_SUCCESS;
    } else if ( solver_par->init_res > solver_par->final_res ) {
        if ( solver_par->verbose > 0 ) {
            if ( (solver_par->numiter)%solver_par->verbose==0 ) {
                solver_par->res_vec[(solver_par->numiter)/solver_par->verbose] 
                        = (real_Double_t) betanom;
                solver_par->timing[(solver_par->numiter)/solver_par->verbose] 
                        = (real_Double_t) tempo2-tempo1;
            }
        }
        solver_par->info = MAGMA_SLOW_CONVERGENCE;
    }
    else {
        if ( solver_par->verbose > 0 ) {
            if ( (solver_par->numiter)%solver_par->verbose==0 ) {
                solver_par->res_vec[(solver_par->numiter)/solver_par->verbose] 
                        = (real_Double_t) betanom;
                solver_par->timing[(solver_par->numiter)/solver_par->verbose] 
                        = (real_Double_t) tempo2-tempo1;
            }
        }
        solver_par->info = MAGMA_DIVERGENCE;
    }
    magma_c_vfree(&r, queue );
    magma_c_vfree(&z, queue );
    magma_c_vfree(&d, queue );

    magma_free( d1 );
    magma_free( d2 );
    magma_free( skp );
    magma_free_cpu( skp_h );

    magmablasSetKernelStream( orig_queue );
    return MAGMA_SUCCESS;
}   /* magma_ccg_merge */


