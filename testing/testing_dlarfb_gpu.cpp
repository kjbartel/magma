/*
 *  -- MAGMA (version 1.2.0) --
 *     Univ. of Tennessee, Knoxville
 *     Univ. of California, Berkeley
 *     Univ. of Colorado, Denver
 *     May 2012
 *
 * @author Mark Gates
 * @generated d Tue May 15 18:18:23 2012
 *
 **/

// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cublas.h>
#include <assert.h>

#include <algorithm>  // std::swap

// includes, project
#include "magma.h"
#include "magma_lapack.h"
#include "testings.h"

/* ////////////////////////////////////////////////////////////////////////////
   -- Testing dlarfb_gpu
*/
int main( int argc, char** argv )
{
    TESTING_CUDA_INIT();
    
    double c_zero    = MAGMA_D_ZERO;
    double c_one     = MAGMA_D_ONE;
    double c_neg_one = MAGMA_D_NEG_ONE;
    magma_int_t ione =  1;
    
    printf( "\nUsage: %s -M m -N n -K k\n\n", argv[0] );

    magma_int_t m = 500;
    magma_int_t n = 300;
    magma_int_t k = 32;
    for( int i = 1; i < argc; i++ ) {
        if (strcmp("-M", argv[i]) == 0 and ++i < argc) {
            m = atoi( argv[i] );
        }
        else if (strcmp("-N", argv[i]) == 0 and ++i < argc) {
            n = atoi( argv[i] );
        }
        else if (strcmp("-K", argv[i]) == 0 and ++i < argc) {
            k = atoi( argv[i] );
        }
        else {
            printf( "invalid argument: %s\n", argv[i] );
            exit(1);
        }
    }
    if ( k <= 0 or k > m or k > n ) {
        printf( "requires 0 < k <= min(m,n)\n" );
        exit(1);
    }
    
    magma_int_t ldc = m;
    magma_int_t ldv = max(m,n);
    magma_int_t ldt = k;
    magma_int_t ldw = max(m,n);
    magma_int_t nv;
    ldc = ((ldc+31)/32)*32;
    ldv = ((ldv+31)/32)*32;
    ldt = ((ldt+31)/32)*32; 
    ldw = ((ldw+31)/32)*32;
    
    // Allocate memory for matrices
    double *C, *R, *V, *T, *W;
    TESTING_MALLOC( C, double, ldc*n );
    TESTING_MALLOC( R, double, ldc*n );
    TESTING_MALLOC( V, double, ldv*k );
    TESTING_MALLOC( T, double, ldt*k );
    TESTING_MALLOC( W, double, ldw*k );
    
    double *dC, *dV, *dT, *dW;
    TESTING_DEVALLOC( dC, double, ldc*n );
    TESTING_DEVALLOC( dV, double, ldv*k );
    TESTING_DEVALLOC( dT, double, ldt*k );
    TESTING_DEVALLOC( dW, double, ldw*k );
    
    magma_int_t size;
    magma_int_t idist    = 1;
    magma_int_t iseed[4] = { 1, 2, 3, 4 };
    double error, work[1];
    
    // test all combinations of input parameters
    const char* side[]   = { "Left",    "Right"     };
    const char* trans[]  = { "NoTrans", "ConjTrans" };
    const char* direct[] = { "Forward", "Backward"  };
    const char* storev[] = { "Colwise", "Rowwise"   };

    printf("    M     N     K  storev     side       direct     trans       ||R||_F / ||HC||_F\n");
    printf("==================================================================================\n");
    for( int iv = 0; iv < 2; ++iv ) {
    for( int is = 0; is < 2; ++is ) {
    for( int id = 0; id < 2; ++id ) {
    for( int it = 0; it < 2; ++it ) {
        //printf( "# ----------\n" );
        //printf( "# %-10s %-10s %-10s %-10s\n", storev[iv], side[is], direct[id], trans[it] );
        
        // C is full
        size = ldc*n;
        lapackf77_dlarnv( &idist, iseed, &size, C );
        //printf( "C=" );  magma_dprint( m, n, C, ldc );
        
        // V is ldv x nv. See larfb docs for description.
        ldv  = (*side[is] == 'L' ? m : n);
        nv   = k;
        size = ldv*nv;
        lapackf77_dlarnv( &idist, iseed, &size, V );
        if ( *storev[iv] == 'C' ) {
            if ( *direct[id] == 'F' ) {
                lapackf77_dlaset( "Upper", &k, &k, &c_zero, &c_one, V, &ldv );
            }
            else {
                lapackf77_dlaset( "Lower", &k, &k, &c_zero, &c_one, &V[(ldv-k)], &ldv );
            }
        }
        else {
            // rowwise, swap V's dimensions
            std::swap( ldv, nv );
            if ( *direct[id] == 'F' ) {
                lapackf77_dlaset( "Lower", &k, &k, &c_zero, &c_one, V, &ldv );
            }
            else {
                lapackf77_dlaset( "Upper", &k, &k, &c_zero, &c_one, &V[(nv-k)*ldv], &ldv );
            }
        }
        //printf( "# ldv %d, nv %d\n", ldv, nv );
        //printf( "V=" );  magma_dprint( ldv, nv, V, ldv );
        
        // T is upper triangular for forward, and lower triangular for backward
        magma_int_t k1 = k-1;
        size = ldt*k;
        lapackf77_dlarnv( &idist, iseed, &size, T );
        if ( *direct[id] == 'F' ) {
            lapackf77_dlaset( "Lower", &k1, &k1, &c_zero, &c_zero, &T[1], &ldt );
        }
        else {
            lapackf77_dlaset( "Upper", &k1, &k1, &c_zero, &c_zero, &T[1*ldt], &ldt );
        }
        //printf( "T=" );  magma_dprint( k, k, T, ldt );
        
        magma_dsetmatrix( m, n, C, ldc, dC, ldc );
        magma_dsetmatrix( ldv, nv, V, ldv, dV, ldv );
        magma_dsetmatrix( k, k, T, ldt, dT, ldt );
        
        lapackf77_dlarfb( side[is], trans[it], direct[id], storev[iv],
                          &m, &n, &k,
                          V, &ldv, T, &ldt, C, &ldc, W, &ldw );
        //printf( "HC=" );  magma_dprint( m, n, C, ldc );
        
        magma_dlarfb_gpu( *side[is], *trans[it], *direct[id], *storev[iv],
                          m, n, k,
                          dV, ldv, dT, ldt, dC, ldc, dW, ldw );
        magma_dgetmatrix( m, n, dC, ldc, R, ldc );
        //printf( "dHC=" );  magma_dprint( m, n, R, ldc );
        
        // compute relative error |HC_magma - HC_lapack| / |HC_lapack|
        error = lapackf77_dlange( "Fro", &m, &n, C, &ldc, work );
        size = ldc*n;
        blasf77_daxpy( &size, &c_neg_one, C, &ione, R, &ione );
        error = lapackf77_dlange( "Fro", &m, &n, R, &ldc, work ) / error;
        printf( "%5d %5d %5d  %-10s %-10s %-10s %-10s  %8.2e\n",
                m, n, k, storev[iv], side[is], direct[id], trans[it], error );
    }}}}
    
    // Memory clean up
    TESTING_FREE( C );
    TESTING_FREE( R );
    TESTING_FREE( V );
    TESTING_FREE( T );
    TESTING_FREE( W );
    
    TESTING_DEVFREE( dC );
    TESTING_DEVFREE( dV );
    TESTING_DEVFREE( dT );
    TESTING_DEVFREE( dW );
    
    // Shutdown
    TESTING_CUDA_FINALIZE();
    return 0;
}
