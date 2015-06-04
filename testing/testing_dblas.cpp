/*
 *  -- MAGMA (version 1.2.1) --
 *     Univ. of Tennessee, Knoxville
 *     Univ. of California, Berkeley
 *     Univ. of Colorado, Denver
 *     June 2012
 *
 * @author Mark Gates
 * @generated d Thu Jun 28 12:31:51 2012
 *
 **/
#include <stdlib.h>
#include <stdio.h>

// make sure that asserts are enabled
#ifdef NDEBUG
#undef NDEBUG
#endif
#include <assert.h>

// includes, project
#include "magma.h"
#include "magma_lapack.h"
#include "testings.h"

#define A(i,j)  &A[  (i) + (j)*ld ]
#define dA(i,j) &dA[ (i) + (j)*ld ]
#define C2(i,j) &C2[ (i) + (j)*ld ]

int main( int argc, char** argv )
{
    TESTING_CUDA_INIT();
    
    double c_neg_one = MAGMA_D_NEG_ONE;
    magma_int_t ione = 1;
    const char trans[] = { 'N', 'C', 'T' };
    const char uplo[]  = { 'L', 'U' };
    const char diag[]  = { 'U', 'N' };
    const char side[]  = { 'L', 'R' };
    
    double  *A,  *B,  *C,   *C2;
    double *dA, *dB, *dC1, *dC2;
    double alpha = MAGMA_D_MAKE( 0.5, 0.1 );
    double beta  = MAGMA_D_MAKE( 0.7, 0.2 );
    double dalpha = 0.6;
    double dbeta  = 0.8;
    double work[1], error;
    magma_int_t m = 50;
    magma_int_t n = 35;
    magma_int_t k = 40;
    magma_int_t ISEED[4] = { 0, 1, 2, 3 };
    magma_int_t size, maxn, ld, info;
    magma_int_t *piv;
    magma_err_t err;
    
    printf( "\n" );
    
    // allocate matrices
    // over-allocate so they can be any combination of {m,n,k} x {m,n,k}.
    maxn = max( max( m, n ), k );
    ld = maxn;
    size = maxn*maxn;
    piv = (magma_int_t*) malloc( maxn * sizeof(magma_int_t) );
    assert( piv != NULL );
    err = magma_dmalloc_pinned( &A , size );  assert( err == 0 );
    err = magma_dmalloc_pinned( &B , size );  assert( err == 0 );
    err = magma_dmalloc_pinned( &C , size );  assert( err == 0 );
    err = magma_dmalloc_pinned( &C2, size );  assert( err == 0 );
    err = magma_dmalloc( &dA,  size );      assert( err == 0 );
    err = magma_dmalloc( &dB,  size );      assert( err == 0 );
    err = magma_dmalloc( &dC1, size );      assert( err == 0 );
    err = magma_dmalloc( &dC2, size );      assert( err == 0 );
    
    // initialize matrices
    size = maxn*maxn;
    lapackf77_dlarnv( &ione, ISEED, &size, A  );
    lapackf77_dlarnv( &ione, ISEED, &size, B  );
    lapackf77_dlarnv( &ione, ISEED, &size, C  );
    
    printf( "========== Level 1 BLAS ==========\n" );
    
    // ----- test DSWAP
    // swap 2nd and 3rd columns of dA, then copy to C2 and compare with A
    printf( "\ntesting dswap\n" );
    assert( k >= 4 );
    magma_dsetmatrix( m, k, A, ld, dA, ld );
    magma_dswap( m, dA(0,1), 1, dA(0,2), 1 );
    magma_dgetmatrix( m, k, dA, ld, C2, ld );
    blasf77_daxpy( &m, &c_neg_one, A(0,0), &ione, C2(0,0), &ione );
    blasf77_daxpy( &m, &c_neg_one, A(0,1), &ione, C2(0,2), &ione );  // swapped
    blasf77_daxpy( &m, &c_neg_one, A(0,2), &ione, C2(0,1), &ione );  // swapped
    blasf77_daxpy( &m, &c_neg_one, A(0,3), &ione, C2(0,3), &ione );
    size = 4;
    error = lapackf77_dlange( "F", &m, &size, C2, &ld, work );
    printf( "dswap diff %.2g\n", error );
    
    // ----- test IDAMAX
    // get argmax of column of A
    printf( "\ntesting idamax\n" );
    magma_dsetmatrix( m, k, A, ld, dA, ld );
    for( int j = 0; j < k; ++j ) {
        int i1 = magma_idamax( m, dA(0,j), 1 );
        int i2 = cublasIdamax( m, dA(0,j), 1 );
        assert( i1 == i2 );
        printf( "i1 %4d, i2 %4d, diff %d\n", i1, i2, i1-i2 );
    }
    
    printf( "\n========== Level 2 BLAS ==========\n" );
    
    // ----- test DGEMV
    // c = alpha*A*b + beta*c,  with A m*n; b,c m or n-vectors
    // try no-trans/trans
    printf( "\ntesting dgemv\n" );
    for( int ia = 0; ia < 3; ++ia ) {
        magma_dsetmatrix( m, n, A,  ld, dA,  ld );
        magma_dsetvector( maxn, B, 1, dB,  1 );
        magma_dsetvector( maxn, C, 1, dC1, 1 );
        magma_dsetvector( maxn, C, 1, dC2, 1 );
        magma_dgemv( trans[ia], m, n, alpha, dA, ld, dB, 1, beta, dC1, 1 );
        cublasDgemv( trans[ia], m, n, alpha, dA, ld, dB, 1, beta, dC2, 1 );
        
        // check results, storing diff between magma and cuda call in C2
        size = (trans[ia] == 'N' ? m : n);
        cublasDaxpy( size, c_neg_one, dC1, 1, dC2, 1 );
        magma_dgetvector( size, dC2, 1, C2, 1 );
        error = lapackf77_dlange( "F", &size, &ione, C2, &ld, work );
        printf( "dgemv( %c ) diff %.2g\n", trans[ia], error );
    }
    
    // ----- test DSYMV
    // c = alpha*A*b + beta*c,  with A m*m symmetric; b,c m-vectors
    // try upper/lower
    printf( "\ntesting dsymv\n" );
    for( int iu = 0; iu < 2; ++iu ) {
        magma_dsetmatrix( m, m, A, ld, dA, ld );
        magma_dsetvector( m, B, 1, dB,  1 );
        magma_dsetvector( m, C, 1, dC1, 1 );
        magma_dsetvector( m, C, 1, dC2, 1 );
        magma_dsymv( uplo[iu], m, alpha, dA, ld, dB, 1, beta, dC1, 1 );
        cublasDsymv( uplo[iu], m, alpha, dA, ld, dB, 1, beta, dC2, 1 );
                                      
        // check results, storing diff between magma and cuda call in C2
        cublasDaxpy( m, c_neg_one, dC1, 1, dC2, 1 );
        magma_dgetvector( m, dC2, 1, C2, 1 );
        error = lapackf77_dlange( "F", &m, &ione, C2, &ld, work );
        printf( "dsymv( %c ) diff %.2g\n", uplo[iu], error );
    }
    
    // ----- test DTRSV
    // solve A*c = c,  with A m*m triangular; c m-vector
    // try upper/lower, no-trans/trans, unit/non-unit diag
    printf( "\ntesting dtrsv\n" );
    // Factor A into LU to get well-conditioned triangles, else solve yields garbage.
    // Still can give garbage if solves aren't consistent with LU factors,
    // e.g., using unit diag for U.
    lapackf77_dgetrf( &m, &m, A, &ld, piv, &info );
    for( int iu = 0; iu < 2; ++iu ) {
    for( int it = 0; it < 3; ++it ) {
    for( int id = 0; id < 2; ++id ) {
        magma_dsetmatrix( m, m, A, ld, dA, ld );
        magma_dsetvector( m, C, 1, dC1, 1 );
        magma_dsetvector( m, C, 1, dC2, 1 );
        magma_dtrsv( uplo[iu], trans[it], diag[id], m, dA, ld, dC1, 1 );
        cublasDtrsv( uplo[iu], trans[it], diag[id], m, dA, ld, dC2, 1 );
                                      
        // check results, storing diff between magma and cuda call in C2
        cublasDaxpy( m, c_neg_one, dC1, 1, dC2, 1 );
        magma_dgetvector( m, dC2, 1, C2, 1 );
        error = lapackf77_dlange( "F", &m, &ione, C2, &ld, work );
        printf( "dtrsv( %c, %c, %c ) diff %.2g\n", uplo[iu], trans[it], diag[id], error );
    }}}
    
    printf( "\n========== Level 3 BLAS ==========\n" );
    
    // ----- test DGEMM
    // C = alpha*A*B + beta*C,  with A m*k or k*m; B k*n or n*k; C m*n
    // try combinations of no-trans/trans
    printf( "\ntesting dgemm\n" );
    for( int ia = 0; ia < 3; ++ia ) {
    for( int ib = 0; ib < 3; ++ib ) {
        bool nta = (trans[ia] == 'N');
        bool ntb = (trans[ib] == 'N');
        magma_dsetmatrix( (nta ? m : k), (nta ? m : k), A, ld, dA,  ld );
        magma_dsetmatrix( (ntb ? k : n), (ntb ? n : k), B, ld, dB,  ld );
        magma_dsetmatrix( m, n, C, ld, dC1, ld );
        magma_dsetmatrix( m, n, C, ld, dC2, ld );
        magma_dgemm( trans[ia], trans[ib], m, n, k, alpha, dA, ld, dB, ld, beta, dC1, ld );
        cublasDgemm( trans[ia], trans[ib], m, n, k, alpha, dA, ld, dB, ld, beta, dC2, ld );
        
        // check results, storing diff between magma and cuda call in C2
        cublasDaxpy( ld*n, c_neg_one, dC1, 1, dC2, 1 );
        magma_dgetmatrix( m, n, dC2, ld, C2, ld );
        error = lapackf77_dlange( "F", &m, &n, C2, &ld, work );
        printf( "dgemm( %c, %c ) diff %.2g\n", trans[ia], trans[ib], error );
    }}
    
    // ----- test DSYMM
    // C = alpha*A*B + beta*C  (left)  with A m*m symmetric; B,C m*n; or
    // C = alpha*B*A + beta*C  (right) with A n*n symmetric; B,C m*n
    // try left/right, upper/lower
    printf( "\ntesting dsymm\n" );
    for( int is = 0; is < 2; ++is ) {
    for( int iu = 0; iu < 2; ++iu ) {
        magma_dsetmatrix( m, m, A, ld, dA,  ld );
        magma_dsetmatrix( m, n, B, ld, dB,  ld );
        magma_dsetmatrix( m, n, C, ld, dC1, ld );
        magma_dsetmatrix( m, n, C, ld, dC2, ld );
        magma_dsymm( side[is], uplo[iu], m, n, alpha, dA, ld, dB, ld, beta, dC1, ld );
        cublasDsymm( side[is], uplo[iu], m, n, alpha, dA, ld, dB, ld, beta, dC2, ld );
        
        // check results, storing diff between magma and cuda call in C2
        cublasDaxpy( ld*n, c_neg_one, dC1, 1, dC2, 1 );
        magma_dgetmatrix( m, n, dC2, ld, C2, ld );
        error = lapackf77_dlange( "F", &m, &n, C2, &ld, work );
        printf( "dsymm( %c, %c ) diff %.2g\n", side[is], uplo[iu], error );
    }}
    
    // ----- test DSYRK
    // C = alpha*A*A^H + beta*C  (no-trans) with A m*k and C m*m symmetric; or
    // C = alpha*A^H*A + beta*C  (trans)    with A k*m and C m*m symmetric
    // try upper/lower, no-trans/trans
    printf( "\ntesting dsyrk\n" );
    for( int iu = 0; iu < 2; ++iu ) {
    for( int it = 0; it < 3; ++it ) {
        magma_dsetmatrix( n, k, A, ld, dA,  ld );
        magma_dsetmatrix( n, n, C, ld, dC1, ld );
        magma_dsetmatrix( n, n, C, ld, dC2, ld );
        magma_dsyrk( uplo[iu], trans[it], n, k, dalpha, dA, ld, dbeta, dC1, ld );
        cublasDsyrk( uplo[iu], trans[it], n, k, dalpha, dA, ld, dbeta, dC2, ld );
        
        // check results, storing diff between magma and cuda call in C2
        cublasDaxpy( ld*n, c_neg_one, dC1, 1, dC2, 1 );
        magma_dgetmatrix( n, n, dC2, ld, C2, ld );
        error = lapackf77_dlange( "F", &n, &n, C2, &ld, work );
        printf( "dsyrk( %c, %c ) diff %.2g\n", uplo[iu], trans[it], error );
    }}
    
    // ----- test DSYR2K
    // C = alpha*A*B^H + ^alpha*B*A^H + beta*C  (no-trans) with A,B n*k; C n*n symmetric; or
    // C = alpha*A^H*B + ^alpha*B^H*A + beta*C  (trans)    with A,B k*n; C n*n symmetric
    // try upper/lower, no-trans/trans
    printf( "\ntesting dsyr2k\n" );
    for( int iu = 0; iu < 2; ++iu ) {
    for( int it = 0; it < 3; ++it ) {
        bool nt = (trans[it] == 'N');
        magma_dsetmatrix( (nt ? n : k), (nt ? n : k), A, ld, dA,  ld );
        magma_dsetmatrix( n, n, C, ld, dC1, ld );
        magma_dsetmatrix( n, n, C, ld, dC2, ld );
        magma_dsyr2k( uplo[iu], trans[it], n, k, alpha, dA, ld, dB, ld, dbeta, dC1, ld );
        cublasDsyr2k( uplo[iu], trans[it], n, k, alpha, dA, ld, dB, ld, dbeta, dC2, ld );
        
        // check results, storing diff between magma and cuda call in C2
        cublasDaxpy( ld*n, c_neg_one, dC1, 1, dC2, 1 );
        magma_dgetmatrix( n, n, dC2, ld, C2, ld );
        error = lapackf77_dlange( "F", &n, &n, C2, &ld, work );
        printf( "dsyr2k( %c, %c ) diff %.2g\n", uplo[iu], trans[it], error );
    }}
    
    // ----- test DTRMM
    // C = alpha*A*C  (left)  with A m*m triangular; C m*n; or
    // C = alpha*C*A  (right) with A n*n triangular; C m*n
    // try left/right, upper/lower, no-trans/trans, unit/non-unit
    printf( "\ntesting dtrmm\n" );
    for( int is = 0; is < 2; ++is ) {
    for( int iu = 0; iu < 2; ++iu ) {
    for( int it = 0; it < 3; ++it ) {
    for( int id = 0; id < 2; ++id ) {
        bool left = (side[is] == 'L');
        magma_dsetmatrix( (left ? m : n), (left ? m : n), A, ld, dA,  ld );
        magma_dsetmatrix( m, n, C, ld, dC1, ld );
        magma_dsetmatrix( m, n, C, ld, dC2, ld );
        magma_dtrmm( side[is], uplo[iu], trans[it], diag[id], m, n, alpha, dA, ld, dC1, ld );
        cublasDtrmm( side[is], uplo[iu], trans[it], diag[id], m, n, alpha, dA, ld, dC2, ld );
        
        // check results, storing diff between magma and cuda call in C2
        cublasDaxpy( ld*n, c_neg_one, dC1, 1, dC2, 1 );
        magma_dgetmatrix( m, n, dC2, ld, C2, ld );
        error = lapackf77_dlange( "F", &n, &n, C2, &ld, work );
        printf( "dtrmm( %c, %c ) diff %.2g\n", uplo[iu], trans[it], error );
    }}}}
    
    // ----- test DTRSM
    // solve A*X = alpha*B  (left)  with A m*m triangular; B m*n; or
    // solve X*A = alpha*B  (right) with A n*n triangular; B m*n
    // try left/right, upper/lower, no-trans/trans, unit/non-unit
    printf( "\ntesting dtrsm\n" );
    for( int is = 0; is < 2; ++is ) {
    for( int iu = 0; iu < 2; ++iu ) {
    for( int it = 0; it < 3; ++it ) {
    for( int id = 0; id < 2; ++id ) {
        bool left = (side[is] == 'L');
        magma_dsetmatrix( (left ? m : n), (left ? m : n), A, ld, dA,  ld );
        magma_dsetmatrix( m, n, C, ld, dC1, ld );
        magma_dsetmatrix( m, n, C, ld, dC2, ld );
        magma_dtrsm( side[is], uplo[iu], trans[it], diag[id], m, n, alpha, dA, ld, dC1, ld );
        cublasDtrsm( side[is], uplo[iu], trans[it], diag[id], m, n, alpha, dA, ld, dC2, ld );
        
        // check results, storing diff between magma and cuda call in C2
        cublasDaxpy( ld*n, c_neg_one, dC1, 1, dC2, 1 );
        magma_dgetmatrix( m, n, dC2, ld, C2, ld );
        error = lapackf77_dlange( "F", &n, &n, C2, &ld, work );
        printf( "dtrsm( %c, %c ) diff %.2g\n", uplo[iu], trans[it], error );
    }}}}
    
    // cleanup
    magma_free_pinned( A  );
    magma_free_pinned( B  );
    magma_free_pinned( C  );
    magma_free_pinned( C2 );
    magma_free( dA  );
    magma_free( dB  );
    magma_free( dC1 );
    magma_free( dC2 );
    
    TESTING_CUDA_FINALIZE();
    return 0;
}
