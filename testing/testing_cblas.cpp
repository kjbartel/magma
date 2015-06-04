/*
 *  -- MAGMA (version 1.2.1) --
 *     Univ. of Tennessee, Knoxville
 *     Univ. of California, Berkeley
 *     Univ. of Colorado, Denver
 *     June 2012
 *
 * @author Mark Gates
 * @generated c Thu Jun 28 12:31:51 2012
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
    
    cuFloatComplex c_neg_one = MAGMA_C_NEG_ONE;
    magma_int_t ione = 1;
    const char trans[] = { 'N', 'C', 'T' };
    const char uplo[]  = { 'L', 'U' };
    const char diag[]  = { 'U', 'N' };
    const char side[]  = { 'L', 'R' };
    
    cuFloatComplex  *A,  *B,  *C,   *C2;
    cuFloatComplex *dA, *dB, *dC1, *dC2;
    cuFloatComplex alpha = MAGMA_C_MAKE( 0.5, 0.1 );
    cuFloatComplex beta  = MAGMA_C_MAKE( 0.7, 0.2 );
    float dalpha = 0.6;
    float dbeta  = 0.8;
    float work[1], error;
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
    err = magma_cmalloc_pinned( &A , size );  assert( err == 0 );
    err = magma_cmalloc_pinned( &B , size );  assert( err == 0 );
    err = magma_cmalloc_pinned( &C , size );  assert( err == 0 );
    err = magma_cmalloc_pinned( &C2, size );  assert( err == 0 );
    err = magma_cmalloc( &dA,  size );      assert( err == 0 );
    err = magma_cmalloc( &dB,  size );      assert( err == 0 );
    err = magma_cmalloc( &dC1, size );      assert( err == 0 );
    err = magma_cmalloc( &dC2, size );      assert( err == 0 );
    
    // initialize matrices
    size = maxn*maxn;
    lapackf77_clarnv( &ione, ISEED, &size, A  );
    lapackf77_clarnv( &ione, ISEED, &size, B  );
    lapackf77_clarnv( &ione, ISEED, &size, C  );
    
    printf( "========== Level 1 BLAS ==========\n" );
    
    // ----- test CSWAP
    // swap 2nd and 3rd columns of dA, then copy to C2 and compare with A
    printf( "\ntesting cswap\n" );
    assert( k >= 4 );
    magma_csetmatrix( m, k, A, ld, dA, ld );
    magma_cswap( m, dA(0,1), 1, dA(0,2), 1 );
    magma_cgetmatrix( m, k, dA, ld, C2, ld );
    blasf77_caxpy( &m, &c_neg_one, A(0,0), &ione, C2(0,0), &ione );
    blasf77_caxpy( &m, &c_neg_one, A(0,1), &ione, C2(0,2), &ione );  // swapped
    blasf77_caxpy( &m, &c_neg_one, A(0,2), &ione, C2(0,1), &ione );  // swapped
    blasf77_caxpy( &m, &c_neg_one, A(0,3), &ione, C2(0,3), &ione );
    size = 4;
    error = lapackf77_clange( "F", &m, &size, C2, &ld, work );
    printf( "cswap diff %.2g\n", error );
    
    // ----- test ICAMAX
    // get argmax of column of A
    printf( "\ntesting icamax\n" );
    magma_csetmatrix( m, k, A, ld, dA, ld );
    for( int j = 0; j < k; ++j ) {
        int i1 = magma_icamax( m, dA(0,j), 1 );
        int i2 = cublasIcamax( m, dA(0,j), 1 );
        assert( i1 == i2 );
        printf( "i1 %4d, i2 %4d, diff %d\n", i1, i2, i1-i2 );
    }
    
    printf( "\n========== Level 2 BLAS ==========\n" );
    
    // ----- test CGEMV
    // c = alpha*A*b + beta*c,  with A m*n; b,c m or n-vectors
    // try no-trans/trans
    printf( "\ntesting cgemv\n" );
    for( int ia = 0; ia < 3; ++ia ) {
        magma_csetmatrix( m, n, A,  ld, dA,  ld );
        magma_csetvector( maxn, B, 1, dB,  1 );
        magma_csetvector( maxn, C, 1, dC1, 1 );
        magma_csetvector( maxn, C, 1, dC2, 1 );
        magma_cgemv( trans[ia], m, n, alpha, dA, ld, dB, 1, beta, dC1, 1 );
        cublasCgemv( trans[ia], m, n, alpha, dA, ld, dB, 1, beta, dC2, 1 );
        
        // check results, storing diff between magma and cuda call in C2
        size = (trans[ia] == 'N' ? m : n);
        cublasCaxpy( size, c_neg_one, dC1, 1, dC2, 1 );
        magma_cgetvector( size, dC2, 1, C2, 1 );
        error = lapackf77_clange( "F", &size, &ione, C2, &ld, work );
        printf( "cgemv( %c ) diff %.2g\n", trans[ia], error );
    }
    
    // ----- test CHEMV
    // c = alpha*A*b + beta*c,  with A m*m symmetric; b,c m-vectors
    // try upper/lower
    printf( "\ntesting chemv\n" );
    for( int iu = 0; iu < 2; ++iu ) {
        magma_csetmatrix( m, m, A, ld, dA, ld );
        magma_csetvector( m, B, 1, dB,  1 );
        magma_csetvector( m, C, 1, dC1, 1 );
        magma_csetvector( m, C, 1, dC2, 1 );
        magma_chemv( uplo[iu], m, alpha, dA, ld, dB, 1, beta, dC1, 1 );
        cublasChemv( uplo[iu], m, alpha, dA, ld, dB, 1, beta, dC2, 1 );
                                      
        // check results, storing diff between magma and cuda call in C2
        cublasCaxpy( m, c_neg_one, dC1, 1, dC2, 1 );
        magma_cgetvector( m, dC2, 1, C2, 1 );
        error = lapackf77_clange( "F", &m, &ione, C2, &ld, work );
        printf( "chemv( %c ) diff %.2g\n", uplo[iu], error );
    }
    
    // ----- test CTRSV
    // solve A*c = c,  with A m*m triangular; c m-vector
    // try upper/lower, no-trans/trans, unit/non-unit diag
    printf( "\ntesting ctrsv\n" );
    // Factor A into LU to get well-conditioned triangles, else solve yields garbage.
    // Still can give garbage if solves aren't consistent with LU factors,
    // e.g., using unit diag for U.
    lapackf77_cgetrf( &m, &m, A, &ld, piv, &info );
    for( int iu = 0; iu < 2; ++iu ) {
    for( int it = 0; it < 3; ++it ) {
    for( int id = 0; id < 2; ++id ) {
        magma_csetmatrix( m, m, A, ld, dA, ld );
        magma_csetvector( m, C, 1, dC1, 1 );
        magma_csetvector( m, C, 1, dC2, 1 );
        magma_ctrsv( uplo[iu], trans[it], diag[id], m, dA, ld, dC1, 1 );
        cublasCtrsv( uplo[iu], trans[it], diag[id], m, dA, ld, dC2, 1 );
                                      
        // check results, storing diff between magma and cuda call in C2
        cublasCaxpy( m, c_neg_one, dC1, 1, dC2, 1 );
        magma_cgetvector( m, dC2, 1, C2, 1 );
        error = lapackf77_clange( "F", &m, &ione, C2, &ld, work );
        printf( "ctrsv( %c, %c, %c ) diff %.2g\n", uplo[iu], trans[it], diag[id], error );
    }}}
    
    printf( "\n========== Level 3 BLAS ==========\n" );
    
    // ----- test CGEMM
    // C = alpha*A*B + beta*C,  with A m*k or k*m; B k*n or n*k; C m*n
    // try combinations of no-trans/trans
    printf( "\ntesting cgemm\n" );
    for( int ia = 0; ia < 3; ++ia ) {
    for( int ib = 0; ib < 3; ++ib ) {
        bool nta = (trans[ia] == 'N');
        bool ntb = (trans[ib] == 'N');
        magma_csetmatrix( (nta ? m : k), (nta ? m : k), A, ld, dA,  ld );
        magma_csetmatrix( (ntb ? k : n), (ntb ? n : k), B, ld, dB,  ld );
        magma_csetmatrix( m, n, C, ld, dC1, ld );
        magma_csetmatrix( m, n, C, ld, dC2, ld );
        magma_cgemm( trans[ia], trans[ib], m, n, k, alpha, dA, ld, dB, ld, beta, dC1, ld );
        cublasCgemm( trans[ia], trans[ib], m, n, k, alpha, dA, ld, dB, ld, beta, dC2, ld );
        
        // check results, storing diff between magma and cuda call in C2
        cublasCaxpy( ld*n, c_neg_one, dC1, 1, dC2, 1 );
        magma_cgetmatrix( m, n, dC2, ld, C2, ld );
        error = lapackf77_clange( "F", &m, &n, C2, &ld, work );
        printf( "cgemm( %c, %c ) diff %.2g\n", trans[ia], trans[ib], error );
    }}
    
    // ----- test CHEMM
    // C = alpha*A*B + beta*C  (left)  with A m*m symmetric; B,C m*n; or
    // C = alpha*B*A + beta*C  (right) with A n*n symmetric; B,C m*n
    // try left/right, upper/lower
    printf( "\ntesting chemm\n" );
    for( int is = 0; is < 2; ++is ) {
    for( int iu = 0; iu < 2; ++iu ) {
        magma_csetmatrix( m, m, A, ld, dA,  ld );
        magma_csetmatrix( m, n, B, ld, dB,  ld );
        magma_csetmatrix( m, n, C, ld, dC1, ld );
        magma_csetmatrix( m, n, C, ld, dC2, ld );
        magma_chemm( side[is], uplo[iu], m, n, alpha, dA, ld, dB, ld, beta, dC1, ld );
        cublasChemm( side[is], uplo[iu], m, n, alpha, dA, ld, dB, ld, beta, dC2, ld );
        
        // check results, storing diff between magma and cuda call in C2
        cublasCaxpy( ld*n, c_neg_one, dC1, 1, dC2, 1 );
        magma_cgetmatrix( m, n, dC2, ld, C2, ld );
        error = lapackf77_clange( "F", &m, &n, C2, &ld, work );
        printf( "chemm( %c, %c ) diff %.2g\n", side[is], uplo[iu], error );
    }}
    
    // ----- test CHERK
    // C = alpha*A*A^H + beta*C  (no-trans) with A m*k and C m*m symmetric; or
    // C = alpha*A^H*A + beta*C  (trans)    with A k*m and C m*m symmetric
    // try upper/lower, no-trans/trans
    printf( "\ntesting cherk\n" );
    for( int iu = 0; iu < 2; ++iu ) {
    for( int it = 0; it < 3; ++it ) {
        magma_csetmatrix( n, k, A, ld, dA,  ld );
        magma_csetmatrix( n, n, C, ld, dC1, ld );
        magma_csetmatrix( n, n, C, ld, dC2, ld );
        magma_cherk( uplo[iu], trans[it], n, k, dalpha, dA, ld, dbeta, dC1, ld );
        cublasCherk( uplo[iu], trans[it], n, k, dalpha, dA, ld, dbeta, dC2, ld );
        
        // check results, storing diff between magma and cuda call in C2
        cublasCaxpy( ld*n, c_neg_one, dC1, 1, dC2, 1 );
        magma_cgetmatrix( n, n, dC2, ld, C2, ld );
        error = lapackf77_clange( "F", &n, &n, C2, &ld, work );
        printf( "cherk( %c, %c ) diff %.2g\n", uplo[iu], trans[it], error );
    }}
    
    // ----- test CHER2K
    // C = alpha*A*B^H + ^alpha*B*A^H + beta*C  (no-trans) with A,B n*k; C n*n symmetric; or
    // C = alpha*A^H*B + ^alpha*B^H*A + beta*C  (trans)    with A,B k*n; C n*n symmetric
    // try upper/lower, no-trans/trans
    printf( "\ntesting cher2k\n" );
    for( int iu = 0; iu < 2; ++iu ) {
    for( int it = 0; it < 3; ++it ) {
        bool nt = (trans[it] == 'N');
        magma_csetmatrix( (nt ? n : k), (nt ? n : k), A, ld, dA,  ld );
        magma_csetmatrix( n, n, C, ld, dC1, ld );
        magma_csetmatrix( n, n, C, ld, dC2, ld );
        magma_cher2k( uplo[iu], trans[it], n, k, alpha, dA, ld, dB, ld, dbeta, dC1, ld );
        cublasCher2k( uplo[iu], trans[it], n, k, alpha, dA, ld, dB, ld, dbeta, dC2, ld );
        
        // check results, storing diff between magma and cuda call in C2
        cublasCaxpy( ld*n, c_neg_one, dC1, 1, dC2, 1 );
        magma_cgetmatrix( n, n, dC2, ld, C2, ld );
        error = lapackf77_clange( "F", &n, &n, C2, &ld, work );
        printf( "cher2k( %c, %c ) diff %.2g\n", uplo[iu], trans[it], error );
    }}
    
    // ----- test CTRMM
    // C = alpha*A*C  (left)  with A m*m triangular; C m*n; or
    // C = alpha*C*A  (right) with A n*n triangular; C m*n
    // try left/right, upper/lower, no-trans/trans, unit/non-unit
    printf( "\ntesting ctrmm\n" );
    for( int is = 0; is < 2; ++is ) {
    for( int iu = 0; iu < 2; ++iu ) {
    for( int it = 0; it < 3; ++it ) {
    for( int id = 0; id < 2; ++id ) {
        bool left = (side[is] == 'L');
        magma_csetmatrix( (left ? m : n), (left ? m : n), A, ld, dA,  ld );
        magma_csetmatrix( m, n, C, ld, dC1, ld );
        magma_csetmatrix( m, n, C, ld, dC2, ld );
        magma_ctrmm( side[is], uplo[iu], trans[it], diag[id], m, n, alpha, dA, ld, dC1, ld );
        cublasCtrmm( side[is], uplo[iu], trans[it], diag[id], m, n, alpha, dA, ld, dC2, ld );
        
        // check results, storing diff between magma and cuda call in C2
        cublasCaxpy( ld*n, c_neg_one, dC1, 1, dC2, 1 );
        magma_cgetmatrix( m, n, dC2, ld, C2, ld );
        error = lapackf77_clange( "F", &n, &n, C2, &ld, work );
        printf( "ctrmm( %c, %c ) diff %.2g\n", uplo[iu], trans[it], error );
    }}}}
    
    // ----- test CTRSM
    // solve A*X = alpha*B  (left)  with A m*m triangular; B m*n; or
    // solve X*A = alpha*B  (right) with A n*n triangular; B m*n
    // try left/right, upper/lower, no-trans/trans, unit/non-unit
    printf( "\ntesting ctrsm\n" );
    for( int is = 0; is < 2; ++is ) {
    for( int iu = 0; iu < 2; ++iu ) {
    for( int it = 0; it < 3; ++it ) {
    for( int id = 0; id < 2; ++id ) {
        bool left = (side[is] == 'L');
        magma_csetmatrix( (left ? m : n), (left ? m : n), A, ld, dA,  ld );
        magma_csetmatrix( m, n, C, ld, dC1, ld );
        magma_csetmatrix( m, n, C, ld, dC2, ld );
        magma_ctrsm( side[is], uplo[iu], trans[it], diag[id], m, n, alpha, dA, ld, dC1, ld );
        cublasCtrsm( side[is], uplo[iu], trans[it], diag[id], m, n, alpha, dA, ld, dC2, ld );
        
        // check results, storing diff between magma and cuda call in C2
        cublasCaxpy( ld*n, c_neg_one, dC1, 1, dC2, 1 );
        magma_cgetmatrix( m, n, dC2, ld, C2, ld );
        error = lapackf77_clange( "F", &n, &n, C2, &ld, work );
        printf( "ctrsm( %c, %c ) diff %.2g\n", uplo[iu], trans[it], error );
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
