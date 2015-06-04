/*
 *  -- MAGMA (version 1.3.0) --
 *     Univ. of Tennessee, Knoxville
 *     Univ. of California, Berkeley
 *     Univ. of Colorado, Denver
 *     November 2012
 *
 * @generated d Wed Nov 14 22:54:17 2012
 *
 **/
// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <cuda_runtime_api.h>
#include <cublas.h>

// includes, project
#include "flops.h"
#include "magma.h"
#include "magma_lapack.h"
#include "testings.h"


// Initialize matrix to random.
// Having this in separate function ensures the same ISEED is always used,
// so we can re-generate the identical matrix.
void init_matrix( int m, int n, double *h_A, magma_int_t lda )
{
    magma_int_t ione = 1;
    magma_int_t ISEED[4] = {0,0,0,1};
    magma_int_t n2 = lda*n;
    lapackf77_dlarnv( &ione, ISEED, &n2, h_A );
}


// On input, A and ipiv is LU factorization of A. On output, A is overwritten.
// Requires m == n.
// Uses init_matrix() to re-generate original A as needed.
// Generates random RHS b and solves Ax=b.
// Returns residual, |Ax - b| / (n |A| |x|).
double get_residual(
    magma_int_t m, magma_int_t n,
    double *A, magma_int_t lda,
    magma_int_t *ipiv )
{
    if ( m != n ) {
        printf( "\nERROR: residual check defined only for square matrices\n" );
        return -1;
    }
    
    const double c_one     = MAGMA_D_ONE;
    const double c_neg_one = MAGMA_D_NEG_ONE;
    const magma_int_t ione = 1;
    
    // this seed should be DIFFERENT than used in init_matrix
    // (else x is column of A, so residual can be exactly zero)
    magma_int_t ISEED[4] = {0,0,0,2};
    magma_int_t info = 0;
    double *x, *b;
    
    // initialize RHS
    TESTING_MALLOC( x, double, n );
    TESTING_MALLOC( b, double, n );
    lapackf77_dlarnv( &ione, ISEED, &n, b );
    blasf77_dcopy( &n, b, &ione, x, &ione );
        
    // solve Ax = b
    lapackf77_dgetrs( "Notrans", &n, &ione, A, &lda, ipiv, x, &n, &info );
    if ( info != 0 )
        printf( "lapackf77_dgetrs returned error %d\n", info );
    
    // reset to original A
    init_matrix( m, n, A, lda );
    
    // compute r = Ax - b, saved in b
    blasf77_dgemv( "Notrans", &m, &n, &c_one, A, &lda, x, &ione, &c_neg_one, b, &ione );
    
    // compute residual |Ax - b| / (n*|A|*|x|)
    double norm_x, norm_A, norm_r, work[1];
    norm_A = lapackf77_dlange( "F", &m, &n, A, &lda, work );
    norm_r = lapackf77_dlange( "F", &n, &ione, b, &n, work );
    norm_x = lapackf77_dlange( "F", &n, &ione, x, &n, work );
    
    //printf( "r=\n" ); magma_dprint( 1, n, b, 1 );
    
    TESTING_FREE( x );
    TESTING_FREE( b );
    
    //printf( "r=%.2e, A=%.2e, x=%.2e, n=%d\n", norm_r, norm_A, norm_x, n );
    return norm_r / (n * norm_A * norm_x);
}


// On input, LU and ipiv is LU factorization of A. On output, LU is overwritten.
// Works for any m, n.
// Uses init_matrix() to re-generate original A as needed.
// Returns error in factorization, |PA - LU| / (n |A|)
// This allocates 3 more matrices to store A, L, and U.
double get_LU_error(magma_int_t M, magma_int_t N, 
                    double *LU, magma_int_t lda,
                    magma_int_t *ipiv)
{
    magma_int_t min_mn = min(M,N);
    magma_int_t ione   = 1;
    magma_int_t i, j;
    double alpha = MAGMA_D_ONE;
    double beta  = MAGMA_D_ZERO;
    double *A, *L, *U;
    double work[1], matnorm, residual;
    
    TESTING_MALLOC( A, double, lda*N    );
    TESTING_MALLOC( L, double, M*min_mn );
    TESTING_MALLOC( U, double, min_mn*N );
    memset( L, 0, M*min_mn*sizeof(double) );
    memset( U, 0, min_mn*N*sizeof(double) );

    // set to original A
    init_matrix( M, N, A, lda );
    lapackf77_dlaswp( &N, A, &lda, &ione, &min_mn, ipiv, &ione);
    
    // copy LU to L and U, and set diagonal to 1
    lapackf77_dlacpy( MagmaLowerStr, &M, &min_mn, LU, &lda, L, &M      );
    lapackf77_dlacpy( MagmaUpperStr, &min_mn, &N, LU, &lda, U, &min_mn );
    for(j=0; j<min_mn; j++)
        L[j+j*M] = MAGMA_D_MAKE( 1., 0. );
    
    matnorm = lapackf77_dlange("f", &M, &N, A, &lda, work);

    blasf77_dgemm("N", "N", &M, &N, &min_mn,
                  &alpha, L, &M, U, &min_mn, &beta, LU, &lda);

    for( j = 0; j < N; j++ ) {
        for( i = 0; i < M; i++ ) {
            LU[i+j*lda] = MAGMA_D_SUB( LU[i+j*lda], A[i+j*lda] );
        }
    }
    residual = lapackf77_dlange("f", &M, &N, LU, &lda, work);

    TESTING_FREE(A);
    TESTING_FREE(L);
    TESTING_FREE(U);

    return residual / (matnorm * N);
}


/* ////////////////////////////////////////////////////////////////////////////
   -- Testing dgetrf
*/
int main( int argc, char** argv)
{
    TESTING_CUDA_INIT();

    real_Double_t   gflops, gpu_perf, gpu_time, cpu_perf, cpu_time;
    double          error;
    double *h_A;
    magma_int_t     *ipiv;

    /* Matrix size */
    magma_int_t M = 0, N = 0, n2, lda, ldda;
    const int MAXTESTS = 10;
    magma_int_t msize[MAXTESTS] = { 1024, 2048, 3072, 4032, 5184, 6016, 7040, 8064, 9088, 10112 };
    magma_int_t nsize[MAXTESTS] = { 1024, 2048, 3072, 4032, 5184, 6016, 7040, 8064, 9088, 10112 };

    magma_int_t i, info, min_mn, nb;

    int lapack   = getenv("MAGMA_RUN_LAPCAK")     != NULL;
    int checkres = getenv("MAGMA_TESTINGS_CHECK") != NULL;
    
    // process command line arguments
    printf( "\nUsage: %s -N <m,n> -c -c2 -l\n"
            "  -N  can be repeated up to %d times. If only m is given, then m=n.\n"
            "  -c  or setting $MAGMA_TESTINGS_CHECK checks result, PA - LU.\n"
            "  -c2 for square matrices, solves one RHS and checks residual, Ax - b.\n"
            "  -l  or setting $MAGMA_RUN_LAPACK runs LAPACK.\n\n",
            argv[0], MAXTESTS );
    
    int ntest = 0;
    for( int i = 1; i < argc; ++i ) {
        if ( strcmp("-N", argv[i]) == 0 && i+1 < argc ) {
            magma_assert( ntest < MAXTESTS, "error: -N repeated more than maximum %d tests\n", MAXTESTS );
            int m, n;
            info = sscanf( argv[++i], "%d,%d", &m, &n );
            if ( info == 2 && m > 0 && n > 0 ) {
                msize[ ntest ] = m;
                nsize[ ntest ] = n;
            }
            else if ( info == 1 && m > 0 ) {
                msize[ ntest ] = m;
                nsize[ ntest ] = m;  // implicitly
            }
            else {
                printf( "error: -N %s is invalid; ensure m > 0, n > 0.\n", argv[i] );
                exit(1);
            }
            M = max( M, msize[ ntest ] );
            N = max( N, nsize[ ntest ] );
            ntest++;
        }
        else if ( strcmp("-M", argv[i]) == 0 ) {
            printf( "-M has been replaced in favor of -N m,n to allow -N to be repeated.\n\n" );
            exit(1);
        }
        else if ( strcmp("-c", argv[i]) == 0 ) {
            checkres = 1;
        }
        else if ( strcmp("-c2", argv[i]) == 0 ) {
            checkres = 2;
        }
        else if ( strcmp("-l", argv[i]) == 0 ) {
            lapack = true;
        }
        else {
            printf( "invalid argument: %s\n", argv[i] );
            exit(1);
        }
    }
    if ( ntest == 0 ) {
        ntest = MAXTESTS;
        M = msize[ntest-1];
        N = nsize[ntest-1];
    }
    
    ldda   = ((M+31)/32)*32;
    n2     = M * N;
    min_mn = min(M, N);
    nb     = magma_get_dgetrf_nb(min_mn);

    /* Allocate memory for the matrix */
    TESTING_MALLOC(ipiv, magma_int_t, min_mn);
    TESTING_HOSTALLOC( h_A, double, n2 );

    if ( checkres == 2 ) {
        printf("    M     N   CPU GFlop/s (sec)   GPU GFlop/s (sec)   |Ax-b|/(N*|A|*|x|)\n");
    }
    else {
        printf("    M     N   CPU GFlop/s (sec)   GPU GFlop/s (sec)   |PA-LU|/(N*|A|)\n");
    }
    printf("=========================================================================\n");
    for( i = 0; i < ntest; ++i ) {
        M = msize[i];
        N = nsize[i];
        min_mn = min(M, N);
        lda    = M;
        n2     = lda*N;
        ldda   = ((M+31)/32)*32;
        gflops = FLOPS_DGETRF( M, N ) / 1e9;

        /* =====================================================================
           Performs operation using LAPACK
           =================================================================== */
        if ( lapack ) {
            init_matrix( M, N, h_A, lda );
            
            cpu_time = magma_wtime();
            lapackf77_dgetrf(&M, &N, h_A, &lda, ipiv, &info);
            cpu_time = magma_wtime() - cpu_time;
            cpu_perf = gflops / cpu_time;
            if (info != 0)
                printf("lapackf77_dgetrf returned error %d.\n", (int) info);
        }

        /* ====================================================================
           Performs operation using MAGMA
           =================================================================== */
        init_matrix( M, N, h_A, lda );
        
        gpu_time = magma_wtime();
        magma_dgetrf( M, N, h_A, lda, ipiv, &info);
        gpu_time = magma_wtime() - gpu_time;
        gpu_perf = gflops / gpu_time;
        if (info != 0)
            printf("magma_dgetrf returned error %d.\n", (int) info);

        /* =====================================================================
           Check the factorization
           =================================================================== */
        printf("%5d %5d", (int) M, (int) N );
        if ( lapack ) {
            printf("   %7.2f (%7.2f)", cpu_perf, cpu_time );
        }
        else {
            printf("     ---   (  ---  )");
        }
        printf("   %7.2f (%7.2f)", gpu_perf, gpu_time );
        if ( checkres == 2 ) {
            error = get_residual( M, N, h_A, lda, ipiv );
            printf("   %8.2e\n", error );            
        }
        else if ( checkres ) {
            error = get_LU_error(M, N, h_A, lda, ipiv );
            printf("   %8.2e\n", error );
        }
        else {
            printf("     ---   \n");
        }
    }

    /* Memory clean up */
    TESTING_FREE( ipiv );
    TESTING_HOSTFREE( h_A );

    TESTING_CUDA_FINALIZE();
    return 0;
}
