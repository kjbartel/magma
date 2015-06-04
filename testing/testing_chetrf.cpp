/*
    -- MAGMA (version 1.6.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date November 2014

       @generated from testing_zhetrf.cpp normal z -> c, Sat Nov 15 19:54:18 2014
       @author Ichitaro Yamazaki
*/
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

/* ================================================================================================== */

// Initialize matrix to random.
// Having this in separate function ensures the same ISEED is always used,
// so we can re-generate the identical matrix.
void init_matrix( int nopiv, int m, int n, magmaFloatComplex *h_A, magma_int_t lda )
{
    magma_int_t ione = 1;
    magma_int_t ISEED[4] = {0,0,0,1};
    magma_int_t n2 = lda*n;
    //float *A = (float*)malloc(n2*sizeof(float));
    //lapackf77_slarnv( &ione, ISEED, &n2, A );
    //for (int i=0; i<n; i++) for (int j=0; j<=i; j++) h_A[j+i*lda] = MAGMA_C_MAKE(A[j+i*lda],0.0);
    //free(A);
    //
    lapackf77_clarnv( &ione, ISEED, &n2, h_A );
    // symmetrize
    for (int i=0; i<n; i++) for (int j=0; j<i; j++) h_A[i+j*lda] = MAGMA_C_CNJG(h_A[j+i*lda]);
    if (nopiv) for (int i=0; i<n; i++) h_A[i+i*lda] = MAGMA_C_MAKE(MAGMA_C_REAL(h_A[i+i*lda]) + n, 0.0);
    else       for (int i=0; i<n; i++) h_A[i+i*lda] = MAGMA_C_MAKE(MAGMA_C_REAL(h_A[i+i*lda]), 0.0);
}


// On input, A and ipiv is LU factorization of A. On output, A is overwritten.
// Requires m == n.
// Uses init_matrix() to re-generate original A as needed.
// Generates random RHS b and solves Ax=b.
// Returns residual, |Ax - b| / (n |A| |x|).
float get_residual(
    int nopiv, magma_uplo_t uplo, magma_int_t n,
    magmaFloatComplex *A, magma_int_t lda,
    magma_int_t *ipiv )
{
    const magmaFloatComplex c_one     = MAGMA_C_ONE;
    const magmaFloatComplex c_neg_one = MAGMA_C_NEG_ONE;
    const magma_int_t ione = 1;
    magma_int_t upper = (uplo == MagmaUpper);
    
    // this seed should be DIFFERENT than used in init_matrix
    // (else x is column of A, so residual can be exactly zero)
    magma_int_t ISEED[4] = {0,0,0,2};
    magma_int_t info = 0;
    magmaFloatComplex *x, *b;
    
    // initialize RHS
    TESTING_MALLOC_CPU( x, magmaFloatComplex, n );
    TESTING_MALLOC_CPU( b, magmaFloatComplex, n );
    lapackf77_clarnv( &ione, ISEED, &n, b );
    blasf77_ccopy( &n, b, &ione, x, &ione );
    
    // solve Ax = b
    if (nopiv) {
        if (upper) {
            blasf77_ctrsm( MagmaLeftStr, MagmaUpperStr, 
                           MagmaConjTransStr, MagmaUnitStr, 
                           &n, &ione, &c_one,
                           A, &lda, x, &n );
            for (int i=0; i<n; i++) x[i] = MAGMA_C_DIV( x[i], A[i+i*lda] );
            blasf77_ctrsm( MagmaLeftStr, MagmaUpperStr, 
                           MagmaNoTransStr, MagmaUnitStr, 
                           &n, &ione, &c_one,
                           A, &lda, x, &n );
        } else {
            blasf77_ctrsm( MagmaLeftStr, MagmaLowerStr, 
                           MagmaNoTransStr, MagmaUnitStr, 
                           &n, &ione, &c_one,
                           A, &lda, x, &n );
            for (int i=0; i<n; i++) x[i] = MAGMA_C_DIV( x[i], A[i+i*lda] );
            blasf77_ctrsm( MagmaLeftStr, MagmaLowerStr, 
                           MagmaConjTransStr, MagmaUnitStr, 
                           &n, &ione, &c_one,
                           A, &lda, x, &n );
        }
    }else {
        lapackf77_chetrs( (upper ? MagmaUpperStr: MagmaLowerStr), &n, &ione, A, &lda, ipiv, x, &n, &info );
    }
    if (info != 0)
        printf("lapackf77_chetrs returned error %d: %s.\n",
               (int) info, magma_strerror( info ));
    // reset to original A
    init_matrix( nopiv, n, n, A, lda );
    
    // compute r = Ax - b, saved in b
    blasf77_cgemv( "Notrans", &n, &n, &c_one, A, &lda, x, &ione, &c_neg_one, b, &ione );
    
    // compute residual |Ax - b| / (n*|A|*|x|)
    float norm_x, norm_A, norm_r, work[1];
    norm_A = lapackf77_clange( MagmaFullStr, &n, &n, A, &lda, work );
    norm_r = lapackf77_clange( MagmaFullStr, &n, &ione, b, &n, work );
    norm_x = lapackf77_clange( MagmaFullStr, &n, &ione, x, &n, work );
    
    //printf( "r=\n" ); magma_cprint( 1, n, b, 1 );
    
    TESTING_FREE_CPU( x );
    TESTING_FREE_CPU( b );
    
    //printf( "r=%.2e, A=%.2e, x=%.2e, n=%d\n", norm_r, norm_A, norm_x, n );
    return norm_r / (n * norm_A * norm_x);
}


// On input, LU and ipiv is LU factorization of A. On output, LU is overwritten.
// Works for any m, n.
// Uses init_matrix() to re-generate original A as needed.
// Returns error in factorization, |PA - LU| / (n |A|)
// This allocates 3 more matrices to store A, L, and U.
float get_LDLt_error(int nopiv, magma_uplo_t uplo, magma_int_t N,
                      magmaFloatComplex *LD, magma_int_t lda,
                      magma_int_t *ipiv)
{
    magma_int_t i, j;
    magmaFloatComplex c_one  = MAGMA_C_ONE;
    magmaFloatComplex c_zero = MAGMA_C_ZERO;
    magmaFloatComplex *A, *L, *D;
    float work[1], matnorm, residual;
    #define LD(i,j) (LD[(i) + (j)*lda])
    #define  A(i,j) ( A[(i) + (j)*N])
    #define  L(i,j) ( L[(i) + (j)*N])
    #define  D(i,j) ( D[(i) + (j)*N])

    TESTING_MALLOC_CPU( A, magmaFloatComplex, N*N );
    TESTING_MALLOC_CPU( L, magmaFloatComplex, N*N );
    TESTING_MALLOC_CPU( D, magmaFloatComplex, N*N );
    memset( L, 0, N*N*sizeof(magmaFloatComplex) );
    memset( D, 0, N*N*sizeof(magmaFloatComplex) );

    // set to original A, and apply pivoting
    init_matrix( nopiv, N, N, A, N );
    if (uplo == MagmaUpper) {
        for (j=N-1; j>=0; j--) {
            int piv = (nopiv ? j+1 : ipiv[j]);
            if (piv < 0) {
                piv = -(piv+1);
                // extract 2-by-2 pivot
                D(j,j)     = LD(j,j);
                D(j,j-1)   = MAGMA_C_CNJG(LD(j-1,j));
                D(j-1,j)   = LD(j-1,j);
                D(j-1,j-1) = LD(j-1,j-1);
                // exract L
                L(j,j) = c_one;
                for (i=0; i<j-1; i++) L(i,j) = LD(i,j);
                j--;
                L(j,j) = c_one;
                for (i=0; i<j; i++) L(i,j) = LD(i,j);
                if (piv != j) {
                    // apply row-pivoting to previous L
                    for (i=j+2; i<N; i++) {
                        magmaFloatComplex val = L(j,i);
                        L(j,i) = L(piv,i);
                        L(piv,i) = val;
                    }
                    // apply row-pivoting to A
                    for (i=0; i<N; i++) {
                        magmaFloatComplex val = A(j,i);
                        A(j,i) = A(piv,i);
                        A(piv,i) = val;
                    }
                    // apply col-pivoting to A
                    for (i=0; i<N; i++) {
                        magmaFloatComplex val = A(i,j);
                        A(i,j) = A(i,piv);
                        A(i,piv) = val;
                    }
                }
            } else {
                piv = piv-1;
                // extract 1-by-1 pivot
                D(j,j) = LD(j,j);
                // exract L
                L(j,j) = c_one;
                for (i=0; i<j; i++) L(i,j) = LD(i,j);
                if (piv != j) {
                    // apply row-pivoting to previous L
                    for (i=j+1; i<N; i++) {
                        magmaFloatComplex val = L(j,i);
                        L(j,i) = L(piv,i);
                        L(piv,i) = val;
                    }
                    // apply row-pivoting to A
                    for (i=0; i<N; i++) {
                        magmaFloatComplex val = A(j,i);
                        A(j,i) = A(piv,i);
                        A(piv,i) = val;
                    }
                    // apply col-pivoting to A
                    for (i=0; i<N; i++) {
                        magmaFloatComplex val = A(i,j);
                        A(i,j) = A(i,piv);
                        A(i,piv) = val;
                    }
                }
            }
        }
        if (nopiv) {
            // compute W = D*U
            blasf77_cgemm(MagmaNoTransStr, MagmaNoTransStr, &N, &N, &N,
                          &c_one, D, &N, L, &N, &c_zero, LD, &lda);
            // compute D = U'*W
            blasf77_cgemm(MagmaConjTransStr, MagmaNoTransStr, &N, &N, &N,
                          &c_one, L, &N, LD, &lda, &c_zero, D, &N);
        } else {
            // compute W = U*D
            blasf77_cgemm(MagmaNoTransStr, MagmaNoTransStr, &N, &N, &N,
                          &c_one, L, &N, D, &N, &c_zero, LD, &lda);
            // compute D = W*U'
            blasf77_cgemm(MagmaNoTransStr, MagmaConjTransStr, &N, &N, &N,
                          &c_one, LD, &lda, L, &N, &c_zero, D, &N);
        }
    } else {
        for (j=0; j<N; j++) {
            int piv = (nopiv ? j+1 : ipiv[j]);
            if (piv < 0) {
                piv = -(piv+1);
                // extract 2-by-2 pivot
                D(j,j)     = LD(j,j);
                D(j,j+1)   = MAGMA_C_CNJG(LD(j+1,j));
                D(j+1,j)   = LD(j+1,j);
                D(j+1,j+1) = LD(j+1,j+1);
                // exract L
                L(j,j) = c_one;
                for (i=j+2; i<N; i++) L(i,j) = LD(i,j);
                j++;
                L(j,j) = c_one;
                for (i=j+1; i<N; i++) L(i,j) = LD(i,j);
                if (piv != j) {
                    // apply row-pivoting to previous L
                    for (i=0; i<j-1; i++) {
                        magmaFloatComplex val = L(j,i);
                        L(j,i) = L(piv,i);
                        L(piv,i) = val;
                    }
                    // apply row-pivoting to A
                    for (i=0; i<N; i++) {
                        magmaFloatComplex val = A(j,i);
                        A(j,i) = A(piv,i);
                        A(piv,i) = val;
                    }
                    // apply col-pivoting to A
                    for (i=0; i<N; i++) {
                        magmaFloatComplex val = A(i,j);
                        A(i,j) = A(i,piv);
                        A(i,piv) = val;
                    }
                }
            } else {
                piv = piv-1;
                // extract 1-by-1 pivot
                D(j,j) = LD(j,j);
                // exract L
                L(j,j) = c_one;
                for (i=j+1; i<N; i++) L(i,j) = LD(i,j);
                if (piv != j) {
                    // apply row-pivoting to previous L
                    for (i=0; i<j; i++) {
                        magmaFloatComplex val = L(j,i);
                        L(j,i) = L(piv,i);
                        L(piv,i) = val;
                    }
                    // apply row-pivoting to A
                    for (i=0; i<N; i++) {
                        magmaFloatComplex val = A(j,i);
                        A(j,i) = A(piv,i);
                        A(piv,i) = val;
                    }
                    // apply col-pivoting to A
                    for (i=0; i<N; i++) {
                        magmaFloatComplex val = A(i,j);
                        A(i,j) = A(i,piv);
                        A(i,piv) = val;
                    }
                }
            }
        }
        // compute W = L*D
        blasf77_cgemm(MagmaNoTransStr, MagmaNoTransStr, &N, &N, &N,
                      &c_one, L, &N, D, &N, &c_zero, LD, &lda);
        // compute D = W*L'
        blasf77_cgemm(MagmaNoTransStr, MagmaConjTransStr, &N, &N, &N,
                      &c_one, LD, &lda, L, &N, &c_zero, D, &N);
    }
    // compute norm of A
    matnorm = lapackf77_clange(MagmaFullStr, &N, &N, A, &lda, work);

    for( j = 0; j < N; j++ ) {
        for( i = 0; i < N; i++ ) {
            D(i,j) = MAGMA_C_SUB( D(i,j), A(i,j) );
        }
    }
    residual = lapackf77_clange(MagmaFullStr, &N, &N, D, &N, work);

    TESTING_FREE_CPU( A );
    TESTING_FREE_CPU( L );
    TESTING_FREE_CPU( D );

    return residual / (matnorm * N);
}


/* ////////////////////////////////////////////////////////////////////////////
   -- Testing chetrf
*/
int main( int argc, char** argv)
{
    TESTING_INIT();

    magmaFloatComplex *h_A, *work, temp;
    real_Double_t   gflops, gpu_perf, gpu_time = 0.0, cpu_perf=0, cpu_time=0;
    float          error, error_lapack = 0.0;
    magma_int_t     *ipiv;
    magma_int_t     N, n2, lda, lwork, info;
    magma_int_t     status = 0;
    magma_int_t     cpu = 0, gpu = 0, nopiv = 0, nopiv_gpu = 0, row = 0;
    
    magma_opts opts;
    parse_opts( argc, argv, &opts );
    switch (opts.version) {
        case 1:
            cpu = 1;
            printf( "\nCPU-Interface to Bunch-Kauffman on GPU" );
            break;
        case 2:
            gpu = 1;
            printf( "\nGPU-Interface to Bunch-Kauffman on GPU" );
            printf( "\n not yet..\n\n" );
            return 0;
            break;
        case 3:
            nopiv = 1;
            printf( "\nCPU-Interface to hybrid Non-pivoted LDLt (A is SPD)" );
            break;
        case 4:
            nopiv_gpu = 1;
            printf( "\nGPU-Interface to hybrid Non-pivoted LDLt (A is SPD)" );
            break;
            break;
        //case 5:
        //    row = 1;
        //    printf( "\n Bunch-Kauffman: GPU-only version (row-major)" );
        //    break;
        default:
        //  printf( " hybrid CPU-GPU version" );
            printf( " version = %d not supported\n\n",opts.version);
            return 0;
    }

    magma_uplo_t uplo = opts.uplo;
    printf( " (%s)\n",(uplo == MagmaUpper ? "upper" : "lower") );
    printf( " (--version: 1 = Bunch-Kauffman (CPU), 2 = Bunch-Kauffman (GPU), 3 = No-piv (CPU), 4 = No-piv (GPU))\n\n" );
    
    magma_int_t upper = (uplo == MagmaUpper);
    float tol = opts.tolerance * lapackf77_slamch("E");

    if ( opts.check == 2 ) {
        printf("    M     N   CPU GFlop/s (sec)   GPU GFlop/s (sec)   |Ax-b|/(N*|A|*|x|)\n");
    }
    else {
        printf("    M     N   CPU GFlop/s (sec)   GPU GFlop/s (sec)   |PAP'-LDL'|/(N*|A|)\n");
    }
    printf("=========================================================================\n");
    for( int itest = 0; itest < opts.ntest; ++itest ) {
        for( int iter = 0; iter < opts.niter; ++iter ) {
            N = opts.nsize[itest];
            lda    = N;
            n2     = lda*N;
            gflops = FLOPS_CPOTRF( N ) / 1e9;
            
            TESTING_MALLOC_CPU( ipiv, magma_int_t, N );
            TESTING_MALLOC_PIN( h_A,  magmaFloatComplex, n2 );
            

            /* =====================================================================
               Performs operation using LAPACK
               =================================================================== */
            if ( opts.lapack ) {
                lwork = -1;
                lapackf77_chetrf((upper ? MagmaUpperStr: MagmaLowerStr), &N, h_A, &lda, ipiv, &temp, &lwork, &info);
                lwork = (int)MAGMA_C_REAL(temp);
                TESTING_MALLOC_CPU( work, magmaFloatComplex, lwork );

                init_matrix( nopiv, N, N, h_A, lda );
                cpu_time = magma_wtime();
                lapackf77_chetrf((upper ? MagmaUpperStr: MagmaLowerStr), &N, h_A, &lda, ipiv, work, &lwork, &info);
                cpu_time = magma_wtime() - cpu_time;
                cpu_perf = gflops / cpu_time;
                if (info != 0)
                    printf("lapackf77_chetrf returned error %d: %s.\n",
                           (int) info, magma_strerror( info ));
                error_lapack = get_residual( nopiv, uplo, N, h_A, lda, ipiv );

                TESTING_FREE_CPU( work );
            }
           
            /* ====================================================================
               Performs operation using MAGMA
               =================================================================== */
            init_matrix( (nopiv | nopiv_gpu), N, N, h_A, lda );

            if (nopiv) {
                // CPU-interface to non-piv LDLt
                magma_setdevice(0);
                gpu_time = magma_wtime();
                magma_chetrf_nopiv( uplo, N, h_A, lda, &info);
                gpu_time = magma_wtime() - gpu_time;
            } else if (cpu) {
                // CPU-interface to Bunch-Kauffman LDLt
                magma_setdevice(0);
                gpu_time = magma_wtime();
                magma_chetrf( uplo, N, h_A, lda, ipiv, &info);
                gpu_time = magma_wtime() - gpu_time;
            } else if (nopiv_gpu) {
                // GPU-interface to non-piv LDLt
                magma_setdevice(0);
                magma_int_t ldda = 32*((N+31)/32);
                magmaFloatComplex *d_A;
                if (MAGMA_SUCCESS != magma_cmalloc( &d_A, N*ldda  )) {
                    printf( " failed to allocate d_A(%dx%d)\n",N,ldda);
                    return 0;
                }
                magma_csetmatrix(N, N, h_A, lda, d_A, ldda);
                gpu_time = magma_wtime();
                magma_chetrf_nopiv_gpu( uplo, N, d_A, ldda, &info);
                gpu_time = magma_wtime() - gpu_time;
                magma_cgetmatrix(N, N, d_A, ldda, h_A, lda);
                magma_free( d_A );
            } else if (row) {
                magma_setdevice(0);
                //magma_chetrf_gpu_row( uplo, N, h_A, lda, ipiv, work, lwork, &info);
            } else {
                magma_setdevice(0);
                //magma_chetrf_hybrid( uplo, N, h_A, lda, ipiv, work, lwork, &info);
            }
            gpu_perf = gflops / gpu_time;
            if (info != 0)
                printf("magma_chetrf returned error %d: %s.\n",
                       (int) info, magma_strerror( info ));
            
            /* =====================================================================
               Check the factorization
               =================================================================== */
            if ( opts.lapack ) {
                printf("%5d %5d   %7.2f (%7.2f)   %7.2f (%7.2f)",
                       (int) N, (int) N, cpu_perf, cpu_time, gpu_perf, gpu_time );
            }
            else {
                printf("%5d %5d     ---   (  ---  )   %7.2f (%7.2f)",
                       (int) N, (int) N, gpu_perf, gpu_time );
            }
            if ( opts.check == 2 ) {
                error = get_residual( (nopiv | nopiv_gpu), uplo, N, h_A, lda, ipiv );
                printf("   %8.2e   %s", error, (error < tol ? "ok" : "failed"));
                if (opts.lapack)
                    printf(" (lapack rel.res. = %8.2e)", error_lapack);
                printf("\n");
                status += ! (error < tol);
            }
            else if ( opts.check ) {
                error = get_LDLt_error( (nopiv | nopiv_gpu), uplo, N, h_A, lda, ipiv );
                printf("   %8.2e   %s\n", error, (error < tol ? "ok" : "failed"));
                status += ! (error < tol);
            }
            else {
                printf("     ---   \n");
            }
            
            TESTING_FREE_CPU( ipiv );
            TESTING_FREE_PIN( h_A  );
            fflush( stdout );
        }
        if ( opts.niter > 1 ) {
            printf( "\n" );
        }
    }

    TESTING_FINALIZE();
    return status;
}
