/*
    -- MAGMA (version 1.4.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       December 2013

       @generated c Tue Dec 17 13:18:56 2013

*/

// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <cuda_runtime_api.h>
#include <cublas.h>
#include <cblas.h>

// includes, project
#include "flops.h"
#include "magma.h"
#include "magma_lapack.h"
#include "testings.h"

#define PRECISION_c
#define COMPLEX

/* ////////////////////////////////////////////////////////////////////////////
   -- Testing cgeev
*/
int main( int argc, char** argv)
{
    TESTING_INIT();

    real_Double_t   gpu_time, cpu_time;
    magmaFloatComplex *h_A, *h_R, *VL, *VR, *h_work, *w1, *w2;
    magmaFloatComplex  c_neg_one = MAGMA_C_NEG_ONE;
    float *rwork;
    float tnrm, result[8];
    magma_int_t N, n2, lda, nb, lwork, info;
    magma_int_t ione     = 1;
    magma_int_t ISEED[4] = {0,0,0,1};
    magma_int_t status = 0;

    magma_opts opts;
    parse_opts( argc, argv, &opts );
    
    // need slightly looser bound (60*eps instead of 30*eps) for some tests
    opts.tolerance = max( 60., opts.tolerance );
    float tol = opts.tolerance * lapackf77_slamch("E");
    
    printf("    N   CPU Time (sec)   GPU Time (sec)   ||R||_F / ||A||_F\n");
    printf("===========================================================\n");
    for( int i = 0; i < opts.ntest; ++i ) {
        for( int iter = 0; iter < opts.niter; ++iter ) {
            N = opts.nsize[i];
            lda   = N;
            n2    = lda*N;
            nb    = magma_get_cgehrd_nb(N);
            lwork = N*(1 + nb);
            // generous workspace - required by cget22
            lwork = max( lwork, N*(5 + 2*N) );
            
            TESTING_MALLOC_CPU( w1,     magmaFloatComplex, N );
            TESTING_MALLOC_CPU( w2,     magmaFloatComplex, N );
            TESTING_MALLOC_CPU( rwork,  float, 2*N );
            TESTING_MALLOC_CPU( h_A,    magmaFloatComplex, n2 );
            
            TESTING_MALLOC_PIN( h_R,    magmaFloatComplex, n2 );
            TESTING_MALLOC_PIN( VL,     magmaFloatComplex, n2 );
            TESTING_MALLOC_PIN( VR,     magmaFloatComplex, n2 );
            TESTING_MALLOC_PIN( h_work, magmaFloatComplex, lwork );
            
            /* Initialize the matrix */
            lapackf77_clarnv( &ione, ISEED, &n2, h_A );
            lapackf77_clacpy( MagmaUpperLowerStr, &N, &N, h_A, &lda, h_R, &lda );
            
            /* ====================================================================
               Performs operation using MAGMA
               =================================================================== */
            gpu_time = magma_wtime();
            magma_cgeev_m( opts.jobvl, opts.jobvr,
                           N, h_R, lda, w1,
                           VL, lda, VR, lda,
                           h_work, lwork, rwork, &info );
            gpu_time = magma_wtime() - gpu_time;
            if (info != 0)
                printf("magma_cgeev returned error %d: %s.\n",
                       (int) info, magma_strerror( info ));
            
            /* =====================================================================
               Performs operation using LAPACK
               =================================================================== */
            if ( opts.lapack ) {
                cpu_time = magma_wtime();
                lapackf77_cgeev( &opts.jobvl, &opts.jobvr,
                                 &N, h_A, &lda, w2,
                                 VL, &lda, VR, &lda,
                                 h_work, &lwork, rwork, &info );
                cpu_time = magma_wtime() - cpu_time;
                if (info != 0)
                    printf("lapackf77_cgeev returned error %d: %s.\n",
                           (int) info, magma_strerror( info ));
                
                printf("%5d   %7.2f          %7.2f\n",
                       (int) N, cpu_time, gpu_time);
            }
            else {
                printf("%5d     ---            %7.2f\n",
                       (int) N, gpu_time);
            }
            
            /* =====================================================================
               Check the result
               =================================================================== */
            if ( opts.check ) {
                /* ===================================================================
                 * Check the result following LAPACK's [zcds]drvev routine.
                 * The following 7 tests are performed:
                 *     (1)     | A * VR - VR * W | / ( n |A| )
                 *
                 *       Here VR is the matrix of unit right eigenvectors.
                 *       W is a diagonal matrix with diagonal entries W(j).
                 *
                 *     (2)     | A**H * VL - VL * W**H | / ( n |A| )
                 *
                 *       Here VL is the matrix of unit left eigenvectors, A**H is the
                 *       conjugate-transpose of A, and W is as above.
                 *
                 *     (3)     | |VR(i)| - 1 |   and whether largest component real
                 *
                 *       VR(i) denotes the i-th column of VR.
                 *
                 *     (4)     | |VL(i)| - 1 |   and whether largest component real
                 *
                 *       VL(i) denotes the i-th column of VL.
                 *
                 *     (5)     W(full) = W(partial)
                 *
                 *       W(full) denotes the eigenvalues computed when both VR and VL
                 *       are also computed, and W(partial) denotes the eigenvalues
                 *       computed when only W, only W and VR, or only W and VL are
                 *       computed.
                 *
                 *     (6)     VR(full) = VR(partial)
                 *
                 *       VR(full) denotes the right eigenvectors computed when both VR
                 *       and VL are computed, and VR(partial) denotes the result
                 *       when only VR is computed.
                 *
                 *     (7)     VL(full) = VL(partial)
                 *
                 *       VL(full) denotes the left eigenvectors computed when both VR
                 *       and VL are also computed, and VL(partial) denotes the result
                 *       when only VL is computed.
                 ================================================================= */
                float ulp, ulpinv, vmx, vrmx, vtst;
                magmaFloatComplex *LRE, DUM;
                TESTING_MALLOC_PIN( LRE, magmaFloatComplex, n2 );
                
                ulp = lapackf77_slamch( "P" );
                ulpinv = 1./ulp;
                
                // Initialize RESULT
                for( int j = 0; j < 8; ++j )
                    result[j] = -1.;
                
                lapackf77_clarnv( &ione, ISEED, &n2, h_A );
                lapackf77_clacpy( MagmaUpperLowerStr, &N, &N, h_A, &lda, h_R, &lda );
                
                // ----------
                // Compute eigenvalues, left and right eigenvectors, and test them
                magma_cgeev_m( MagmaVec, MagmaVec,
                               N, h_R, lda, w1,
                               VL, lda, VR, lda,
                               h_work, lwork, rwork, &info );
                
                // Do test 1
                lapackf77_cget22( MagmaNoTransStr, MagmaNoTransStr, MagmaNoTransStr,
                                  &N, h_A, &lda, VR, &lda, w1,
                                  h_work, rwork, &result[0] );
                result[0] *= ulp;
                
                // Do test 2
                lapackf77_cget22( MagmaConjTransStr, MagmaNoTransStr, MagmaConjTransStr,
                                  &N, h_A, &lda, VL, &lda, w1,
                                  h_work, rwork, &result[1] );
                result[1] *= ulp;
                
                // Do test 3
                result[2] = -1.;
                for( int j = 0; j < N; ++j ) {
                    tnrm = cblas_scnrm2(N, &VR[j*lda], ione);
                    result[2] = max( result[2], min( ulpinv, fabs(tnrm-1.)/ulp ));
                    
                    vmx  = vrmx = 0.;
                    for( int jj = 0; jj <N; ++jj ) {
                        vtst = cuCabsf(VR[jj + j*lda]);
                        if (vtst > vmx)
                            vmx = vtst;
                        
                        if (MAGMA_C_IMAG(VR[jj + j*lda])==0. &&
                            fabs( MAGMA_C_REAL(VR[jj+j*lda]) ) > vrmx)
                        {
                            vrmx = fabs( MAGMA_C_REAL( VR[jj+j*lda] ) );
                        }
                    }
                    if (vrmx / vmx < 1. - ulp*2.)
                        result[2] = ulpinv;
                }
                result[2] *= ulp;
                
                // Do test 4
                result[3] = -1.;
                for( int j = 0; j < N; ++j ) {
                    tnrm = cblas_scnrm2(N, &VL[j*lda], ione);
                    result[3] = max( result[3], min( ulpinv, fabs(tnrm - 1.)/ ulp ));
                    
                    vmx = vrmx = 0.;
                    for( int jj = 0; jj < N; ++jj ) {
                        vtst = cuCabsf(VL[jj + j*lda]);
                        if (vtst > vmx)
                            vmx = vtst;
                        
                        if (MAGMA_C_IMAG(VL[jj + j*lda])==0. &&
                            fabs( MAGMA_C_REAL( VL[jj + j*lda] ) ) > vrmx)
                        {
                            vrmx = fabs( MAGMA_C_REAL( VL[jj+j*lda]) );
                        }
                    }
                    if (vrmx / vmx < 1. - ulp*2.)
                        result[3] = ulpinv;
                }
                result[3] *= ulp;
                
                // ----------
                // Compute eigenvalues only, and test them
                lapackf77_clacpy( MagmaUpperLowerStr, &N, &N, h_A, &lda, h_R, &lda );
                magma_cgeev_m( MagmaNoVec, MagmaNoVec,
                               N, h_R, lda, w2,
                               &DUM, 1, &DUM, 1,
                               h_work, lwork, rwork, &info );
                
                if (info != 0) {
                    result[0] = ulpinv;
                    printf("magma_cgeev (case N, N) returned error %d: %s.\n",
                           (int) info, magma_strerror( info ));
                }
                
                // Do test 5
                result[4] = 1;
                for( int j = 0; j < N; ++j )
                    if ( ! MAGMA_C_EQUAL( w1[j], w2[j] ))
                        result[4] = 0;
                //if (result[4] == 0) printf("test 5 failed with N N\n");
                
                // ----------
                // Compute eigenvalues and right eigenvectors, and test them
                lapackf77_clacpy( MagmaUpperLowerStr, &N, &N, h_A, &lda, h_R, &lda );
                magma_cgeev_m( MagmaNoVec, MagmaVec,
                               N, h_R, lda, w2,
                               &DUM, 1, LRE, lda,
                               h_work, lwork, rwork, &info );
                
                if (info != 0) {
                    result[0] = ulpinv;
                    printf("magma_cgeev (case N, V) returned error %d: %s.\n",
                           (int) info, magma_strerror( info ));
                }
                
                // Do test 5 again
                result[4] = 1;
                for( int j = 0; j < N; ++j )
                    if ( ! MAGMA_C_EQUAL( w1[j], w2[j] ))
                        result[4] = 0;
                //if (result[4] == 0) printf("test 5 failed with N V\n");
                
                // Do test 6
                result[5] = 1;
                for( int j = 0; j < N; ++j )
                    for( int jj = 0; jj < N; ++jj )
                        if ( ! MAGMA_C_EQUAL( VR[j+jj*lda], LRE[j+jj*lda] ))
                            result[5] = 0;
                
                // ----------
                // Compute eigenvalues and left eigenvectors, and test them
                lapackf77_clacpy( MagmaUpperLowerStr, &N, &N, h_A, &lda, h_R, &lda );
                magma_cgeev_m( MagmaVec, MagmaNoVec,
                               N, h_R, lda, w2,
                               LRE, lda, &DUM, 1,
                               h_work, lwork, rwork, &info );
                
                if (info != 0) {
                    result[0] = ulpinv;
                    printf("magma_cgeev (case V, N) returned error %d: %s.\n",
                           (int) info, magma_strerror( info ));
                }
                
                // Do test 5 again
                result[4] = 1;
                for( int j = 0; j < N; ++j )
                    if ( ! MAGMA_C_EQUAL( w1[j], w2[j] ))
                        result[4] = 0;
                //if (result[4] == 0) printf("test 5 failed with V N\n");
                
                // Do test 7
                result[6] = 1;
                for( int j = 0; j < N; ++j )
                    for( int jj = 0; jj < N; ++jj )
                        if ( ! MAGMA_C_EQUAL( VL[j+jj*lda], LRE[j+jj*lda] ))
                            result[6] = 0;
                
                printf("Test 1: | A * VR - VR * W | / ( n |A| ) = %8.2e%s\n", result[0], (result[0] < tol ? "  ok" : "  failed"));
                printf("Test 2: | A'* VL - VL * W'| / ( n |A| ) = %8.2e%s\n", result[1], (result[1] < tol ? "  ok" : "  failed"));
                printf("Test 3: |  |VR(i)| - 1    |             = %8.2e%s\n", result[2], (result[2] < tol ? "  ok" : "  failed"));
                printf("Test 4: |  |VL(i)| - 1    |             = %8.2e%s\n", result[3], (result[3] < tol ? "  ok" : "  failed"));
                printf("Test 5:   W (full)  ==  W (partial)     = %s\n",                 (result[4] == 1. ? "  ok" : "  failed"));
                printf("Test 6:  VR (full)  == VR (partial)     = %s\n",                 (result[5] == 1. ? "  ok" : "  failed"));
                printf("Test 7:  VL (full)  == VL (partial)     = %s\n\n",               (result[6] == 1. ? "  ok" : "  failed"));
                status |= ! (result[0] < tol);
                status |= ! (result[1] < tol);
                status |= ! (result[2] < tol);
                status |= ! (result[3] < tol);
                status |= ! (result[4] == 1.);
                status |= ! (result[5] == 1.);
                status |= ! (result[6] == 1.);
                
                TESTING_FREE_PIN( LRE );
            }
            
            TESTING_FREE_CPU( w1    );
            TESTING_FREE_CPU( w2    );
            TESTING_FREE_CPU( rwork );
            TESTING_FREE_CPU( h_A   );
            
            TESTING_FREE_PIN( h_R    );
            TESTING_FREE_PIN( VL     );
            TESTING_FREE_PIN( VR     );
            TESTING_FREE_PIN( h_work );
        }
        if ( opts.niter > 1 ) {
            printf( "\n" );
        }
    }

    TESTING_FINALIZE();
    return status;
}
