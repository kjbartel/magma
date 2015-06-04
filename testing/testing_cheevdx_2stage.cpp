/*
    -- MAGMA (version 1.4.0-beta2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2010

    @author Raffaele Solca
    @author Azzam Haidar

    @generated c Fri Jun 28 19:34:03 2013

*/

// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <cuda_runtime_api.h>
#include <cublas.h>

// includes, project
#include "magma.h"
#include "magma_lapack.h"
#include "testings.h"
#include "magma_cbulge.h"

#define PRECISION_c
#define absv(v1) ((v1)>0? (v1): -(v1))


static magma_int_t check_orthogonality(magma_int_t, magma_int_t, magmaFloatComplex*, magma_int_t, float);
static magma_int_t check_reduction(magma_int_t, magma_int_t, magma_int_t, magmaFloatComplex*, float*, magma_int_t, magmaFloatComplex*, float);
static magma_int_t check_solution(magma_int_t, float*, float*, float);
/* ////////////////////////////////////////////////////////////////////////////
   -- Testing chegvdx
*/
int main( int argc, char** argv)
{

    TESTING_INIT();

    real_Double_t gpu_time;

    magmaFloatComplex *h_A, *h_R, *h_B, *h_S, *h_work;

#if defined(PRECISION_z) || defined(PRECISION_c)
    float *rwork;
    magma_int_t lrwork;
#endif

    /* Matrix size */
    float *w1, *w2, result[2];
    magma_int_t *iwork;
    magma_int_t N, n2, info, lwork, liwork;
    magmaFloatComplex c_zero    = MAGMA_C_ZERO;
    magmaFloatComplex c_one     = MAGMA_C_ONE;
    magmaFloatComplex c_neg_one = MAGMA_C_NEG_ONE;
    magma_int_t ione     = 1;
    magma_int_t ISEED[4] = {0,0,0,1};;
    magma_int_t info_ortho     = 0;
    magma_int_t info_solution  = 0;
    magma_int_t info_reduction = 0;

    magma_timestr_t start, end;

    magma_opts opts;
    parse_opts( argc, argv, &opts );

    magma_int_t ngpu = opts.ngpu;
    char jobz = opts.jobz;
    magma_int_t checkres = opts.check;

    char range = 'A';
    char uplo = opts.uplo;
    magma_int_t itype = opts.itype;
    float f = opts.fraction;

    if (f != 1)
        range='I';

    if ( checkres && jobz == MagmaNoVec ) {
        fprintf( stderr, "checking results requires vectors; setting jobz=V (option -JV)\n" );
        jobz = MagmaVec;
    }

    printf("using: itype = %d, jobz = %c,range = %c, uplo = %c, checkres = %d, fraction = %6.4f\n", itype, jobz, range, uplo, checkres, f);

    printf("  N     M     GPU Time(s) \n");
    printf("==========================\n");
    for( magma_int_t i = 0; i < opts.ntest; ++i ) {
        for( magma_int_t iter = 0; iter < opts.niter; ++iter ) {
            N = opts.nsize[i];
            n2     = N*N;
#if defined(PRECISION_z) || defined(PRECISION_c)
            lwork  = magma_cbulge_get_lq2(N) + 2*N + N*N;
            lrwork = 1 + 5*N +2*N*N;
#else
            lwork  = magma_cbulge_get_lq2(N) + 1 + 6*N + 2*N*N;
#endif
            liwork = 3 + 5*N;

            /* Allocate host memory for the matrix */
            TESTING_MALLOC(   h_A, magmaFloatComplex, n2);
            TESTING_MALLOC(    w1, float         ,  N);
            TESTING_MALLOC(    w2, float         ,  N);
            TESTING_HOSTALLOC(h_R, magmaFloatComplex, n2);
            TESTING_HOSTALLOC(h_work, magmaFloatComplex,  lwork);
#if defined(PRECISION_z) || defined(PRECISION_c)
            TESTING_HOSTALLOC( rwork,          float, lrwork);
#endif
            TESTING_MALLOC(    iwork,     magma_int_t, liwork);


            /* Initialize the matrix */
            lapackf77_clarnv( &ione, ISEED, &n2, h_A );
            /* Make diagonal real */
            for(int i=0; i<N; i++) {
                MAGMA_C_SET2REAL( h_A[i*N+i], MAGMA_C_REAL(h_A[i*N+i]) );
            }


            magma_int_t m1 = 0;
            float vl = 0;
            float vu = 0;
            magma_int_t il = 0;
            magma_int_t iu = 0;
            if (range == 'I'){
                il = 1;
                iu = f*N;
            }

            if(opts.warmup){
                // ==================================================================
                // Warmup using MAGMA
                // ==================================================================
                lapackf77_clacpy( MagmaUpperLowerStr, &N, &N, h_A, &N, h_R, &N );
                if(ngpu==1){
                    printf("calling cheevdx_2stage 1 GPU\n");
                    magma_cheevdx_2stage(jobz, range, uplo, N, 
                                    h_R, N, 
                                    vl, vu, il, iu, 
                                    &m1, w1, 
                                    h_work, lwork, 
#if defined(PRECISION_z) || defined(PRECISION_c)
                                    rwork, lrwork, 
#endif         
                                    iwork, liwork, 
                                    &info);
               
                }else{
                    printf("calling cheevdx_2stage_m %d GPU\n", ngpu);
                    magma_cheevdx_2stage_m(ngpu, jobz, range, uplo, N, 
                                    h_R, N, 
                                    vl, vu, il, iu, 
                                    &m1, w1, 
                                    h_work, lwork, 
#if defined(PRECISION_z) || defined(PRECISION_c)
                                    rwork, lrwork, 
#endif         
                                    iwork, liwork, 
                                    &info);
                }
            }


            // ===================================================================
            // Performs operation using MAGMA
            // ===================================================================
            lapackf77_clacpy( MagmaUpperLowerStr, &N, &N, h_A, &N, h_R, &N );
            start = get_current_time();
            if(ngpu==1){
                printf("calling cheevdx_2stage 1 GPU\n");
                magma_cheevdx_2stage(jobz, range, uplo, N, 
                                h_R, N, 
                                vl, vu, il, iu, 
                                &m1, w1, 
                                h_work, lwork, 
#if defined(PRECISION_z) || defined(PRECISION_c)
                                rwork, lrwork, 
#endif     
                                iwork, liwork, 
                                &info);
           
            }else{
                printf("calling cheevdx_2stage_m %d GPU\n", ngpu);
                magma_cheevdx_2stage_m(ngpu, jobz, range, uplo, N, 
                                h_R, N, 
                                vl, vu, il, iu, 
                                &m1, w1, 
                                h_work, lwork, 
#if defined(PRECISION_z) || defined(PRECISION_c)
                                rwork, lrwork, 
#endif     
                                iwork, liwork, 
                                &info);
            }
            end = get_current_time();
            gpu_time = GetTimerValue(start,end)/1000.;

            if ( checkres ) {
                float eps   = lapackf77_slamch("E");
                printf("\n");
                printf("------ TESTS FOR MAGMA CHEEVD ROUTINE -------  \n");
                printf("        Size of the Matrix %d by %d\n", N, N);
                printf("\n");
                printf(" The matrix A is randomly generated for each test.\n");
                printf("============\n");
                printf(" The relative machine precision (eps) is to be %e \n",eps);
                printf(" Computational tests pass if scaled residuals are less than 60.\n");
              
                /* Check the orthogonality, reduction and the eigen solutions */
                if (jobz == MagmaVec) {
                    info_ortho = check_orthogonality(N, N, h_R, N, eps);
                    info_reduction = check_reduction(uplo, N, 1, h_A, w1, N, h_R, eps);
                }
                printf("------ CALLING LAPACK CHEEVD TO COMPUTE only eigenvalue and verify elementswise -------  \n");
                lapackf77_cheevd("N", "L", &N, 
                                h_A, &N, w2, 
                                h_work, &lwork, 
#if defined(PRECISION_z) || defined(PRECISION_c)
                                rwork, &lrwork, 
#endif     
                                iwork, &liwork, 
                                &info);
                info_solution = check_solution(N, w2, w1, eps);
              
                if ( (info_solution == 0) & (info_ortho == 0) & (info_reduction == 0) ) {
                    printf("***************************************************\n");
                    printf(" ---- TESTING CHEEVD ...................... PASSED !\n");
                    printf("***************************************************\n");
                }
                else {
                    printf("************************************************\n");
                    printf(" - TESTING CHEEVD ... FAILED !\n");
                    printf("************************************************\n");
                }
            }


            /* =====================================================================
             Print execution time
             =================================================================== */
            printf("%5d %5d     %6.2f\n",
                   (magma_int_t) N, (magma_int_t) m1, gpu_time);

            TESTING_FREE(       h_A);
            TESTING_FREE(        w1);
            TESTING_FREE(        w2);
#if defined(PRECISION_z) || defined(PRECISION_c)
            TESTING_HOSTFREE( rwork);
#endif
            TESTING_FREE(     iwork);
            TESTING_HOSTFREE(h_work);
            TESTING_HOSTFREE(   h_R);
        }
        if ( opts.niter > 1 ) {
            printf( "\n" );
        }
    }

    /* Shutdown */
    TESTING_FINALIZE();

    return 0;
}



/*-------------------------------------------------------------------
 * Check the orthogonality of Q
 */
static magma_int_t check_orthogonality(magma_int_t M, magma_int_t N, magmaFloatComplex *Q, magma_int_t LDQ, float eps)
{
    float  done  =  1.0;
    float  mdone = -1.0;
    magmaFloatComplex c_zero    = MAGMA_C_ZERO;
    magmaFloatComplex c_one     = MAGMA_C_ONE;
    float  normQ, result;
    magma_int_t     info_ortho;
    magma_int_t     minMN = min(M, N);
    float *work = (float *)malloc(minMN*sizeof(float));

    /* Build the idendity matrix */
    magmaFloatComplex *Id = (magmaFloatComplex *) malloc(minMN*minMN*sizeof(magmaFloatComplex));
    lapackf77_claset("A", &minMN, &minMN, &c_zero, &c_one, Id, &minMN);

    /* Perform Id - Q'Q */
    if (M >= N)
        blasf77_cherk("U", "C", &N, &M, &done, Q, &LDQ, &mdone, Id, &N);
    else
        blasf77_cherk("U", "N", &M, &N, &done, Q, &LDQ, &mdone, Id, &M);

    normQ = lapackf77_clanhe("I", "U", &minMN, Id, &minMN, work);

    result = normQ / (minMN * eps);
    printf(" ======================================================\n");
    printf(" ||Id-Q'*Q||_oo / (minMN*eps)          : %15.3E \n",  result );
    printf(" ======================================================\n");

    if ( isnan(result) || isinf(result) || (result > 60.0) ) {
        printf("-- Orthogonality is suspicious ! \n");
        info_ortho=1;
    }
    else {
        printf("-- Orthogonality is CORRECT ! \n");
        info_ortho=0;
    }
    free(work); free(Id);
    return info_ortho;
}
/*------------------------------------------------------------
 *  Check the reduction 
 */
static magma_int_t check_reduction(magma_int_t uplo, magma_int_t N, magma_int_t bw, magmaFloatComplex *A, float *D, magma_int_t LDA, magmaFloatComplex *Q, float eps )
{
    magmaFloatComplex c_one     = MAGMA_C_ONE;
    magmaFloatComplex c_neg_one = MAGMA_C_NEG_ONE;
    magmaFloatComplex *TEMP     = (magmaFloatComplex *)malloc(N*N*sizeof(magmaFloatComplex));
    magmaFloatComplex *Residual = (magmaFloatComplex *)malloc(N*N*sizeof(magmaFloatComplex));
    float *work = (float *)malloc(N*sizeof(float));
    float Anorm, Rnorm, result;
    magma_int_t info_reduction;
    magma_int_t i;
    magma_int_t ione=1;
    char luplo =  uplo == MagmaLower ? 'L' : 'U';

    /* Compute TEMP =  Q * LAMBDA */
    lapackf77_clacpy("A", &N, &N, Q, &LDA, TEMP, &N);        
    for (i = 0; i < N; i++){
            blasf77_csscal(&N, &D[i], &(TEMP[i*N]), &ione);
    }
    /* Compute Residual = A - Q * LAMBDA * Q^H */
    /* A is Hermetian but both upper and lower 
     * are assumed valable here for checking 
     * otherwise it need to be symetrized before 
     * checking.
     */ 
    lapackf77_clacpy("A", &N, &N, A, &LDA, Residual, &N);        
    blasf77_cgemm("N", "C", &N, &N, &N, &c_neg_one, TEMP, &N, Q, &LDA, &c_one, Residual,     &N);

    Rnorm = lapackf77_clange("1", &N, &N, Residual, &N,   work);
    Anorm = lapackf77_clanhe("1", &luplo, &N, A,        &LDA, work);

    result = Rnorm / ( Anorm * N * eps);
    if ( uplo == MagmaLower ){
        printf(" ======================================================\n");
        printf(" ||A-Q*LAMBDA*Q'||_oo/(||A||_oo.N.eps) : %15.3E \n",  result );
        printf(" ======================================================\n");
    }else{ 
        printf(" ======================================================\n");
        printf(" ||A-Q'*LAMBDA*Q||_oo/(||A||_oo.N.eps) : %15.3E \n",  result );
        printf(" ======================================================\n");
    }

    if ( isnan(result) || isinf(result) || (result > 60.0) ) {
        printf("-- Reduction is suspicious ! \n");
        info_reduction = 1;
    }
    else {
        printf("-- Reduction is CORRECT ! \n");
        info_reduction = 0;
    }

    free(TEMP); free(Residual);
    free(work);

    return info_reduction;
}
/*------------------------------------------------------------
 *  Check the eigenvalues 
 */
static magma_int_t check_solution(magma_int_t N, float *E1, float *E2, float eps)
{
    magma_int_t info_solution, i;
    float resid;
    float maxtmp;
    float maxel = fabs( fabs(E1[0]) - fabs(E2[0]) );
    float maxeig = max( fabs(E1[0]), fabs(E2[0]) );
    for (i = 1; i < N; i++){
        resid   = fabs(fabs(E1[i])-fabs(E2[i]));
        maxtmp  = max(fabs(E1[i]), fabs(E2[i]));

        /* Update */
        maxeig = max(maxtmp, maxeig);
        maxel  = max(resid,  maxel );
    }

    maxel = maxel / (maxeig * N * eps);
    printf(" ======================================================\n");
    printf(" | D - eigcomputed | / (|D| * N * eps) : %15.3E \n",  maxel );
    printf(" ======================================================\n");

    if ( isnan(maxel) || isinf(maxel) || (maxel > 100) ) {
        printf("-- The eigenvalues are suspicious ! \n");
        info_solution = 1;
    }
    else{
        printf("-- The eigenvalues are CORRECT ! \n");
        info_solution = 0;
    }
    return info_solution;
}

