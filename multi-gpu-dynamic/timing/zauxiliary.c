/**
 *
 * @precisions normal z -> c d s
 *
 **/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <cblas.h>
#include <lapacke.h>
#include <plasma.h>
#include <core_blas.h>
#include "auxiliary.h"

/*-------------------------------------------------------------------
 * Check the orthogonality of Q
 */

int zcheck_orthogonality(int M, int N, int LDQ, PLASMA_Complex64_t *Q)
{
    double alpha, beta;
    double normQ;
    int info_ortho;
    int i;
    int minMN = min(M, N);
    double eps;
    double *work = (double *)malloc(minMN*sizeof(double));

    eps = LAPACKE_dlamch_work('e');
    alpha = 1.0;
    beta  = -1.0;

    /* Build the idendity matrix USE DLASET?*/
    PLASMA_Complex64_t *Id = (PLASMA_Complex64_t *) malloc(minMN*minMN*sizeof(PLASMA_Complex64_t));
    memset((void*)Id, 0, minMN*minMN*sizeof(PLASMA_Complex64_t));
    for (i = 0; i < minMN; i++)
        Id[i*minMN+i] = (PLASMA_Complex64_t)1.0;

    /* Perform Id - Q'Q */
    if (M >= N)
        cblas_zherk(CblasColMajor, CblasUpper, CblasConjTrans, N, M, alpha, Q, LDQ, beta, Id, N);
    else
        cblas_zherk(CblasColMajor, CblasUpper, CblasNoTrans, M, N, alpha, Q, LDQ, beta, Id, M);

    normQ = LAPACKE_zlansy_work(LAPACK_COL_MAJOR, lapack_const(PlasmaInfNorm), 'u', minMN, Id, minMN, work);

    printf("============\n");
    printf("Checking the orthogonality of Q \n");
    printf("||Id-Q'*Q||_oo / (N*eps) = %e \n",normQ/(minMN*eps));

    if ( isnan(normQ / (minMN * eps)) || (normQ / (minMN * eps) > 10.0) ) {
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
 *  Check the factorization QR
 */

int zcheck_QRfactorization(int M, int N, PLASMA_Complex64_t *A1, PLASMA_Complex64_t *A2, int LDA, PLASMA_Complex64_t *Q)
{
    double Anorm, Rnorm;
    PLASMA_Complex64_t alpha, beta;
    int info_factorization;
    int i,j;
    double eps;

    eps = LAPACKE_dlamch_work('e');

    PLASMA_Complex64_t *Ql       = (PLASMA_Complex64_t *)malloc(M*N*sizeof(PLASMA_Complex64_t));
    PLASMA_Complex64_t *Residual = (PLASMA_Complex64_t *)malloc(M*N*sizeof(PLASMA_Complex64_t));
    double *work              = (double *)malloc(max(M,N)*sizeof(double));

    alpha=1.0;
    beta=0.0;

    if (M >= N) {
        /* Extract the R */
        PLASMA_Complex64_t *R = (PLASMA_Complex64_t *)malloc(N*N*sizeof(PLASMA_Complex64_t));
        memset((void*)R, 0, N*N*sizeof(PLASMA_Complex64_t));
        LAPACKE_zlacpy_work(LAPACK_COL_MAJOR,'u', M, N, A2, LDA, R, N);

        /* Perform Ql=Q*R */
        memset((void*)Ql, 0, M*N*sizeof(PLASMA_Complex64_t));
        cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, M, N, N, CBLAS_SADDR(alpha), Q, LDA, R, N, CBLAS_SADDR(beta), Ql, M);
        free(R);
    }
    else {
        /* Extract the L */
        PLASMA_Complex64_t *L = (PLASMA_Complex64_t *)malloc(M*M*sizeof(PLASMA_Complex64_t));
        memset((void*)L, 0, M*M*sizeof(PLASMA_Complex64_t));
        LAPACKE_zlacpy_work(LAPACK_COL_MAJOR,'l', M, N, A2, LDA, L, M);

    /* Perform Ql=LQ */
        memset((void*)Ql, 0, M*N*sizeof(PLASMA_Complex64_t));
        cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, M, N, M, CBLAS_SADDR(alpha), L, M, Q, LDA, CBLAS_SADDR(beta), Ql, M);
        free(L);
    }

    /* Compute the Residual */
    for (i = 0; i < M; i++)
        for (j = 0 ; j < N; j++)
            Residual[j*M+i] = A1[j*LDA+i]-Ql[j*M+i];

    Rnorm = LAPACKE_zlange_work(LAPACK_COL_MAJOR, lapack_const(PlasmaInfNorm), M, N, Residual, M, work);
    Anorm = LAPACKE_zlange_work(LAPACK_COL_MAJOR, lapack_const(PlasmaInfNorm), M, N, A2, LDA, work);

    if (M >= N) {
        printf("============\n");
        printf("Checking the QR Factorization \n");
        printf("-- ||A-QR||_oo/(||A||_oo.N.eps) = %e \n",Rnorm/(Anorm*N*eps));
    }
    else {
        printf("============\n");
        printf("Checking the LQ Factorization \n");
        printf("-- ||A-LQ||_oo/(||A||_oo.N.eps) = %e \n",Rnorm/(Anorm*N*eps));
    }

    if (isnan(Rnorm / (Anorm * N *eps)) || (Rnorm / (Anorm * N * eps) > 10.0) ) {
        printf("-- Factorization is suspicious ! \n");
        info_factorization = 1;
    }
    else {
        printf("-- Factorization is CORRECT ! \n");
        info_factorization = 0;
    }

    free(work); free(Ql); free(Residual);

    return info_factorization;
}

/*------------------------------------------------------------------------
 *  Check the factorization of the matrix A2
 */

int zcheck_LLTfactorization(int N, PLASMA_Complex64_t *A1, PLASMA_Complex64_t *A2, int LDA, int uplo)
{
    double Anorm, Rnorm;
    PLASMA_Complex64_t alpha;
    int info_factorization;
    int i,j;
    double eps;

    eps = LAPACKE_dlamch_work('e');

    PLASMA_Complex64_t *Residual = (PLASMA_Complex64_t *)malloc(N*N*sizeof(PLASMA_Complex64_t));
    PLASMA_Complex64_t *L1       = (PLASMA_Complex64_t *)malloc(N*N*sizeof(PLASMA_Complex64_t));
    PLASMA_Complex64_t *L2       = (PLASMA_Complex64_t *)malloc(N*N*sizeof(PLASMA_Complex64_t));
    double *work              = (double *)malloc(N*sizeof(double));

    memset((void*)L1, 0, N*N*sizeof(PLASMA_Complex64_t));
    memset((void*)L2, 0, N*N*sizeof(PLASMA_Complex64_t));

    alpha= 1.0;

    LAPACKE_zlacpy_work(LAPACK_COL_MAJOR,' ', N, N, A1, LDA, Residual, N);

    /* Dealing with L'L or U'U  */
    if (uplo == PlasmaUpper){
        LAPACKE_zlacpy_work(LAPACK_COL_MAJOR,'u', N, N, A2, LDA, L1, N);
        LAPACKE_zlacpy_work(LAPACK_COL_MAJOR,'u', N, N, A2, LDA, L2, N);
        cblas_ztrmm(CblasColMajor, CblasLeft, CblasUpper, CblasConjTrans, CblasNonUnit, N, N, CBLAS_SADDR(alpha), L1, N, L2, N);
    }
    else{
        LAPACKE_zlacpy_work(LAPACK_COL_MAJOR,'l', N, N, A2, LDA, L1, N);
        LAPACKE_zlacpy_work(LAPACK_COL_MAJOR,'l', N, N, A2, LDA, L2, N);
        cblas_ztrmm(CblasColMajor, CblasRight, CblasLower, CblasConjTrans, CblasNonUnit, N, N, CBLAS_SADDR(alpha), L1, N, L2, N);
    }

    /* Compute the Residual || A -L'L|| */
    for (i = 0; i < N; i++)
        for (j = 0; j < N; j++)
           Residual[j*N+i] = L2[j*N+i] - Residual[j*N+i];

    Rnorm = LAPACKE_zlange_work(LAPACK_COL_MAJOR, lapack_const(PlasmaInfNorm), N, N, Residual, N, work);
    Anorm = LAPACKE_zlange_work(LAPACK_COL_MAJOR, lapack_const(PlasmaInfNorm), N, N, A1, LDA, work);

    printf("============\n");
    printf("Checking the Cholesky Factorization \n");
    printf("-- ||L'L-A||_oo/(||A||_oo.N.eps) = %e \n",Rnorm/(Anorm*N*eps));

    if ( isnan(Rnorm/(Anorm*N*eps)) || (Rnorm/(Anorm*N*eps) > 10.0) ){
        printf("-- Factorization is suspicious ! \n");
        info_factorization = 1;
    }
    else{
        printf("-- Factorization is CORRECT ! \n");
        info_factorization = 0;
    }

    free(Residual); free(L1); free(L2); free(work);

    return info_factorization;
}

/*--------------------------------------------------------------
 * Check the gemm
 */
double zcheck_gemm(PLASMA_enum transA, PLASMA_enum transB, int M, int N, int K,
                   PLASMA_Complex64_t alpha, PLASMA_Complex64_t *A, int LDA,
                   PLASMA_Complex64_t *B, int LDB,
                   PLASMA_Complex64_t beta, PLASMA_Complex64_t *Cplasma,
                   PLASMA_Complex64_t *Cref, int LDC,
                   double *Cinitnorm, double *Cplasmanorm, double *Clapacknorm )
{
    PLASMA_Complex64_t beta_const = -1.0;
    double Rnorm;
    double *work = (double *)malloc(max(K,max(M, N))* sizeof(double));

    *Cinitnorm   = LAPACKE_zlange_work(LAPACK_COL_MAJOR, lapack_const(PlasmaInfNorm), M, N, Cref,    LDC, work);
    *Cplasmanorm = LAPACKE_zlange_work(LAPACK_COL_MAJOR, lapack_const(PlasmaInfNorm), M, N, Cplasma, LDC, work);

    cblas_zgemm(CblasColMajor, (CBLAS_TRANSPOSE)transA, (CBLAS_TRANSPOSE)transB, M, N, K,
                CBLAS_SADDR(alpha), A, LDA, B, LDB, CBLAS_SADDR(beta), Cref, LDC);

    *Clapacknorm = LAPACKE_zlange_work(LAPACK_COL_MAJOR, lapack_const(PlasmaInfNorm), M, N, Cref, LDC, work);

    cblas_zaxpy(LDC * N, CBLAS_SADDR(beta_const), Cplasma, 1, Cref, 1);

    Rnorm = LAPACKE_zlange_work(LAPACK_COL_MAJOR, lapack_const(PlasmaInfNorm), M, N, Cref, LDC, work);

    free(work);

    return Rnorm;
}

/*--------------------------------------------------------------
 * Check the gemm
 */
double zcheck_trsm(PLASMA_enum side, PLASMA_enum uplo, PLASMA_enum trans, PLASMA_enum diag,
                   int M, int NRHS, PLASMA_Complex64_t alpha,
                   PLASMA_Complex64_t *A, int LDA,
                   PLASMA_Complex64_t *Bplasma, PLASMA_Complex64_t *Bref, int LDB,
                   double *Binitnorm, double *Bplasmanorm, double *Blapacknorm )
{
    PLASMA_Complex64_t beta_const = -1.0;
    double Rnorm;
    double *work = (double *)malloc(max(M, NRHS)* sizeof(double));
    double eps = LAPACKE_dlamch_work('e');

    *Binitnorm   = LAPACKE_zlange_work(LAPACK_COL_MAJOR, 'i', M, NRHS, Bref,    LDB, work);
    *Bplasmanorm = LAPACKE_zlange_work(LAPACK_COL_MAJOR, 'i', M, NRHS, Bplasma, LDB, work);

    cblas_ztrsm(CblasColMajor, (CBLAS_SIDE)side, (CBLAS_UPLO)uplo,
                (CBLAS_TRANSPOSE)trans, (CBLAS_DIAG)diag, M, NRHS,
                CBLAS_SADDR(alpha), A, LDA, Bref, LDB);

    *Blapacknorm = LAPACKE_zlange_work(LAPACK_COL_MAJOR, 'i', M, NRHS, Bref, LDB, work);

    cblas_zaxpy(LDB * NRHS, CBLAS_SADDR(beta_const), Bplasma, 1, Bref, 1);

    Rnorm = LAPACKE_zlange_work(LAPACK_COL_MAJOR, 'i', M, NRHS, Bref, LDB, work);
    Rnorm = Rnorm / (*Blapacknorm * max(M,NRHS) * eps);

    free(work);

    return Rnorm;
}

/*--------------------------------------------------------------
 * Check the solution
 */

double zcheck_solution(int M, int N, int NRHS, PLASMA_Complex64_t *A, int LDA,
                      PLASMA_Complex64_t *B,  PLASMA_Complex64_t *X, int LDB,
                      double *anorm, double *bnorm, double *xnorm )
{
/*     int info_solution; */
    double Rnorm = -1.00;
    PLASMA_Complex64_t zone  =  1.0;
    PLASMA_Complex64_t mzone = -1.0;
    double *work = (double *)malloc(max(M, N)* sizeof(double));

    *anorm = LAPACKE_zlange_work(LAPACK_COL_MAJOR, lapack_const(PlasmaInfNorm), M, N,    A, LDA, work);
    *xnorm = LAPACKE_zlange_work(LAPACK_COL_MAJOR, lapack_const(PlasmaInfNorm), M, NRHS, X, LDB, work);
    *bnorm = LAPACKE_zlange_work(LAPACK_COL_MAJOR, lapack_const(PlasmaInfNorm), N, NRHS, B, LDB, work);

    cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, M, NRHS, N, CBLAS_SADDR(zone), A, LDA, X, LDB, CBLAS_SADDR(mzone), B, LDB);

    Rnorm = LAPACKE_zlange_work(LAPACK_COL_MAJOR, lapack_const(PlasmaInfNorm), N, NRHS, B, M, work);

    free(work);

    return Rnorm;
}
