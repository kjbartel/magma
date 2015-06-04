/**
 *
 * @precisions normal z -> c d s
 *
 **/
#ifndef ZAUXILIARY_H
#define ZAUXILIARY_H

int    zcheck_orthogonality   (int M, int N, int LDQ, PLASMA_Complex64_t *Q);
int    zcheck_QRfactorization (int M, int N, PLASMA_Complex64_t *A1, PLASMA_Complex64_t *A2, int LDA, PLASMA_Complex64_t *Q);
int    zcheck_LLTfactorization(int N, PLASMA_Complex64_t *A1, PLASMA_Complex64_t *A2, int LDA, int uplo);
double zcheck_gemm(PLASMA_enum transA, PLASMA_enum transB, int M, int N, int K,
                   PLASMA_Complex64_t alpha, PLASMA_Complex64_t *A, int LDA,
                   PLASMA_Complex64_t *B, int LDB,
                   PLASMA_Complex64_t beta, PLASMA_Complex64_t *Cplasma,
                   PLASMA_Complex64_t *Cref, int LDC,
                   double *Cinitnorm, double *Cplasmanorm, double *Clapacknorm );

double zcheck_trsm(PLASMA_enum side, PLASMA_enum uplo, PLASMA_enum trans, PLASMA_enum diag,
           int M, int NRHS, PLASMA_Complex64_t alpha,
           PLASMA_Complex64_t *A, int LDA,
                   PLASMA_Complex64_t *Bplasma, PLASMA_Complex64_t *Bref, int LDB,
                   double *Binitnorm, double *Bplasmanorm, double *Blapacknorm );

double zcheck_solution(int M, int N, int NRHS,
                      PLASMA_Complex64_t *A1, int LDA,
                      PLASMA_Complex64_t *B1, PLASMA_Complex64_t *B2, int LDB,
                      double *anorm, double *bnorm, double *xnorm);

#endif /* ZAUXILIARY_H */
