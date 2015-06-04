/**
 *
 * @file magma_morse_z.h
 *
 *  PLASMA header file for double _Complex routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.4.5
 * @author Jakub Kurzak
 * @author Hatem Ltaief
 * @author Mathieu Faverge
 * @author Azzam Haidar
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 **/
#ifndef _MAGMA_MORSE_Z_H_
#define _MAGMA_MORSE_Z_H_

#ifdef __cplusplus
extern "C" {
#endif

#define PRECISION_z

/** ****************************************************************************
 *  Declarations of math functions (LAPACK layout) - alphabetical order
 **/
int MAGMA_zgebrd(int M, int N, PLASMA_Complex64_t *A, int LDA, double *D, double *E, magma_desc_t *T);
int MAGMA_zgeev(PLASMA_enum jobvl, PLASMA_enum jobvr, int N, PLASMA_Complex64_t *A, int LDA, PLASMA_Complex64_t *W, PLASMA_Complex64_t *VL, int LDVL, PLASMA_Complex64_t *VR, int LDVR, magma_desc_t *T);
int MAGMA_zgehrd(int N, int ILO, int IHI, PLASMA_Complex64_t *A, int LDA, magma_desc_t *T);
int MAGMA_zgelqf(int M, int N, PLASMA_Complex64_t *A, int LDA, magma_desc_t *T);
int MAGMA_zgelqs(int M, int N, int NRHS, PLASMA_Complex64_t *A, int LDA, magma_desc_t *T, PLASMA_Complex64_t *B, int LDB);
int MAGMA_zgels(PLASMA_enum trans, int M, int N, int NRHS, PLASMA_Complex64_t *A, int LDA, magma_desc_t *T, PLASMA_Complex64_t *B, int LDB);
int MAGMA_zgemm(PLASMA_enum transA, PLASMA_enum transB, int M, int N, int K, PLASMA_Complex64_t alpha, PLASMA_Complex64_t *A, int LDA, PLASMA_Complex64_t *B, int LDB, PLASMA_Complex64_t beta, PLASMA_Complex64_t *C, int LDC);
int MAGMA_zgeqrf(int M, int N, PLASMA_Complex64_t *A, int LDA, magma_desc_t *T);
int MAGMA_zgeqrs(int M, int N, int NRHS, PLASMA_Complex64_t *A, int LDA, magma_desc_t *T, PLASMA_Complex64_t *B, int LDB);
int MAGMA_zgesv(int N, int NRHS, PLASMA_Complex64_t *A, int LDA, int *IPIV, PLASMA_Complex64_t *B, int LDB);
int MAGMA_zgesv_incpiv(int N, int NRHS, PLASMA_Complex64_t *A, int LDA, magma_desc_t *L, int *IPIV, PLASMA_Complex64_t *B, int LDB);
int MAGMA_zgesvd(PLASMA_enum jobu, PLASMA_enum jobvt, int M, int N, PLASMA_Complex64_t *A, int LDA, double *S, PLASMA_Complex64_t *U, int LDU, PLASMA_Complex64_t *VT, int LDVT, magma_desc_t *T);
int MAGMA_zgetrf(int M, int N, PLASMA_Complex64_t *A, int LDA, int *IPIV);
int MAGMA_zgetrf_incpiv(int M, int N, PLASMA_Complex64_t *A, int LDA, magma_desc_t *L, int *IPIV);
int MAGMA_zgetri(int N, PLASMA_Complex64_t *A, int LDA, int *IPIV);
int MAGMA_zgetrs(PLASMA_enum trans, int N, int NRHS, PLASMA_Complex64_t *A, int LDA, int *IPIV, PLASMA_Complex64_t *B, int LDB);
int MAGMA_zgetrs_incpiv(PLASMA_enum trans, int N, int NRHS, PLASMA_Complex64_t *A, int LDA, magma_desc_t *L, int *IPIV, PLASMA_Complex64_t *B, int LDB);
#if defined(PRECISION_z) || defined(PRECISION_c)
int MAGMA_zhemm(PLASMA_enum side, PLASMA_enum uplo, int M, int N, PLASMA_Complex64_t alpha, PLASMA_Complex64_t *A, int LDA, PLASMA_Complex64_t *B, int LDB, PLASMA_Complex64_t beta, PLASMA_Complex64_t *C, int LDC);
int MAGMA_zherk(PLASMA_enum uplo, PLASMA_enum trans, int N, int K, double alpha, PLASMA_Complex64_t *A, int LDA, double beta, PLASMA_Complex64_t *C, int LDC);
int MAGMA_zher2k(PLASMA_enum uplo, PLASMA_enum trans, int N, int K, PLASMA_Complex64_t alpha, PLASMA_Complex64_t *A, int LDA, PLASMA_Complex64_t *B, int LDB, double beta, PLASMA_Complex64_t *C, int LDC);
#endif
int MAGMA_zheev(PLASMA_enum jobz, PLASMA_enum uplo, int N, PLASMA_Complex64_t *A, int LDA, double *W, magma_desc_t *T, PLASMA_Complex64_t *Q, int LDQ);
int MAGMA_zhegv(PLASMA_enum itype, PLASMA_enum jobz, PLASMA_enum uplo, int N, PLASMA_Complex64_t *A, int LDA, PLASMA_Complex64_t *B, int LDB, double *W, magma_desc_t *T, PLASMA_Complex64_t *Q, int LDQ);
int MAGMA_zhegst(PLASMA_enum itype, PLASMA_enum uplo, int N, PLASMA_Complex64_t *A, int LDA, PLASMA_Complex64_t *B, int LDB);
int MAGMA_zhetrd(PLASMA_enum uplo, int N, PLASMA_Complex64_t *A, int LDA, double *D, double *E, magma_desc_t *T);
int MAGMA_zlacpy(PLASMA_enum uplo, int M, int N, PLASMA_Complex64_t *A, int LDA, PLASMA_Complex64_t *B, int LDB);
double MAGMA_zlange(PLASMA_enum norm, int M, int N, PLASMA_Complex64_t *A, int LDA, double *work);
#if defined(PRECISION_z) || defined(PRECISION_c)
double MAGMA_zlanhe(PLASMA_enum norm, PLASMA_enum uplo, int N, PLASMA_Complex64_t *A, int LDA, double *work);
#endif
double MAGMA_zlansy(PLASMA_enum norm, PLASMA_enum uplo, int N, PLASMA_Complex64_t *A, int LDA, double *work);
int MAGMA_zlaset(PLASMA_enum uplo, int M, int N, PLASMA_Complex64_t alpha, PLASMA_Complex64_t beta, PLASMA_Complex64_t *A, int LDA);
int MAGMA_zlaswp(int N, PLASMA_Complex64_t *A, int LDA, int K1, int K2, int *IPIV, int INCX);
int MAGMA_zlaswpc(int N, PLASMA_Complex64_t *A, int LDA, int K1, int K2, int *IPIV, int INCX);
int MAGMA_zlauum(PLASMA_enum uplo, int N, PLASMA_Complex64_t *A, int LDA);
#if defined(PRECISION_z) || defined(PRECISION_c)
int MAGMA_zplghe( double bump, int N, PLASMA_Complex64_t *A, int LDA, unsigned long long int seed );
#endif
int MAGMA_zplgsy( PLASMA_Complex64_t bump, int N, PLASMA_Complex64_t *A, int LDA, unsigned long long int seed );
int MAGMA_zplrnt( int M, int N, PLASMA_Complex64_t *A, int LDA, unsigned long long int seed );
int MAGMA_zposv(PLASMA_enum uplo, int N, int NRHS, PLASMA_Complex64_t *A, int LDA, PLASMA_Complex64_t *B, int LDB);
int MAGMA_zpotrf(PLASMA_enum uplo, int N, PLASMA_Complex64_t *A, int LDA);
int MAGMA_zpotri(PLASMA_enum uplo, int N, PLASMA_Complex64_t *A, int LDA);
int MAGMA_zpotrs(PLASMA_enum uplo, int N, int NRHS, PLASMA_Complex64_t *A, int LDA, PLASMA_Complex64_t *B, int LDB);
int MAGMA_zsymm(PLASMA_enum side, PLASMA_enum uplo, int M, int N, PLASMA_Complex64_t alpha, PLASMA_Complex64_t *A, int LDA, PLASMA_Complex64_t *B, int LDB, PLASMA_Complex64_t beta, PLASMA_Complex64_t *C, int LDC);
int MAGMA_zsyrk(PLASMA_enum uplo, PLASMA_enum trans, int N, int K, PLASMA_Complex64_t alpha, PLASMA_Complex64_t *A, int LDA, PLASMA_Complex64_t beta, PLASMA_Complex64_t *C, int LDC);
int MAGMA_zsyr2k(PLASMA_enum uplo, PLASMA_enum trans, int N, int K, PLASMA_Complex64_t alpha, PLASMA_Complex64_t *A, int LDA, PLASMA_Complex64_t *B, int LDB, PLASMA_Complex64_t beta, PLASMA_Complex64_t *C, int LDC);
int MAGMA_ztrmm(PLASMA_enum side, PLASMA_enum uplo, PLASMA_enum transA, PLASMA_enum diag, int N, int NRHS, PLASMA_Complex64_t alpha, PLASMA_Complex64_t *A, int LDA, PLASMA_Complex64_t *B, int LDB);
int MAGMA_ztrsm(PLASMA_enum side, PLASMA_enum uplo, PLASMA_enum transA, PLASMA_enum diag, int N, int NRHS, PLASMA_Complex64_t alpha, PLASMA_Complex64_t *A, int LDA, PLASMA_Complex64_t *B, int LDB);
int MAGMA_ztrsmpl(int N, int NRHS, PLASMA_Complex64_t *A, int LDA, magma_desc_t *L, int *IPIV, PLASMA_Complex64_t *B, int LDB);
int MAGMA_ztrsmrv(PLASMA_enum side, PLASMA_enum uplo, PLASMA_enum transA, PLASMA_enum diag, int N, int NRHS, PLASMA_Complex64_t alpha, PLASMA_Complex64_t *A, int LDA, PLASMA_Complex64_t *B, int LDB);
int MAGMA_ztrtri(PLASMA_enum uplo, PLASMA_enum diag, int N, PLASMA_Complex64_t *A, int LDA);
int MAGMA_zungbr(PLASMA_enum side, int M, int N, int K, PLASMA_Complex64_t *A, int LDA, magma_desc_t *T, PLASMA_Complex64_t *Q, int LDQ);
int MAGMA_zunghr(int N, int ILO, int IHI, PLASMA_Complex64_t *A, int LDA, magma_desc_t *T, PLASMA_Complex64_t *Q, int LDQ);
int MAGMA_zunglq(int M, int N, int K, PLASMA_Complex64_t *A, int LDA, magma_desc_t *T, PLASMA_Complex64_t *B, int LDB);
int MAGMA_zungqr(int M, int N, int K, PLASMA_Complex64_t *A, int LDA, magma_desc_t *T, PLASMA_Complex64_t *B, int LDB);
int MAGMA_zungtr(PLASMA_enum uplo, int N, PLASMA_Complex64_t *A, int LDA, magma_desc_t *T, PLASMA_Complex64_t *B, int LDB);
int MAGMA_zunmlq(PLASMA_enum side, PLASMA_enum trans, int M, int N, int K, PLASMA_Complex64_t *A, int LDA, magma_desc_t *T, PLASMA_Complex64_t *B, int LDB);
int MAGMA_zunmqr(PLASMA_enum side, PLASMA_enum trans, int M, int N, int K, PLASMA_Complex64_t *A, int LDA, magma_desc_t *T, PLASMA_Complex64_t *B, int LDB);

int MAGMA_zgecfi(int m, int n, PLASMA_Complex64_t *A, PLASMA_enum fin, int imb, int inb, PLASMA_enum fout, int omb, int onb);
int MAGMA_zgetmi(int m, int n, PLASMA_Complex64_t *A, PLASMA_enum fin, int mb,  int nb);

/** ****************************************************************************
 *  Declarations of math functions (tile layout) - alphabetical order
 **/
int MAGMA_zgebrd_Tile(magma_desc_t *A, double *D, double *E, magma_desc_t *T);
int MAGMA_zgeev_Tile(PLASMA_enum jobvl, PLASMA_enum jobvr, magma_desc_t *A, PLASMA_Complex64_t *W, magma_desc_t *VL, magma_desc_t *VR, magma_desc_t *T);
int MAGMA_zgehrd_Tile(magma_desc_t *A, magma_desc_t *T);
int MAGMA_zgelqf_Tile(magma_desc_t *A, magma_desc_t *T);
int MAGMA_zgelqs_Tile(magma_desc_t *A, magma_desc_t *T, magma_desc_t *B);
int MAGMA_zgels_Tile(PLASMA_enum trans, magma_desc_t *A, magma_desc_t *T, magma_desc_t *B);
int MAGMA_zgemm_Tile(PLASMA_enum transA, PLASMA_enum transB, PLASMA_Complex64_t alpha, magma_desc_t *A, magma_desc_t *B, PLASMA_Complex64_t beta, magma_desc_t *C);
int MAGMA_zgeqrf_Tile(magma_desc_t *A, magma_desc_t *T);
int MAGMA_zgeqrs_Tile(magma_desc_t *A, magma_desc_t *T, magma_desc_t *B);
int MAGMA_zgesv_Tile(magma_desc_t *A, int *IPIV, magma_desc_t *B);
int MAGMA_zgesv_incpiv_Tile(magma_desc_t *A, magma_desc_t *L, int *IPIV, magma_desc_t *B);
int MAGMA_zgesvd_Tile(PLASMA_enum jobu, PLASMA_enum jobvt, magma_desc_t *A, double *S, magma_desc_t *U, magma_desc_t *VT, magma_desc_t *T);
int MAGMA_zgetrf_Tile(magma_desc_t *A, int *IPIV);
int MAGMA_zgetrf_incpiv_Tile(magma_desc_t *A, magma_desc_t *L, int *IPIV);
int MAGMA_zgetri_Tile(magma_desc_t *A, int *IPIV);
int MAGMA_zgetrs_Tile(PLASMA_enum trans, magma_desc_t *A, int *IPIV, magma_desc_t *B);
int MAGMA_zgetrs_incpiv_Tile(magma_desc_t *A, magma_desc_t *L, int *IPIV, magma_desc_t *B);
#if defined(PRECISION_z) || defined(PRECISION_c)
int MAGMA_zhemm_Tile(PLASMA_enum side, PLASMA_enum uplo, PLASMA_Complex64_t alpha, magma_desc_t *A, magma_desc_t *B, PLASMA_Complex64_t beta, magma_desc_t *C);
int MAGMA_zherk_Tile(PLASMA_enum uplo, PLASMA_enum trans, double alpha, magma_desc_t *A, double beta, magma_desc_t *C);
int MAGMA_zher2k_Tile(PLASMA_enum uplo, PLASMA_enum trans, PLASMA_Complex64_t alpha, magma_desc_t *A, magma_desc_t *B, double beta, magma_desc_t *C);
#endif
int MAGMA_zheev_Tile(PLASMA_enum jobz, PLASMA_enum uplo, magma_desc_t *A, double *W, magma_desc_t *T, magma_desc_t *Q);
int MAGMA_zhegv_Tile(PLASMA_enum itype, PLASMA_enum jobz, PLASMA_enum uplo, magma_desc_t *A, magma_desc_t *B, double *W, magma_desc_t *T, magma_desc_t *Q);
int MAGMA_zhegst_Tile(PLASMA_enum itype, PLASMA_enum uplo, magma_desc_t *A, magma_desc_t *B);
int MAGMA_zhetrd_Tile(PLASMA_enum uplo, magma_desc_t *A, double *D, double *E, magma_desc_t *T);
int MAGMA_zlacpy_Tile(PLASMA_enum uplo, magma_desc_t *A, magma_desc_t *B);
double MAGMA_zlange_Tile(PLASMA_enum norm, magma_desc_t *A, double *work);
#if defined(PRECISION_z) || defined(PRECISION_c)
double MAGMA_zlanhe_Tile(PLASMA_enum norm, PLASMA_enum uplo, magma_desc_t *A, double *work);
#endif
double MAGMA_zlansy_Tile(PLASMA_enum norm, PLASMA_enum uplo, magma_desc_t *A, double *work);
int MAGMA_zlaset_Tile(PLASMA_enum uplo, PLASMA_Complex64_t alpha, PLASMA_Complex64_t beta, magma_desc_t *A);
int MAGMA_zlaswp_Tile(magma_desc_t *A, int K1, int K2, int *IPIV, int INCX);
int MAGMA_zlaswpc_Tile(magma_desc_t *A, int K1, int K2, int *IPIV, int INCX);
int MAGMA_zlauum_Tile(PLASMA_enum uplo, magma_desc_t *A);
#if defined(PRECISION_z) || defined(PRECISION_c)
int MAGMA_zplghe_Tile(double bump, magma_desc_t *A, unsigned long long int seed );
#endif
int MAGMA_zplgsy_Tile(PLASMA_Complex64_t bump, magma_desc_t *A, unsigned long long int seed );
int MAGMA_zplrnt_Tile(magma_desc_t *A, unsigned long long int seed );
int MAGMA_zposv_Tile(PLASMA_enum uplo, magma_desc_t *A, magma_desc_t *B);
int MAGMA_zpotrf_Tile(PLASMA_enum uplo, magma_desc_t *A);
int MAGMA_zpotri_Tile(PLASMA_enum uplo, magma_desc_t *A);
int MAGMA_zpotrs_Tile(PLASMA_enum uplo, magma_desc_t *A, magma_desc_t *B);
int MAGMA_zsymm_Tile(PLASMA_enum side, PLASMA_enum uplo, PLASMA_Complex64_t alpha, magma_desc_t *A, magma_desc_t *B, PLASMA_Complex64_t beta, magma_desc_t *C);
int MAGMA_zsyrk_Tile(PLASMA_enum uplo, PLASMA_enum trans, PLASMA_Complex64_t alpha, magma_desc_t *A, PLASMA_Complex64_t beta, magma_desc_t *C);
int MAGMA_zsyr2k_Tile(PLASMA_enum uplo, PLASMA_enum trans, PLASMA_Complex64_t alpha, magma_desc_t *A, magma_desc_t *B, PLASMA_Complex64_t beta, magma_desc_t *C);
int MAGMA_ztrmm_Tile(PLASMA_enum side, PLASMA_enum uplo, PLASMA_enum transA, PLASMA_enum diag, PLASMA_Complex64_t alpha, magma_desc_t *A, magma_desc_t *B);
int MAGMA_ztrsm_Tile(PLASMA_enum side, PLASMA_enum uplo, PLASMA_enum transA, PLASMA_enum diag, PLASMA_Complex64_t alpha, magma_desc_t *A, magma_desc_t *B);
int MAGMA_ztrsmpl_Tile(magma_desc_t *A, magma_desc_t *L, int *IPIV, magma_desc_t *B);
int MAGMA_ztrsmrv_Tile(PLASMA_enum side, PLASMA_enum uplo, PLASMA_enum transA, PLASMA_enum diag, PLASMA_Complex64_t alpha, magma_desc_t *A, magma_desc_t *B);
int MAGMA_ztrtri_Tile(PLASMA_enum uplo, PLASMA_enum diag, magma_desc_t *A);
int MAGMA_zungbr_Tile(PLASMA_enum size, magma_desc_t *A, magma_desc_t *T, magma_desc_t *Q);
int MAGMA_zunghr_Tile(magma_desc_t *A, magma_desc_t *T, magma_desc_t *Q);
int MAGMA_zunglq_Tile(magma_desc_t *A, magma_desc_t *T, magma_desc_t *B);
int MAGMA_zungqr_Tile(magma_desc_t *A, magma_desc_t *T, magma_desc_t *B);
int MAGMA_zungtr_Tile(PLASMA_enum uplo, magma_desc_t *A, magma_desc_t *T, magma_desc_t *B);
int MAGMA_zunmlq_Tile(PLASMA_enum side, PLASMA_enum trans, magma_desc_t *A, magma_desc_t *T, magma_desc_t *B);
int MAGMA_zunmqr_Tile(PLASMA_enum side, PLASMA_enum trans, magma_desc_t *A, magma_desc_t *T, magma_desc_t *B);

/** ****************************************************************************
 *  Declarations of math functions (tile layout, asynchronous execution) - alphabetical order
 **/
int MAGMA_zgebrd_Tile_Async(magma_desc_t *A, double *D, double *E, magma_desc_t *T, magma_sequence_t *sequence, magma_request_t *request);
int MAGMA_zgeev_Tile_Async(PLASMA_enum jobvl, PLASMA_enum jobvr, magma_desc_t *A, PLASMA_Complex64_t *W, magma_desc_t *VL, magma_desc_t *VR, magma_desc_t *T, magma_sequence_t *sequence, magma_request_t *request);
int MAGMA_zgehrd_Tile_Async(magma_desc_t *A, magma_desc_t *T, magma_sequence_t *sequence, magma_request_t *request);
int MAGMA_zgelqf_Tile_Async(magma_desc_t *A, magma_desc_t *T, magma_sequence_t *sequence, magma_request_t *request);
int MAGMA_zgelqs_Tile_Async(magma_desc_t *A, magma_desc_t *T, magma_desc_t *B, magma_sequence_t *sequence, magma_request_t *request);
int MAGMA_zgels_Tile_Async(PLASMA_enum trans, magma_desc_t *A, magma_desc_t *T, magma_desc_t *B, magma_sequence_t *sequence, magma_request_t *request);
int MAGMA_zgemm_Tile_Async(PLASMA_enum transA, PLASMA_enum transB, PLASMA_Complex64_t alpha, magma_desc_t *A, magma_desc_t *B, PLASMA_Complex64_t beta, magma_desc_t *C, magma_sequence_t *sequence, magma_request_t *request);
int MAGMA_zgeqrf_Tile_Async(magma_desc_t *A, magma_desc_t *T, magma_sequence_t *sequence, magma_request_t *request);
int MAGMA_zgeqrs_Tile_Async(magma_desc_t *A, magma_desc_t *T, magma_desc_t *B, magma_sequence_t *sequence, magma_request_t *request);
int MAGMA_zgesv_Tile_Async(magma_desc_t *A, int *IPIV, magma_desc_t *B, magma_sequence_t *sequence, magma_request_t *request);
int MAGMA_zgesv_incpiv_Tile_Async(magma_desc_t *A, magma_desc_t *L, int *IPIV, magma_desc_t *B, magma_sequence_t *sequence, magma_request_t *request);
int MAGMA_zgesvd_Tile_Async(PLASMA_enum jobu, PLASMA_enum jobvt, magma_desc_t *A, double *S, magma_desc_t *U, magma_desc_t *VT, magma_desc_t *T, magma_sequence_t *sequence, magma_request_t *request);
int MAGMA_zgetrf_Tile_Async(magma_desc_t *A, int *IPIV, magma_sequence_t *sequence, magma_request_t *request);
int MAGMA_zgetrf_incpiv_Tile_Async(magma_desc_t *A, magma_desc_t *L, int *IPIV, magma_sequence_t *sequence, magma_request_t *request);
int MAGMA_zgetri_Tile_Async(magma_desc_t *A, int *IPIV, magma_desc_t *W, magma_sequence_t *sequence, magma_request_t *request);
int MAGMA_zgetrs_Tile_Async(PLASMA_enum trans, magma_desc_t *A, int *IPIV, magma_desc_t *B, magma_sequence_t *sequence, magma_request_t *request);
int MAGMA_zgetrs_incpiv_Tile_Async(magma_desc_t *A, magma_desc_t *L, int *IPIV, magma_desc_t *B, magma_sequence_t *sequence, magma_request_t *request);
#if defined(PRECISION_z) || defined(PRECISION_c)
int MAGMA_zhemm_Tile_Async(PLASMA_enum side, PLASMA_enum uplo, PLASMA_Complex64_t alpha, magma_desc_t *A, magma_desc_t *B, PLASMA_Complex64_t beta, magma_desc_t *C, magma_sequence_t *sequence, magma_request_t *request);
int MAGMA_zherk_Tile_Async(PLASMA_enum uplo, PLASMA_enum trans, double alpha, magma_desc_t *A, double beta, magma_desc_t *C, magma_sequence_t *sequence, magma_request_t *request);
int MAGMA_zher2k_Tile_Async(PLASMA_enum uplo, PLASMA_enum trans, PLASMA_Complex64_t alpha, magma_desc_t *A, magma_desc_t *B, double beta, magma_desc_t *C, magma_sequence_t *sequence, magma_request_t *request);
#endif
int MAGMA_zheev_Tile_Async(PLASMA_enum jobz, PLASMA_enum uplo, magma_desc_t *A, double *W, magma_desc_t *T, magma_desc_t *Q, magma_sequence_t *sequence, magma_request_t *request);
int MAGMA_zhegv_Tile_Async(PLASMA_enum itype, PLASMA_enum jobz, PLASMA_enum uplo, magma_desc_t *A, magma_desc_t *B, double *W, magma_desc_t *T, magma_desc_t *Q, magma_sequence_t *sequence, magma_request_t *request);
int MAGMA_zhegst_Tile_Async(PLASMA_enum itype, PLASMA_enum uplo, magma_desc_t *A, magma_desc_t *B, magma_sequence_t *sequence, magma_request_t *request);
int MAGMA_zhetrd_Tile_Async(PLASMA_enum uplo, magma_desc_t *A, double *D, double *E, magma_desc_t *T, magma_sequence_t *sequence, magma_request_t *request);
int MAGMA_zlacpy_Tile_Async(PLASMA_enum uplo, magma_desc_t *A, magma_desc_t *B, magma_sequence_t *sequence, magma_request_t *request);
int MAGMA_zlange_Tile_Async(PLASMA_enum norm, magma_desc_t *A, double *work, double *value, magma_sequence_t *sequence, magma_request_t *request);
#if defined(PRECISION_z) || defined(PRECISION_c)
int MAGMA_zlanhe_Tile_Async(PLASMA_enum norm, PLASMA_enum uplo, magma_desc_t *A, double *work, double *value, magma_sequence_t *sequence, magma_request_t *request);
#endif
int MAGMA_zlansy_Tile_Async(PLASMA_enum norm, PLASMA_enum uplo, magma_desc_t *A, double *work, double *value, magma_sequence_t *sequence, magma_request_t *request);
int MAGMA_zlaset_Tile_Async(PLASMA_enum uplo, PLASMA_Complex64_t alpha, PLASMA_Complex64_t beta, magma_desc_t *A, magma_sequence_t *sequence, magma_request_t *request);
int MAGMA_zlaswp_Tile_Async(magma_desc_t *A, int K1, int K2, int *IPIV, int INCX, magma_sequence_t *sequence, magma_request_t *request);
int MAGMA_zlaswpc_Tile_Async(magma_desc_t *A, int K1, int K2, int *IPIV, int INCX, magma_sequence_t *sequence, magma_request_t *request);
int MAGMA_zlauum_Tile_Async(PLASMA_enum uplo, magma_desc_t *A, magma_sequence_t *sequence, magma_request_t *request);
#if defined(PRECISION_z) || defined(PRECISION_c)
int MAGMA_zplghe_Tile_Async(double bump, magma_desc_t *A, unsigned long long int seed, magma_sequence_t *sequence, magma_request_t *request );
#endif
int MAGMA_zplgsy_Tile_Async(PLASMA_Complex64_t bump, magma_desc_t *A, unsigned long long int seed, magma_sequence_t *sequence, magma_request_t *request );
int MAGMA_zplrnt_Tile_Async(magma_desc_t *A, unsigned long long int seed, magma_sequence_t *sequence, magma_request_t *request );
int MAGMA_zposv_Tile_Async(PLASMA_enum uplo, magma_desc_t *A, magma_desc_t *B, magma_sequence_t *sequence, magma_request_t *request);
int MAGMA_zpotrf_Tile_Async(PLASMA_enum uplo, magma_desc_t *A, magma_sequence_t *sequence, magma_request_t *request);
int MAGMA_zpotri_Tile_Async(PLASMA_enum uplo, magma_desc_t *A, magma_sequence_t *sequence, magma_request_t *request);
int MAGMA_zpotrs_Tile_Async(PLASMA_enum uplo, magma_desc_t *A, magma_desc_t *B, magma_sequence_t *sequence, magma_request_t *request);
int MAGMA_zsymm_Tile_Async(PLASMA_enum side, PLASMA_enum uplo, PLASMA_Complex64_t alpha, magma_desc_t *A, magma_desc_t *B, PLASMA_Complex64_t beta, magma_desc_t *C, magma_sequence_t *sequence, magma_request_t *request);
int MAGMA_zsyrk_Tile_Async(PLASMA_enum uplo, PLASMA_enum trans, PLASMA_Complex64_t alpha, magma_desc_t *A, PLASMA_Complex64_t beta, magma_desc_t *C, magma_sequence_t *sequence, magma_request_t *request);
int MAGMA_zsyr2k_Tile_Async(PLASMA_enum uplo, PLASMA_enum trans, PLASMA_Complex64_t alpha, magma_desc_t *A, magma_desc_t *B, PLASMA_Complex64_t beta, magma_desc_t *C, magma_sequence_t *sequence, magma_request_t *request);
int MAGMA_ztrmm_Tile_Async(PLASMA_enum side, PLASMA_enum uplo, PLASMA_enum transA, PLASMA_enum diag, PLASMA_Complex64_t alpha, magma_desc_t *A, magma_desc_t *B, magma_sequence_t *sequence, magma_request_t *request);
int MAGMA_ztrsm_Tile_Async(PLASMA_enum side, PLASMA_enum uplo, PLASMA_enum transA, PLASMA_enum diag, PLASMA_Complex64_t alpha, magma_desc_t *A, magma_desc_t *B, magma_sequence_t *sequence, magma_request_t *request);
int MAGMA_ztrsmpl_Tile_Async(magma_desc_t *A, magma_desc_t *L, int *IPIV, magma_desc_t *B, magma_sequence_t *sequence, magma_request_t *request);
int MAGMA_ztrsmrv_Tile_Async(PLASMA_enum side, PLASMA_enum uplo, PLASMA_enum transA, PLASMA_enum diag, PLASMA_Complex64_t alpha, magma_desc_t *A, magma_desc_t *B, magma_sequence_t *sequence, magma_request_t *request);
int MAGMA_ztrtri_Tile_Async(PLASMA_enum uplo, PLASMA_enum diag, magma_desc_t *A, magma_sequence_t *sequence, magma_request_t *request);
int MAGMA_zungbr_Tile_Async(PLASMA_enum side, magma_desc_t *A, magma_desc_t *T, magma_desc_t *Q, magma_sequence_t *sequence, magma_request_t *request);
int MAGMA_zunghr_Tile_Async(magma_desc_t *A, magma_desc_t *T, magma_desc_t *Q, magma_sequence_t *sequence, magma_request_t *request);
int MAGMA_zunglq_Tile_Async(magma_desc_t *A, magma_desc_t *T, magma_desc_t *B, magma_sequence_t *sequence, magma_request_t *request);
int MAGMA_zungqr_Tile_Async(magma_desc_t *A, magma_desc_t *T, magma_desc_t *B, magma_sequence_t *sequence, magma_request_t *request);
int MAGMA_zungtr_Tile_Async(PLASMA_enum uplo, magma_desc_t *A, magma_desc_t *T, magma_desc_t *B, magma_sequence_t *sequence, magma_request_t *request);
int MAGMA_zunmlq_Tile_Async(PLASMA_enum side, PLASMA_enum trans, magma_desc_t *A, magma_desc_t *T, magma_desc_t *B, magma_sequence_t *sequence, magma_request_t *request);
int MAGMA_zunmqr_Tile_Async(PLASMA_enum side, PLASMA_enum trans, magma_desc_t *A, magma_desc_t *T, magma_desc_t *B, magma_sequence_t *sequence, magma_request_t *request);

int MAGMA_zgecfi_Async(int m, int n, PLASMA_Complex64_t *A, PLASMA_enum f_in, int imb, int inb, PLASMA_enum f_out, int omb, int onb, magma_sequence_t *sequence, magma_request_t *request);
int MAGMA_zgetmi_Async(int m, int n, PLASMA_Complex64_t *A, PLASMA_enum f_in, int mb,  int inb, magma_sequence_t *sequence, magma_request_t *request);

/** ****************************************************************************
 *  Declarations of workspace allocation functions (tile layout) - alphabetical order
 **/
int MAGMA_Alloc_Workspace_zgebrd(       int M, int N, magma_desc_t **T);
int MAGMA_Alloc_Workspace_zgeev(        int N,        magma_desc_t **T);
int MAGMA_Alloc_Workspace_zgehrd(       int N,        magma_desc_t **T);
int MAGMA_Alloc_Workspace_zgelqf(       int M, int N, magma_desc_t **T);
int MAGMA_Alloc_Workspace_zgels(        int M, int N, magma_desc_t **T);
int MAGMA_Alloc_Workspace_zgeqrf(       int M, int N, magma_desc_t **T);
int MAGMA_Alloc_Workspace_zgesv_incpiv( int N,        magma_desc_t **L, int **IPIV);
int MAGMA_Alloc_Workspace_zgesvd(       int M, int N, magma_desc_t **T);
int MAGMA_Alloc_Workspace_zgetrf_incpiv(int M, int N, magma_desc_t **L, int **IPIV);
int MAGMA_Alloc_Workspace_zheev(        int M, int N, magma_desc_t **T);
int MAGMA_Alloc_Workspace_zhegv(        int M, int N, magma_desc_t **T);
int MAGMA_Alloc_Workspace_zhetrd(       int M, int N, magma_desc_t **T);

/** ****************************************************************************
 *  Auxiliary function prototypes
 **/
int MAGMA_zLapack_to_Tile(PLASMA_Complex64_t *Af77, int LDA, magma_desc_t *A);
int MAGMA_zTile_to_Lapack(magma_desc_t *A, PLASMA_Complex64_t *Af77, int LDA);
int MAGMA_zLapack_to_Tile_Async(PLASMA_Complex64_t *Af77, int LDA, magma_desc_t *A, magma_sequence_t *sequence, magma_request_t *request);
int MAGMA_zTile_to_Lapack_Async(magma_desc_t *A, PLASMA_Complex64_t *Af77, int LDA, magma_sequence_t *sequence, magma_request_t *request);

#undef PRECISION_z

#ifdef __cplusplus
}
#endif

#endif /* _MAGMA_MORSE_Z_H_ */
