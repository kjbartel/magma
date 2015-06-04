/**
 *
 *  @file morse_z.h
 *
 *  MAGMA codelets kernel
 *  MAGMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver,
 *  and INRIA Bordeaux Sud-Ouest
 *
 *  @version 1.1.0
 *  @author Mathieu Faverge
 *  @author Cedric Augonnet
 *  @date 2011-06-01
 *  @precisions normal z -> c d s
 *
 **/
#ifndef _MORSE_Z_H_
#define _MORSE_Z_H_

#define PRECISION_z

void MORSE_zgemm( MorseOption_t *options, 
                  int transA, int transB,
                  int m, int n, int k, 
                  PLASMA_Complex64_t alpha, magma_desc_t *A, int Am, int An,
                                            magma_desc_t *B, int Bm, int Bn,
                  PLASMA_Complex64_t beta,  magma_desc_t *C, int Cm, int Cn);
void MORSE_zgeqrt( MorseOption_t *options, 
                   int m, int n, int ib, 
                   magma_desc_t *A, int Am, int An,
                   magma_desc_t *T, int Tm, int Tn);
void MORSE_zgessm( MorseOption_t *options, 
                   int m, int n, int k, int ib,
                   int *IPIV,
                   magma_desc_t *L, int Lm, int Ln,
                   magma_desc_t *D, int Dm, int Dn,
                   magma_desc_t *A, int Am, int An);
void MORSE_zgetrl( MorseOption_t *options, 
                   int m, int n, int ib,
                   magma_desc_t *A, int Am, int An,
                   magma_desc_t *L, int Lm, int Ln,
                   int *IPIV,
                   PLASMA_bool check, int iinfo);
void MORSE_zherk( MorseOption_t *options, 
                  int uplo, int trans,
                  int n, int k, 
                  double alpha, magma_desc_t *A, int Am, int An,
                  double beta,  magma_desc_t *C, int Cm, int Cn);
void MORSE_zlacpy( MorseOption_t *options, 
                   PLASMA_enum uplo, int m, int n,
                   magma_desc_t *A, int Am, int An,
                   magma_desc_t *B, int Bm, int Bn);
#if defined(PRECISION_z) || defined(PRECISION_c)
void MORSE_zplghe( MorseOption_t *options, 
                   double bump, int m, int n, magma_desc_t *A, int Am, int An,
                   int bigM, int m0, int n0, unsigned long long int seed );
#endif
void MORSE_zplgsy( MorseOption_t *options, 
                   PLASMA_Complex64_t bump, int m, int n, magma_desc_t *A, int Am, int An,
                   int bigM, int m0, int n0, unsigned long long int seed );
void MORSE_zplrnt( MorseOption_t *options, 
                   int m, int n, magma_desc_t *A, int Am, int An,
                   int bigM, int m0, int n0, unsigned long long int seed );
void MORSE_zpotrf( MorseOption_t *option, 
                   PLASMA_enum uplo, int n, 
                   magma_desc_t *A, int Am, int An, int iinfo);
void MORSE_zssssm( MorseOption_t *options, 
                   int m1, int n1, int m2, int n2, int k, int ib,
                   magma_desc_t *A1, int A1m, int A1n,
                   magma_desc_t *A2, int A2m, int A2n,
                   magma_desc_t *L1, int L1m, int L1n,
                   magma_desc_t *L2, int L2m, int L2n,
                   int *IPIV);
void MORSE_ztrsm( MorseOption_t *options, 
                  int side, int uplo, int transA, int diag,
                  int m, int n, 
                  PLASMA_Complex64_t alpha, 
                  magma_desc_t *A, int Am, int An,
                  magma_desc_t *B, int Bm, int Bn);
void MORSE_ztsmqr( MorseOption_t *options, 
                   int side, int trans,
                   int m1, int n1, int m2, int n2, int k, int ib,
                   magma_desc_t *A1, int A1m, int A1n,
                   magma_desc_t *A2, int A2m, int A2n,
                   magma_desc_t *V,  int Vm, int Vn,
                   magma_desc_t *T,  int Tm, int Tn);
void MORSE_ztsqrt( MorseOption_t *options, 
                   int m, int n, int ib,
                   magma_desc_t *A1, int A1m, int A1n,
                   magma_desc_t *A2, int A2m, int A2n,
                   magma_desc_t *T,  int Tm, int Tn);
void MORSE_ztstrf( MorseOption_t *options, 
                   int m, int n, int ib, int nb,
                   magma_desc_t *U, int Um, int Un,
                   magma_desc_t *A, int Am, int An,
                   magma_desc_t *L, int Lm, int Ln,
                   int  *IPIV,
                   PLASMA_bool check, int iinfo);
void MORSE_zunmqr( MorseOption_t *options, 
                   int side, int trans,
                   int m, int n, int k, int ib, 
                   magma_desc_t *A, int Am, int An,
                   magma_desc_t *T, int Tm, int Tn,
                   magma_desc_t *C, int Cm, int Cn);
#undef PRECISION_z

#endif /* _MORSE_Z_H_ */
