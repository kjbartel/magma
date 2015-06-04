/*
    -- MAGMA (version 1.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2011

       @generated c Sun Nov 13 20:48:38 2011

*/
#include "common_magma.h"

static __global__ void 
clacpy_generic(int M, int N, cuFloatComplex *A, int LDA, cuFloatComplex *B, int LDB ) 
{ 
    int ibx = blockIdx.x * 64;
    int tx  = threadIdx.x;
    int ty  = threadIdx.y;
    int idt = ty * 16 + tx;
    
    if( (ibx+idt) >=M ){
        A += M-1;
        B += M-1;
    }
    else{
        A += ibx+idt;
        B += ibx+idt;
    }
    const cuFloatComplex * Aend = A+LDA*N;
    cuFloatComplex Ap[1]={A[0]};
    do {
        A += LDA;
        B[0] = Ap[0];
        Ap[0]= A[0];        
        B += LDB;
    } while (A < Aend);
    B[0] = Ap[0];
}

static __global__ void 
clacpy_special(int M, int N, cuFloatComplex *A, int LDA, cuFloatComplex *B, int LDB ) 
{ 
    int ibx = blockIdx.x * 64;
    int tx  = threadIdx.x;
    int ty  = threadIdx.y;
    int idt = ty * 16 + tx;
    if( (ibx+idt) >=M ){
        A += M-1;
        B += M-1;
    }
    else{
        A += ibx+idt;
        B += ibx+idt;
    }
    cuFloatComplex Ap[1]={A[0]};
    B[0] = Ap[0];
}

extern "C" void 
magmablas_clacpy_64_64_16_4_v2(int M, int N, cuFloatComplex *A, int LDA, cuFloatComplex *B, int LDB )
{
        dim3 threads( 16, 4 );
        dim3 grid( M/64+(M%64!=0), 1 );
        if( N == 1 ) 
            clacpy_special<<< grid, threads, 0, magma_stream >>> ( M, N, A, LDA, B, LDB ) ;
        else        
            clacpy_generic<<< grid, threads, 0, magma_stream >>> ( M, N, A, LDA, B, LDB ) ;
}

extern "C" void 
magmablas_clacpy(char uplo, int M, int N, cuFloatComplex *A, int LDA, cuFloatComplex *B, int LDB)
{
/*
    Note
  ========
  - UPLO Parameter is disabled
  - Do we want to provide a generic function to the user with all the options?

  Purpose
  =======

  ZLACPY copies all or part of a two-dimensional matrix A to another
  matrix B.

  Arguments
  =========

  UPLO    (input) CHARACTER*1
          Specifies the part of the matrix A to be copied to B.
          = 'U':      Upper triangular part
          = 'L':      Lower triangular part
          Otherwise:  All of the matrix A

  M       (input) INTEGER
          The number of rows of the matrix A.  M >= 0.

  N       (input) INTEGER
          The number of columns of the matrix A.  N >= 0.

  A       (input) COMPLEX DOUBLE PRECISION array, dimension (LDA,N)
          The m by n matrix A.  If UPLO = 'U', only the upper triangle
          or trapezoid is accessed; if UPLO = 'L', only the lower
          triangle or trapezoid is accessed.

  LDA     (input) INTEGER
          The leading dimension of the array A.  LDA >= max(1,M).

  B       (output) COMPLEX DOUBLE PRECISION array, dimension (LDB,N)
          On exit, B = A in the locations specified by UPLO.

  LDB     (input) INTEGER
          The leading dimension of the array B.  LDB >= max(1,M).

  =====================================================================   */

    if ( (uplo == 'U') || (uplo == 'u') ) {
        fprintf(stderr, "Lacpy upper is not implemented\n");
    } else if ( (uplo == 'L') || (uplo == 'l') ) {
        fprintf(stderr, "Lacpy lower is not implemented\n");
    } else {
        magmablas_clacpy_64_64_16_4_v2(M, N, A, LDA, B, LDB);
    }
}
