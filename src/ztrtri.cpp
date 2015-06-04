/*
    -- MAGMA (version 1.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2011

       @precisions normal z -> s d c

*/
#include "common_magma.h"

// === Define what BLAS to use ============================================
 #define PRECISION_z
 #if (defined(PRECISION_s) || defined(PRECISION_d))
   #define cublasZgemm magmablas_zgemm
     #define cublasZtrsm magmablas_ztrsm
     #endif

     #if (GPUSHMEM >= 200)
     #if (defined(PRECISION_s))
          #undef  cublasSgemm
               #define cublasSgemm magmablas_sgemm_fermi80
                 #endif
                 #endif
// === End defining what BLAS to use ======================================

                 #define A(i, j)  (a   +(j)*lda  + (i))
                 #define dA(i, j) (work+(j)*ldda + (i))

                 
extern "C" magma_int_t
magma_ztrtri(char uplo, char diag, magma_int_t n,
              cuDoubleComplex *a, magma_int_t lda, magma_int_t *info)
{
/*  -- MAGMA (version 1.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2011

    Purpose
    =======

        DTRTRI computes the inverse of a real upper or lower triangular
        matrix A.

        This is the Level 3 BLAS version of the algorithm.

        Arguments
        =========

        UPLO    (input) CHARACTER*1
                        = 'U':  A is upper triangular;
                        = 'L':  A is lower triangular.

        DIAG    (input) CHARACTER*1
                        = 'N':  A is non-unit triangular;
                        = 'U':  A is unit triangular.

        N       (input) INTEGER
                        The order of the matrix A.  N >= 0.

        A       (input/output) COMPLEX_16 array, dimension (LDA,N)
                        On entry, the triangular matrix A.  If UPLO = 'U', the
                        leading N-by-N upper triangular part of the array A contains
                        the upper triangular matrix, and the strictly lower
                        triangular part of A is not referenced.  If UPLO = 'L', the
                        leading N-by-N lower triangular part of the array A contains
                        the lower triangular matrix, and the strictly upper
                        triangular part of A is not referenced.  If DIAG = 'U', the
                        diagonal elements of A are also not referenced and are
                        assumed to be 1.
                        On exit, the (triangular) inverse of the original matrix, in
                        the same storage format.

        LDA     (input) INTEGER
                        The leading dimension of the array A.  LDA >= max(1,N).
        INFO    (output) INTEGER
                        = 0: successful exit
                        < 0: if INFO = -i, the i-th argument had an illegal value
                        > 0: if INFO = i, A(i,i) is exactly zero.  The triangular
                                matrix is singular and its inverse can not be computed.

        ===================================================================== */


        /* Local variables */
        char uplo_[2] = {uplo, 0};
        char diag_[2] = {diag, 0};
        magma_int_t        ldda, nb, nn;
        static magma_int_t j, jb;
        cuDoubleComplex    zone   = MAGMA_Z_ONE;
        cuDoubleComplex    mzone  = MAGMA_Z_NEG_ONE;
        cuDoubleComplex    *work;
        
        long int    upper  = lapackf77_lsame(uplo_, "U");
        long int    nounit = lapackf77_lsame(diag_, "N");
        
        *info = 0;

        if ((! upper) && (! lapackf77_lsame(uplo_, "L")))
                *info = -1;
        else if ((! nounit) && (! lapackf77_lsame(diag_, "U")))
                *info = -2;
        else if (n < 0)
                *info = -3;
        else if (lda < max(1,n))
                *info = -5;

        if (*info != 0) {
                magma_xerbla( __func__, -(*info) );
                return MAGMA_ERR_ILLEGAL_VALUE;
        }

        /* Quick return */
        if ( n == 0 )
                return MAGMA_SUCCESS;


        /*  Check for singularity if non-unit */
        if (nounit)
        { 
                for (*info=0; *info < n; *info=*info+1)
                {
                        if(A(*info,*info)==0)
                                return MAGMA_ERR_ILLEGAL_VALUE;
                }
                *info=0;
        }


        /* Determine the block size for this environment */
        ldda = ((n+31)/32)*32;

        if (CUBLAS_STATUS_SUCCESS != cublasAlloc((n)*ldda, sizeof( cuDoubleComplex), (void**)&work))
        {
                *info = -6;
                return MAGMA_ERR_CUBLASALLOC;
        }  

        static cudaStream_t stream[2];
        cudaStreamCreate(&stream[0]);
        cudaStreamCreate(&stream[1]);

        nb = magma_get_zpotrf_nb(n);
        
        if (nb <= 1 || nb >= n)
                lapackf77_ztrtri(uplo_, diag_, &n, a, &lda, info);
        else
        {
                if (upper)
                {
                        /* Compute inverse of upper triangular matrix */
                        for (j=0; j<n; j =j+ nb)
                        {
                                jb = min(nb, (n-j));
                                cublasSetMatrix(jb, (n-j), sizeof( cuDoubleComplex), A(j, j), lda, dA(j, j), ldda);

                                /* Compute rows 1:j-1 of current block column */
                                cublasZtrmm(MagmaLeft, MagmaUpper,
                                                        MagmaNoTrans, MagmaNonUnit, j, jb,
                                                        zone, dA(0,0), ldda, dA(0, j),ldda);

                                cublasZtrsm(MagmaRight, MagmaUpper,
                                                        MagmaNoTrans, MagmaNonUnit, j, jb,
                                                        mzone, dA(j,j), ldda, dA(0, j),ldda);

                                //cublasGetMatrix(j ,jb, sizeof( cuDoubleComplex),
                                //dA(0, j), ldda, A(0, j), lda);

                                cudaMemcpy2DAsync( A(j, j), lda *sizeof( cuDoubleComplex),
                                                   dA(j, j), ldda*sizeof( cuDoubleComplex),
                                                   sizeof( cuDoubleComplex)*jb, jb,
                                                   cudaMemcpyDeviceToHost, stream[1]);


                                cudaMemcpy2DAsync( A(0, j), lda *sizeof( cuDoubleComplex),
                                                   dA(0, j), ldda*sizeof( cuDoubleComplex),
                                                   sizeof( cuDoubleComplex)*j, jb,
                                                   cudaMemcpyDeviceToHost, stream[0]);

                                cudaStreamSynchronize(stream[1]);
                        
                                /* Compute inverse of current diagonal block */
                                lapackf77_ztrtri(MagmaUpperStr, diag_, &jb, A(j,j), &lda, info);

                                cublasSetMatrix(jb, jb, sizeof( cuDoubleComplex),
                                                A(j, j), lda, dA(j, j), ldda);
                        }

                }
                else
                {
                        /* Compute inverse of lower triangular matrix */
                        nn=((n-1)/nb)*nb+1;

                        for(j=nn-1; j>=0; j=j-nb)
                        {

                                jb=min(nb,(n-j));

                                if((j+jb) < n)
                                {
                                        cublasSetMatrix((n-j), jb, sizeof( cuDoubleComplex),
                                                A(j, j), lda, dA(j, j), ldda);

                                        /* Compute rows j+jb:n of current block column */
                                        cublasZtrmm(MagmaLeft, MagmaLower,
                                                        MagmaNoTrans, MagmaNonUnit, (n-j-jb), jb,
                                                        zone, dA(j+jb,j+jb), ldda, dA(j+jb, j), ldda);

                                        cublasZtrsm(MagmaRight, MagmaLower,
                                                        MagmaNoTrans, MagmaNonUnit, (n-j-jb), jb,
                                                        mzone, dA(j,j), ldda, dA(j+jb, j), ldda);

                                        //cublasGetMatrix((n-j), jb, sizeof( cuDoub
                                        //leComplex),dA(j, j), ldda, A(j, j), lda);

                                        cudaMemcpy2DAsync( A(j+jb, j),  lda *sizeof( cuDoubleComplex),
                                                           dA(j+jb, j), ldda*sizeof( cuDoubleComplex),
                                                           sizeof( cuDoubleComplex)*(n-j-jb),jb,
                                                           cudaMemcpyDeviceToHost,stream[1]);

                                        cudaMemcpy2DAsync( A(j,j),  lda *sizeof( cuDoubleComplex),
                                                           dA(j,j), ldda*sizeof( cuDoubleComplex),
                                                           sizeof( cuDoubleComplex)*jb, jb,
                                                           cudaMemcpyDeviceToHost,stream[0]);

                                        cudaStreamSynchronize(stream[0]);
                                }

                                /* Compute inverse of current diagonal block */
                                lapackf77_ztrtri(MagmaLowerStr, diag_, &jb, A(j,j), &lda, info);

                                cublasSetMatrix(jb, jb, sizeof( cuDoubleComplex),
                                                A(j, j), lda, dA(j, j), ldda);
                        }
                }
        }

        cudaStreamDestroy(stream[0]);
        cudaStreamDestroy(stream[1]);

        cublasFree(work);

        return MAGMA_SUCCESS;
}
