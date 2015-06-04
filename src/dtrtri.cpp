/*
    -- MAGMA (version 1.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2011

       @generated d Sun Nov 13 20:48:32 2011

*/
#include "common_magma.h"

// === Define what BLAS to use ============================================
 #define PRECISION_d
 #if (defined(PRECISION_s) || defined(PRECISION_d))
   #define cublasDgemm magmablas_dgemm
     #define cublasDtrsm magmablas_dtrsm
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
magma_dtrtri(char uplo, char diag, magma_int_t n,
              double *a, magma_int_t lda, magma_int_t *info)
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

        A       (input/output) DOUBLE_PRECISION array, dimension (LDA,N)
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
        double    zone   = MAGMA_D_ONE;
        double    mzone  = MAGMA_D_NEG_ONE;
        double    *work;
        
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

        if (CUBLAS_STATUS_SUCCESS != cublasAlloc((n)*ldda, sizeof( double), (void**)&work))
        {
                *info = -6;
                return MAGMA_ERR_CUBLASALLOC;
        }  

        static cudaStream_t stream[2];
        cudaStreamCreate(&stream[0]);
        cudaStreamCreate(&stream[1]);

        nb = magma_get_dpotrf_nb(n);
        
        if (nb <= 1 || nb >= n)
                lapackf77_dtrtri(uplo_, diag_, &n, a, &lda, info);
        else
        {
                if (upper)
                {
                        /* Compute inverse of upper triangular matrix */
                        for (j=0; j<n; j =j+ nb)
                        {
                                jb = min(nb, (n-j));
                                cublasSetMatrix(jb, (n-j), sizeof( double), A(j, j), lda, dA(j, j), ldda);

                                /* Compute rows 1:j-1 of current block column */
                                cublasDtrmm(MagmaLeft, MagmaUpper,
                                                        MagmaNoTrans, MagmaNonUnit, j, jb,
                                                        zone, dA(0,0), ldda, dA(0, j),ldda);

                                cublasDtrsm(MagmaRight, MagmaUpper,
                                                        MagmaNoTrans, MagmaNonUnit, j, jb,
                                                        mzone, dA(j,j), ldda, dA(0, j),ldda);

                                //cublasGetMatrix(j ,jb, sizeof( double),
                                //dA(0, j), ldda, A(0, j), lda);

                                cudaMemcpy2DAsync( A(j, j), lda *sizeof( double),
                                                   dA(j, j), ldda*sizeof( double),
                                                   sizeof( double)*jb, jb,
                                                   cudaMemcpyDeviceToHost, stream[1]);


                                cudaMemcpy2DAsync( A(0, j), lda *sizeof( double),
                                                   dA(0, j), ldda*sizeof( double),
                                                   sizeof( double)*j, jb,
                                                   cudaMemcpyDeviceToHost, stream[0]);

                                cudaStreamSynchronize(stream[1]);
                        
                                /* Compute inverse of current diagonal block */
                                lapackf77_dtrtri(MagmaUpperStr, diag_, &jb, A(j,j), &lda, info);

                                cublasSetMatrix(jb, jb, sizeof( double),
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
                                        cublasSetMatrix((n-j), jb, sizeof( double),
                                                A(j, j), lda, dA(j, j), ldda);

                                        /* Compute rows j+jb:n of current block column */
                                        cublasDtrmm(MagmaLeft, MagmaLower,
                                                        MagmaNoTrans, MagmaNonUnit, (n-j-jb), jb,
                                                        zone, dA(j+jb,j+jb), ldda, dA(j+jb, j), ldda);

                                        cublasDtrsm(MagmaRight, MagmaLower,
                                                        MagmaNoTrans, MagmaNonUnit, (n-j-jb), jb,
                                                        mzone, dA(j,j), ldda, dA(j+jb, j), ldda);

                                        //cublasGetMatrix((n-j), jb, sizeof( cuDoub
                                        //leComplex),dA(j, j), ldda, A(j, j), lda);

                                        cudaMemcpy2DAsync( A(j+jb, j),  lda *sizeof( double),
                                                           dA(j+jb, j), ldda*sizeof( double),
                                                           sizeof( double)*(n-j-jb),jb,
                                                           cudaMemcpyDeviceToHost,stream[1]);

                                        cudaMemcpy2DAsync( A(j,j),  lda *sizeof( double),
                                                           dA(j,j), ldda*sizeof( double),
                                                           sizeof( double)*jb, jb,
                                                           cudaMemcpyDeviceToHost,stream[0]);

                                        cudaStreamSynchronize(stream[0]);
                                }

                                /* Compute inverse of current diagonal block */
                                lapackf77_dtrtri(MagmaLowerStr, diag_, &jb, A(j,j), &lda, info);

                                cublasSetMatrix(jb, jb, sizeof( double),
                                                A(j, j), lda, dA(j, j), ldda);
                        }
                }
        }

        cudaStreamDestroy(stream[0]);
        cudaStreamDestroy(stream[1]);

        cublasFree(work);

        return MAGMA_SUCCESS;
}
