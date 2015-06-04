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

#define dA(i, j) (dA+(j)*ldda + (i))

extern "C" magma_int_t
magma_ztrtri_gpu(char uplo, char diag, magma_int_t n,
             cuDoubleComplex *dA, magma_int_t ldda, magma_int_t *info)
{
/*  -- MAGMA (version 1.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2011

    Purpose
    =======

        DTRTRI computes the inverse of a real upper or lower triangular
        matrix dA.

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

        dA       (input/output) DOUBLE PRECISION array ON THE GPU, dimension (LDDA,N)
                        On entry, the triangular matrix A.  If UPLO = 'U', the
                        leading N-by-N upper triangular part of the array dA contains
                        the upper triangular matrix, and the strictly lower
                        triangular part of A is not referenced.  If UPLO = 'L', the
                        leading N-by-N lower triangular part of the array dA contains
                        the lower triangular matrix, and the strictly upper
                        triangular part of A is not referenced.  If DIAG = 'U', the
                        diagonal elements of A are also not referenced and are
                        assumed to be 1.
                        On exit, the (triangular) inverse of the original matrix, in
                        the same storage format.

        LDDA     (input) INTEGER
                        The leading dimension of the array dA.  LDDA >= max(1,N).
        INFO    (output) INTEGER
                        = 0: successful exit
                        < 0: if INFO = -i, the i-th argument had an illegal value
                        > 0: if INFO = i, dA(i,i) is exactly zero.  The triangular
                                matrix is singular and its inverse can not be computed.

        ===================================================================== */


        /* Local variables */
        char uplo_[2] = {uplo, 0};
        char diag_[2] = {diag, 0};
        magma_int_t         nb, nn, j, jb;
        cuDoubleComplex     zone   = MAGMA_Z_ONE;
        cuDoubleComplex     mzone  = MAGMA_Z_NEG_ONE;
        cuDoubleComplex     *work;

        long int        upper  = lapackf77_lsame(uplo_, "U");
        long int    nounit = lapackf77_lsame(diag_, "N");

        *info = 0;

        if ((! upper) && (! lapackf77_lsame(uplo_, "L")))
                *info = -1;
        else if ((! nounit) && (! lapackf77_lsame(diag_, "U")))
                *info = -2;
        else if (n < 0)
                *info = -3;
        else if (ldda < max(1,n))
                *info = -5;

        if (*info != 0) {
                magma_xerbla( __func__, -(*info) );
                return MAGMA_ERR_ILLEGAL_VALUE;
        }


        /*  Check for singularity if non-unit */
        if (nounit)
        { 
                for (*info=0; *info < n; *info=*info+1)
                {
                        if(dA(*info,*info)==0)
                                return MAGMA_ERR_ILLEGAL_VALUE;
                }
                *info=0;
        }

        nb = magma_get_zpotrf_nb(n);
        
        if (cudaSuccess != cudaMallocHost( (void**)&work, nb*nb*sizeof(cuDoubleComplex) ) ) 
        {
                *info = -6;
                return MAGMA_ERR_HOSTALLOC;
        }
        
        static cudaStream_t stream[2];
        cudaStreamCreate(&stream[0]);
        cudaStreamCreate(&stream[1]);

        
        if (nb <= 1 || nb >= n)
        {
                cublasGetMatrix(n, n, sizeof(cuDoubleComplex), dA, ldda, work, n);
                lapackf77_ztrtri(uplo_, diag_, &n, work, &n, info);
                cublasSetMatrix(n, n, sizeof(cuDoubleComplex), work, n, dA, ldda);
        }
        else
        {
                if (upper)
                {
                        /* Compute inverse of upper triangular matrix */
                        for (j=0; j<n; j =j+ nb)
                        {
                                jb = min(nb, (n-j));

                                /* Compute rows 1:j-1 of current block column */
                                cublasZtrmm(MagmaLeft, MagmaUpper,
                                                        MagmaNoTrans, MagmaNonUnit, j, jb,
                                                        zone, dA(0,0), ldda, dA(0, j),ldda);

                                cublasZtrsm(MagmaRight, MagmaUpper,
                                                        MagmaNoTrans, MagmaNonUnit, j, jb,
                                                        mzone, dA(j,j), ldda, dA(0, j),ldda);

                        
                                //cublasGetMatrix(jb ,jb, sizeof(cuDoubleComplex),
                                //                dA(j, j), ldda, work, jb);

                                cudaMemcpy2DAsync( work,     jb  *sizeof(cuDoubleComplex),
                                                   dA(j, j), ldda*sizeof(cuDoubleComplex),
                                                   jb*sizeof(cuDoubleComplex), jb,
                                                   cudaMemcpyDeviceToHost,stream[1]);

                                cudaStreamSynchronize(stream[1]);

                                /* Compute inverse of current diagonal block */
                                lapackf77_ztrtri(MagmaUpperStr, diag_, &jb, work, &jb, info);

                                //cublasSetMatrix(jb, jb, sizeof(cuDoubleComplex),
                                //                work, jb, dA(j, j), ldda);

                                cudaMemcpy2DAsync( dA(j, j), ldda*sizeof(cuDoubleComplex),
                                                   work,     jb  *sizeof(cuDoubleComplex),
                                                   sizeof(cuDoubleComplex)*jb, jb,
                                                   cudaMemcpyHostToDevice,stream[0]);
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

                                        /* Compute rows j+jb:n of current block column */
                                        cublasZtrmm(MagmaLeft, MagmaLower,
                                                        MagmaNoTrans, MagmaNonUnit, (n-j-jb), jb,
                                                        zone, dA(j+jb,j+jb), ldda, dA(j+jb, j), ldda);

                                        cublasZtrsm(MagmaRight, MagmaLower,
                                                        MagmaNoTrans, MagmaNonUnit, (n-j-jb), jb,
                                                        mzone, dA(j,j), ldda, dA(j+jb, j), ldda);
                                }

                                //cublasGetMatrix(jb, jb, sizeof(cuDoubleComplex),
                                //               dA(j, j), ldda, work, jb);
                                
                                cudaMemcpy2DAsync( work,     jb  *sizeof(cuDoubleComplex),
                                                   dA(j, j), ldda*sizeof(cuDoubleComplex),
                                                   jb*sizeof(cuDoubleComplex), jb,
                                                   cudaMemcpyDeviceToHost,stream[1]);

                                cudaStreamSynchronize(stream[1]);

                                /* Compute inverse of current diagonal block */
                                lapackf77_ztrtri(MagmaLowerStr, diag_, &jb, work, &jb, info);
                        
                                //cublasSetMatrix(jb, jb, sizeof(cuDoubleComplex),
                                //                work, jb, dA(j, j), ldda);

                                cudaMemcpy2DAsync( dA(j, j), ldda*sizeof(cuDoubleComplex),
                                                   work,     jb  *sizeof(cuDoubleComplex),
                                                   sizeof(cuDoubleComplex)*jb, jb,
                                                   cudaMemcpyHostToDevice,stream[0]);
                        }
                }
        }

        cudaStreamDestroy(stream[0]);
        cudaStreamDestroy(stream[1]);

        cudaFreeHost(work);

        return MAGMA_SUCCESS;
}
