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
magma_zlauum(char uplo, magma_int_t n,
             cuDoubleComplex *a, magma_int_t lda, magma_int_t *info)
{
/*  -- MAGMA (version 1.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2011

        Purpose
        =======

        ZLAUUM computes the product U * U' or L' * L, where the triangular
        factor U or L is stored in the upper or lower triangular part of
        the array A.

        If UPLO = 'U' or 'u' then the upper triangle of the result is stored,
        overwriting the factor U in A.
        If UPLO = 'L' or 'l' then the lower triangle of the result is stored,
        overwriting the factor L in A.
        This is the blocked form of the algorithm, calling Level 3 BLAS.

        Arguments
        =========

        UPLO    (input) CHARACTER*1
                        Specifies whether the triangular factor stored in the array A
                        is upper or lower triangular:
                        = 'U':  Upper triangular
                        = 'L':  Lower triangular

        N       (input) INTEGER
                        The order of the triangular factor U or L.  N >= 0.

        A       (input/output) COPLEX_16 array, dimension (LDA,N)
                        On entry, the triangular factor U or L.
                        On exit, if UPLO = 'U', the upper triangle of A is
                        overwritten with the upper triangle of the product U * U';
                        if UPLO = 'L', the lower triangle of A is overwritten with
                        the lower triangle of the product L' * L.

        LDA     (input) INTEGER
                        The leading dimension of the array A.  LDA >= max(1,N).

        INFO    (output) INTEGER
                        = 0: successful exit
                        < 0: if INFO = -k, the k-th argument had an illegal value

        ===================================================================== */


        /* Local variables */
        char uplo_[2] = {uplo, 0};
        magma_int_t     ldda, nb;
        static magma_int_t i, ib;
        cuDoubleComplex    zone  = MAGMA_Z_ONE;
        double             done  = MAGMA_D_ONE;
        cuDoubleComplex    *work;
        long int           upper = lapackf77_lsame(uplo_, "U");

        *info = 0;
        if ((! upper) && (! lapackf77_lsame(uplo_, "L")))
                *info = -1;
        else if (n < 0)
                *info = -2;
        else if (lda < max(1,n))
                *info = -4;

        if (*info != 0) {
                magma_xerbla( __func__, -(*info) );
                return MAGMA_ERR_ILLEGAL_VALUE;
        }

        /* Quick return */
        if ( n == 0 )
                return MAGMA_SUCCESS;

        ldda = ((n+31)/32)*32;

        if (CUBLAS_STATUS_SUCCESS != cublasAlloc((n)*ldda, sizeof(cuDoubleComplex), (void**)&work))
        {
                *info = -6;
                return MAGMA_ERR_CUBLASALLOC;
        }

        static cudaStream_t stream[2];
        cudaStreamCreate(&stream[0]);
        cudaStreamCreate(&stream[1]);

        nb = magma_get_zpotrf_nb(n);

        if (nb <= 1 || nb >= n)
                lapackf77_zlauum(uplo_, &n, a, &lda, info);
        else
        {
                if (upper)
                {
                        /* Compute the product U * U'. */
                        for (i=0; i<n; i=i+nb)
                        {
                                ib=min(nb,n-i);

                                //cublasSetMatrix(ib, (n-i), sizeof(cuDoubleComplex), A(i, i), lda, dA(i, i), ldda);
                                
                                cudaMemcpy2DAsync( dA(i, i), ldda *sizeof(cuDoubleComplex),
                                                   A(i,i), lda*sizeof(cuDoubleComplex),
                                                   sizeof(cuDoubleComplex)*ib, ib,
                                                   cudaMemcpyHostToDevice,stream[1]);

                                cudaMemcpy2DAsync( dA(i,i+ib),  ldda *sizeof(cuDoubleComplex),
                                                   A(i,i+ib), lda*sizeof(cuDoubleComplex),
                                                   sizeof(cuDoubleComplex)*ib, (n-i-ib),
                                                   cudaMemcpyHostToDevice,stream[0]);

                                cudaStreamSynchronize(stream[1]);

                                cublasZtrmm( MagmaRight, MagmaUpper,
                                             MagmaConjTrans, MagmaNonUnit, i, ib,
                                             zone, dA(i,i), ldda, dA(0, i),ldda);


                                lapackf77_zlauum(MagmaUpperStr, &ib, A(i,i), &lda, info);

                                cudaMemcpy2DAsync( dA(i, i), ldda * sizeof(cuDoubleComplex),
                                                   A(i, i), lda  * sizeof(cuDoubleComplex),
                                                   sizeof(cuDoubleComplex)*ib, ib,
                                                   cudaMemcpyHostToDevice,stream[0]);

                                if (i+ib < n)
                                {
                                        cublasZgemm( MagmaNoTrans, MagmaConjTrans,
                                                     i, ib, (n-i-ib), zone, dA(0,i+ib),
                                                     ldda, dA(i, i+ib),ldda, zone,
                                                     dA(0,i), ldda);

                                        cudaStreamSynchronize(stream[0]);

                                        cublasZherk( MagmaUpper, MagmaNoTrans, ib,(n-i-ib),
                                                     done, dA(i, i+ib), ldda,
                                                     done,  dA(i, i), ldda);
                                }
                                
                                cublasGetMatrix( i+ib,ib, sizeof(cuDoubleComplex),
                                                 dA(0, i), ldda, A(0, i), lda);
                        }
                }
                else
                {
                        /* Compute the product L' * L. */
                        for(i=0; i<n; i=i+nb)
                        {
                                ib=min(nb,n-i);
                                //cublasSetMatrix((n-i), ib, sizeof(cuDoubleComplex),
                                //                A(i, i), lda, dA(i, i), ldda);

                                cudaMemcpy2DAsync( dA(i, i), ldda *sizeof(cuDoubleComplex),
                                                   A(i,i), lda*sizeof(cuDoubleComplex),
                                                   sizeof(cuDoubleComplex)*ib, ib,
                                                   cudaMemcpyHostToDevice,stream[1]);

                                cudaMemcpy2DAsync( dA(i+ib, i),  ldda *sizeof(cuDoubleComplex),
                                                   A(i+ib, i), lda*sizeof(cuDoubleComplex),
                                                   sizeof(cuDoubleComplex)*(n-i-ib), ib,
                                                   cudaMemcpyHostToDevice,stream[0]);

                                cudaStreamSynchronize(stream[1]);

                                cublasZtrmm( MagmaLeft, MagmaLower,
                                             MagmaConjTrans, MagmaNonUnit, ib,
                                             i, zone, dA(i,i), ldda,
                                             dA(i, 0),ldda);


                                lapackf77_zlauum(MagmaLowerStr, &ib, A(i,i), &lda, info);

                                //cublasSetMatrix(ib, ib, sizeof(cuDoubleComplex),
                                //                A(i, i), lda, dA(i, i), ldda);

                                cudaMemcpy2DAsync( dA(i, i), ldda * sizeof(cuDoubleComplex),
                                                   A(i, i), lda  * sizeof(cuDoubleComplex),
                                                   sizeof(cuDoubleComplex)*ib, ib,
                                                   cudaMemcpyHostToDevice,stream[0]);

                                if (i+ib < n)
                                {
                                        cublasZgemm(MagmaConjTrans, MagmaNoTrans,
                                                        ib, i, (n-i-ib), zone, dA( i+ib,i),
                                                        ldda, dA(i+ib, 0),ldda, zone,
                                                        dA(i,0), ldda);

                                        cudaStreamSynchronize(stream[0]);
                                        
                                        cublasZherk(MagmaLower, MagmaConjTrans, ib, (n-i-ib),
                                                        done, dA(i+ib, i), ldda,
                                                        done,  dA(i, i), ldda);
                                }
                                cublasGetMatrix(ib, i+ib, sizeof(cuDoubleComplex),
                                        dA(i, 0), ldda, A(i, 0), lda);
                        }
                }
        }
        cudaStreamDestroy(stream[0]);
        cudaStreamDestroy(stream[1]);

        cublasFree(work);

        return MAGMA_SUCCESS;

}
