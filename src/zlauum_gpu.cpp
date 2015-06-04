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
magma_zlauum_gpu(char uplo, magma_int_t n,
             cuDoubleComplex  *dA, magma_int_t ldda, magma_int_t *info)
{


/*  -- MAGMA (version 1.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2011

        Purpose
        =======

        DLAUUM computes the product U * U' or L' * L, where the triangular
        factor U or L is stored in the upper or lower triangular part of
        the array dA.

        If UPLO = 'U' or 'u' then the upper triangle of the result is stored,
        overwriting the factor U in dA.
        If UPLO = 'L' or 'l' then the lower triangle of the result is stored,
        overwriting the factor L in dA.
        This is the blocked form of the algorithm, calling Level 3 BLAS.

        Arguments
        =========

        UPLO    (input) CHARACTER*1
                        Specifies whether the triangular factor stored in the array dA
                        is upper or lower triangular:
                        = 'U':  Upper triangular
                        = 'L':  Lower triangular

        N       (input) INTEGER
                        The order of the triangular factor U or L.  N >= 0.

        dA       (input/output) DOUBLE PRECISION array on the GPU, dimension (LDDA,N)
                        On entry, the triangular factor U or L.
                        On exit, if UPLO = 'U', the upper triangle of dA is
                        overwritten with the upper triangle of the product U * U';
                        if UPLO = 'L', the lower triangle of dA is overwritten with
                        the lower triangle of the product L' * L.

        LDDA     (input) INTEGER
                        The leading dimension of the array A.  LDDA >= max(1,N).

        INFO    (output) INTEGER
                        = 0: successful exit
                        < 0: if INFO = -k, the k-th argument had an illegal value

        ===================================================================== */



        /* Local variables */
        char uplo_[2] = {uplo, 0};
        magma_int_t         nb, i, ib;
        double              done   = MAGMA_D_ONE;
        cuDoubleComplex     zone = MAGMA_Z_ONE;
        cuDoubleComplex     *work;

        long int upper  = lapackf77_lsame(uplo_, "U");

        *info = 0;

        if ((! upper) && (! lapackf77_lsame(uplo_, "L")))
                *info = -1;
        else if (n < 0)
                *info = -2;
        else if (ldda < max(1,n))
                *info = -4;

        if (*info != 0) {
                magma_xerbla( __func__, -(*info) );
                return MAGMA_ERR_ILLEGAL_VALUE;
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
                lapackf77_zlauum(uplo_, &n, work, &n, info);
                cublasSetMatrix(n, n, sizeof(cuDoubleComplex), work, n, dA, ldda);
        }
        else
        {
                if (upper)
                {
                        /* Compute inverse of upper triangular matrix */
                        for (i=0; i<n; i =i+ nb)
                        {
                                ib = min(nb, (n-i));

                                /* Compute the product U * U'. */
                                cublasZtrmm( MagmaRight, MagmaUpper,
                                             MagmaConjTrans, MagmaNonUnit, i, ib,
                                             zone, dA(i,i), ldda, dA(0, i),ldda);

                                cublasGetMatrix( ib ,ib, sizeof(cuDoubleComplex),
                                                 dA(i, i), ldda, work, ib);
                                                
                                lapackf77_zlauum(MagmaUpperStr, &ib, work, &ib, info);

                                cublasSetMatrix( ib, ib, sizeof(cuDoubleComplex),
                                                 work, ib, dA(i, i), ldda);

                                if(i+ib < n)
                                {
                                        cublasZgemm( MagmaNoTrans, MagmaConjTrans,
                                                     i, ib, (n-i-ib), zone, dA(0,i+ib),
                                                     ldda, dA(i, i+ib),ldda, zone,
                                                     dA(0,i), ldda);


                                        cublasZherk( MagmaUpper, MagmaNoTrans, ib,(n-i-ib),
                                                     done, dA(i, i+ib), ldda,
                                                     done,  dA(i, i), ldda);
                                }
                        }

                }
                else
                {
                        /* Compute the product L' * L. */
                        for(i=0; i<n; i=i+nb)
                        {
                                ib=min(nb,(n-i));

                                cublasZtrmm( MagmaLeft, MagmaLower,
                                             MagmaConjTrans, MagmaNonUnit, ib,
                                             i, zone, dA(i,i), ldda,
                                             dA(i, 0),ldda);
                                
                                cublasGetMatrix( ib ,ib, sizeof(cuDoubleComplex),
                                                 dA(i, i), ldda, work, ib);

                                lapackf77_zlauum(MagmaLowerStr, &ib, work, &ib, info);

                                cublasSetMatrix( ib, ib, sizeof(cuDoubleComplex),
                                                 work, ib, dA(i, i), ldda);
                                

                                if((i+ib) < n)
                                {
                                        cublasZgemm( MagmaConjTrans, MagmaNoTrans,
                                                     ib, i, (n-i-ib), zone, dA( i+ib,i),
                                                     ldda, dA(i+ib, 0),ldda, zone,
                                                     dA(i,0), ldda);
                                        cublasZherk( MagmaLower, MagmaConjTrans, ib, (n-i-ib),
                                                     done, dA(i+ib, i), ldda,
                                                     done,  dA(i, i), ldda);
                                }
                        }
                }
        }

        cudaStreamDestroy(stream[0]);
        cudaStreamDestroy(stream[1]);

        cudaFreeHost(work);

        return MAGMA_SUCCESS;
}
