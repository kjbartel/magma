/*
    -- MAGMA (version 1.6.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2013
       
       @author Azzam Haidar

       @generated from zgesv_rbt_batched.cpp normal z -> c, Sat Nov 15 19:54:10 2014
*/
#include "common_magma.h"
#include "batched_kernel_param.h"
/**
    Purpose
    -------
    Solves a system of linear equations
      A * X = B,  A**T * X = B,  or  A**H * X = B
    with a general N-by-N matrix A using the LU factorization computed by CGETRF_GPU.

    Arguments
    ---------
    @param[in]
    trans   magma_trans_t
            Specifies the form of the system of equations:
      -     = MagmaNoTrans:    A    * X = B  (No transpose)
      -     = MagmaTrans:      A**T * X = B  (Transpose)
      -     = MagmaConjTrans:  A**H * X = B  (Conjugate transpose)

    @param[in]
    n       INTEGER
            The order of the matrix A.  N >= 0.

    @param[in]
    nrhs    INTEGER
            The number of right hand sides, i.e., the number of columns
            of the matrix B.  NRHS >= 0.

    @param[in]
    dA      COMPLEX array on the GPU, dimension (LDA,N)
            The factors L and U from the factorization A = P*L*U as computed
            by CGETRF_GPU.

    @param[in]
    ldda    INTEGER
            The leading dimension of the array A.  LDA >= max(1,N).

    @param[in]
    ipiv    INTEGER array, dimension (N)
            The pivot indices from CGETRF; for 1 <= i <= N, row i of the
            matrix was interchanged with row IPIV(i).

    @param[in,out]
    dB      COMPLEX array on the GPU, dimension (LDB,NRHS)
            On entry, the right hand side matrix B.
            On exit, the solution matrix X.

    @param[in]
    lddb    INTEGER
            The leading dimension of the array B.  LDB >= max(1,N).

    @param[out]
    info    INTEGER
      -     = 0:  successful exit
      -     < 0:  if INFO = -i, the i-th argument had an illegal value

    @ingroup magma_cgesv_comp
    ********************************************************************/
extern "C" magma_int_t
magma_cgesv_rbt_batched(
                  magma_int_t n, magma_int_t nrhs,
                  magmaFloatComplex **dA_array, magma_int_t ldda,
                  magmaFloatComplex **dB_array, magma_int_t lddb,
                  magma_int_t *info_array,
                  magma_int_t batchCount)
{
    /* Local variables */
    
    magma_int_t info;
    info = 0;
    if (n < 0) {
        info = -1;
    } else if (nrhs < 0) {
        info = -2;
    } else if (ldda < max(1,n)) {
        info = -4;
    } else if (lddb < max(1,n)) {
        info = -6;
    }
    if (info != 0) {
        magma_xerbla( __func__, -(info) );
        return info;
    }


    /* Quick return if possible */
    if (n == 0 || nrhs == 0) {
        return info;
    }

    magmaFloatComplex *hu, *hv;
    if (MAGMA_SUCCESS != magma_cmalloc_cpu( &hu, 2*n )) {
        info = MAGMA_ERR_HOST_ALLOC;
        return info;
    }

    if (MAGMA_SUCCESS != magma_cmalloc_cpu( &hv, 2*n )) {
        info = MAGMA_ERR_HOST_ALLOC;
        return info;
    }



    info = magma_cgerbt_batched(MagmaTrue, n, nrhs, dA_array, n, dB_array, n, hu, hv, &info, batchCount);
    if (info != MAGMA_SUCCESS)  {
        return info;
    }


    info = magma_cgetrf_nopiv_batched( n, n, dA_array, ldda, info_array, batchCount);
    if ( (info != MAGMA_SUCCESS) ){
        return info;
    }

#ifdef CHECK_INFO
    // check correctness of results throught "dinfo_magma" and correctness of argument throught "info"
    magma_int_t *cpu_info = (magma_int_t*) malloc(batchCount*sizeof(magma_int_t));
    magma_getvector( batchCount, sizeof(magma_int_t), dinfo_array, 1, cpu_info, 1);
    for(int i=0; i<batchCount; i++)
    {
        if(cpu_info[i] != 0 ){
            printf("magma_cgetrf_batched matrix %d returned error %d\n",i, (int)cpu_info[i] );
            info = cpu_info[i];
            free (cpu_info);
            return info;
        }
    }
    free (cpu_info);
#endif

    info = magma_cgetrs_nopiv_batched( MagmaNoTrans, n, nrhs, dA_array, ldda, dB_array, lddb, info_array, batchCount );


    /* The solution of A.x = b is Vy computed on the GPU */
    magmaFloatComplex *dv;

    if (MAGMA_SUCCESS != magma_cmalloc( &dv, 2*n )) {
        info = MAGMA_ERR_DEVICE_ALLOC;
        return info;
    }

    magma_csetvector(2*n, hv, 1, dv, 1);

    for(int i = 0; i < nrhs; i++)
        magmablas_cprbt_mv_batched(n, dv, dB_array+(i), batchCount);

 //   magma_cgetmatrix(n, nrhs, db, nn, B, ldb);


    return info;
}
