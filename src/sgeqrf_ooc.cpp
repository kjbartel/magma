/*
    -- MAGMA (version 1.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2011

       @generated s Sun Nov 13 20:48:19 2011

*/
#include "common_magma.h"

extern "C" magma_int_t
magma_sgeqrf_ooc(magma_int_t m, magma_int_t n, 
                 float *a,    magma_int_t lda, float *tau, 
                 float *work, magma_int_t lwork,
                 magma_int_t *info )
{
/*  -- MAGMA (version 1.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2011

    Purpose
    =======
    SGEQRF_OOC computes a QR factorization of a REAL M-by-N matrix A:
    A = Q * R. This version does not require work space on the GPU
    passed as input. GPU memory is allocated in the routine.
    This is an out-of-core (ooc) version that is similar to magma_sgeqrf but
    the difference is that this version can use a GPU even if the matrix
    does not fit into the GPU memory at once.

    Arguments
    =========
    M       (input) INTEGER
            The number of rows of the matrix A.  M >= 0.

    N       (input) INTEGER
            The number of columns of the matrix A.  N >= 0.

    A       (input/output) REAL array, dimension (LDA,N)
            On entry, the M-by-N matrix A.
            On exit, the elements on and above the diagonal of the array
            contain the min(M,N)-by-N upper trapezoidal matrix R (R is
            upper triangular if m >= n); the elements below the diagonal,
            with the array TAU, represent the orthogonal matrix Q as a
            product of min(m,n) elementary reflectors (see Further
            Details).

            Higher performance is achieved if A is in pinned memory, e.g.
            allocated using cudaMallocHost.

    LDA     (input) INTEGER
            The leading dimension of the array A.  LDA >= max(1,M).

    TAU     (output) REAL array, dimension (min(M,N))
            The scalar factors of the elementary reflectors (see Further
            Details).

    WORK    (workspace/output) REAL array, dimension (MAX(1,LWORK))
            On exit, if INFO = 0, WORK(1) returns the optimal LWORK.

            Higher performance is achieved if WORK is in pinned memory, e.g.
            allocated using cudaMallocHost.

    LWORK   (input) INTEGER
            The dimension of the array WORK.  LWORK >= N*NB,
            where NB can be obtained through magma_get_sgeqrf_nb(M).

            If LWORK = -1, then a workspace query is assumed; the routine
            only calculates the optimal size of the WORK array, returns
            this value as the first entry of the WORK array, and no error
            message related to LWORK is issued.

    INFO    (output) INTEGER
            = 0:  successful exit
            < 0:  if INFO = -i, the i-th argument had an illegal value
                  if INFO = -8, the GPU memory allocation failed

    Further Details
    ===============
    The matrix Q is represented as a product of elementary reflectors

       Q = H(1) H(2) . . . H(k), where k = min(m,n).

    Each H(i) has the form

       H(i) = I - tau * v * v'

    where tau is a real scalar, and v is a real vector with
    v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i),
    and tau in TAU(i).
    =====================================================================    */

    #define  a_ref(a_1,a_2) ( a+(a_2)*(lda) + (a_1))
    #define da_ref(a_1,a_2) (da+(a_2)*ldda  + (a_1))

    float *da, *dwork;
    float c_one = MAGMA_S_ONE;

    int  k, lddwork, ldda;

    *info = 0;
    int nb = magma_get_sgeqrf_nb(min(m, n));

    int lwkopt = n * nb;
    work[0] = MAGMA_S_MAKE( (float)lwkopt, 0 );
    long int lquery = (lwork == -1);
    if (m < 0) {
        *info = -1;
    } else if (n < 0) {
        *info = -2;
    } else if (lda < max(1,m)) {
        *info = -4;
    } else if (lwork < max(1,n) && ! lquery) {
        *info = -7;
    }
    if (*info != 0) {
        magma_xerbla( __func__, -(*info) );
        return MAGMA_ERR_ILLEGAL_VALUE;
    }
    else if (lquery)
        return MAGMA_SUCCESS;

    /* Check how much memory do we have */
    #if CUDA_VERSION > 3010
        size_t totalMem;
    #else
        unsigned int totalMem;
    #endif

    CUdevice dev;
    cuDeviceGet( &dev, 0);
    cuDeviceTotalMem( &totalMem, dev );
    totalMem /= sizeof(float);

    magma_int_t IB, NB = (magma_int_t)(0.8*totalMem/m);
    NB = (NB / nb) * nb;

    if (NB >= n)
      return magma_sgeqrf(m, n, a, lda, tau, work, lwork, info);

    k = min(m,n);
    if (k == 0) {
        work[0] = c_one;
        return MAGMA_SUCCESS;
    }

    lddwork = ((NB+31)/32)*32+nb;
    ldda    = ((m+31)/32)*32;

    if (CUBLAS_STATUS_SUCCESS != cublasAlloc((NB + nb)*ldda + nb*lddwork, 
                                             sizeof(float), (void**)&da) ) {
        *info = -8;
        return MAGMA_ERR_CUBLASALLOC;
    }

    static cudaStream_t stream[2];
    cudaStreamCreate(&stream[0]);
    cudaStreamCreate(&stream[1]);

    //   cublasSetKernelStream(stream[1]);

    float *ptr = da + ldda * NB;
    dwork = da + ldda*(NB + nb);

    /* start the main loop over the blocks that fit in the GPU memory */
    for(int i=0; i<n; i+=NB)
      { 
        IB = min(n-i, NB);
        //printf("Processing %5d columns -- %5d to %5d ... \n", IB, i, i+IB);

        /* 1. Copy the next part of the matrix to the GPU */
        cudaMemcpy2DAsync(da_ref(0,0), ldda*sizeof(float),
                          a_ref(0,i), lda *sizeof(float),
                          sizeof(float)*(m), IB,
                          cudaMemcpyHostToDevice,stream[0]);
        cudaStreamSynchronize(stream[0]);

        /* 2. Update it with the previous transformations */
        for(int j=0; j<min(i,k); j+=nb)
          {
            magma_int_t ib = min(k-j, nb);

            /* Get a panel in ptr.                                           */
            //   1. Form the triangular factor of the block reflector
            //   2. Send it to the GPU.
            //   3. Put 0s in the upper triangular part of V.
            //   4. Send V to the GPU in ptr.
            //   5. Update the matrix.
            //   6. Restore the upper part of V.
            int rows = m-j;
            lapackf77_slarft( MagmaForwardStr, MagmaColumnwiseStr,
                              &rows, &ib, a_ref(j,j), &lda, tau+j, work, &ib);
            cudaMemcpy2DAsync(dwork, lddwork *sizeof(float),
                              work,  ib      *sizeof(float),
                              sizeof(float)*ib, ib,
                              cudaMemcpyHostToDevice,stream[1]);

            spanel_to_q(MagmaUpper, ib, a_ref(j,j), lda, work+ib*ib);
            cudaMemcpy2DAsync(ptr,        rows *sizeof(float),
                              a_ref(j,j), lda  *sizeof(float),
                              sizeof(float)*rows, ib, 
                              cudaMemcpyHostToDevice,stream[1]);
            cudaStreamSynchronize(stream[1]);

            magma_slarfb_gpu( MagmaLeft, MagmaTrans, MagmaForward, MagmaColumnwise,
                              rows, IB, ib,
                              ptr, rows, dwork,    lddwork,
                              da_ref(j, 0), ldda, dwork+ib, lddwork);

            sq_to_panel(MagmaUpper, ib, a_ref(j,j), lda, work+ib*ib);
          }

        /* 3. Do a QR on the current part */
        if (i<k)
          magma_sgeqrf2_gpu(m-i, IB, da_ref(i,0), ldda, tau+i, info);

        /* 4. Copy the current part back to the CPU */
        cudaMemcpy2DAsync( a_ref(0,i), lda *sizeof(float),
                          da_ref(0,0), ldda*sizeof(float),
                          sizeof(float)*(m), IB,
                          cudaMemcpyDeviceToHost,stream[0]);
      }

    cudaStreamSynchronize(stream[0]);

    cudaStreamDestroy( stream[0] );
    cudaStreamDestroy( stream[1] );
    cublasFree( da );

    return MAGMA_SUCCESS;
} /* magma_sgeqrf_ooc */

