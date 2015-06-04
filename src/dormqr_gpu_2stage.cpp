/*
    -- MAGMA (version 1.5.0-beta3) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date July 2014

       @author Azzam Haidar
       @author Stan Tomov
       @author Raffaele Solca

       @generated from zunmqr_gpu_2stage.cpp normal z -> d, Fri Jul 18 17:34:19 2014

*/
#include "common_magma.h"

/**
    Purpose
    -------
    DORMQR_GPU overwrites the general real M-by-N matrix C with

    @verbatim
                    SIDE = MagmaLeft     SIDE = MagmaRight
    TRANS = MagmaNoTrans:      Q * C          C * Q
    TRANS = MagmaTrans:      Q**T * C       C * Q**T
    @endverbatim

    where Q is a real unitary matrix defined as the product of k
    elementary reflectors

        Q = H(1) H(2) . . . H(k)

    as returned by DGEQRF. Q is of order M if SIDE = MagmaLeft and of order N
    if SIDE = MagmaRight.

    Arguments
    ---------
    @param[in]
    side    magma_side_t
      -      = MagmaLeft:      apply Q or Q**T from the Left;
      -      = MagmaRight:     apply Q or Q**T from the Right.

    @param[in]
    trans   magma_trans_t
      -     = MagmaNoTrans:    No transpose, apply Q;
      -     = MagmaTrans:      Transpose, apply Q**T.

    @param[in]
    m       INTEGER
            The number of rows of the matrix C. M >= 0.

    @param[in]
    n       INTEGER
            The number of columns of the matrix C. N >= 0.

    @param[in]
    k       INTEGER
            The number of elementary reflectors whose product defines
            the matrix Q.
            If SIDE = MagmaLeft,  M >= K >= 0;
            if SIDE = MagmaRight, N >= K >= 0.

    @param[in]
    dA      DOUBLE_PRECISION array on the GPU, dimension (LDDA,K)
            The i-th column must contain the vector which defines the
            elementary reflector H(i), for i = 1,2,...,k, as returned by
            DGEQRF in the first k columns of its array argument DA.
            DA is modified by the routine but restored on exit.

    @param[in]
    ldda    INTEGER
            The leading dimension of the array DA.
            If SIDE = MagmaLeft,  LDDA >= max(1,M);
            if SIDE = MagmaRight, LDDA >= max(1,N).

    @param[in,out]
    dC      DOUBLE_PRECISION array on the GPU, dimension (LDDC,N)
            On entry, the M-by-N matrix C.
            On exit, C is overwritten by Q*C or Q**T * C or C * Q**T or C*Q.

    @param[in]
    lddc     INTEGER
            The leading dimension of the array DC. LDDC >= max(1,M).

    @param[in]
    dT      DOUBLE_PRECISION array on the GPU that is the output
            (the 9th argument) of magma_dgeqrf_gpu.

    @param[in]
    nb      INTEGER
            This is the blocking size that was used in pre-computing DT, e.g.,
            the blocking size used in magma_dgeqrf_gpu.

    @param[out]
    info    INTEGER
      -     = 0:  successful exit
      -     < 0:  if INFO = -i, the i-th argument had an illegal value

    @ingroup magma_dsyev_2stage
    ********************************************************************/
extern "C" magma_int_t
magma_dormqr_gpu_2stages(magma_side_t side, magma_trans_t trans,
                         magma_int_t m, magma_int_t n, magma_int_t k,
                         double *dA,   magma_int_t ldda,
                         double *dC,    magma_int_t lddc,
                         double *dT,    magma_int_t nb,
                         magma_int_t *info)
{
    #define dA(i_,j_) (dA + (i_) + (j_)*ldda)
    #define dC(i_,j_) (dC + (i_) + (j_)*lddc)
    
    double *dwork;

    magma_int_t i, i1, i2, step, ib, ic, jc, mi, ni, nq, nw, ret;
    int left, notran;
    //magma_int_t lwkopt;

    *info = 0;
    left   = (side == MagmaLeft);
    notran = (trans == MagmaNoTrans);

    /* NQ is the order of Q and NW is the minimum dimension of WORK */
    if (left) {
        nq = m;
        nw = n;
    } else {
        nq = n;
        nw = m;
    }
    if ( ! left && side != MagmaRight ) {
        *info = -1;
    } else if ( ! notran && trans != MagmaTrans ) {
        *info = -2;
    } else if (m < 0) {
        *info = -3;
    } else if (n < 0) {
        *info = -4;
    } else if (k < 0 || k > nq) {
        *info = -5;
    } else if (ldda < max(1,nq)) {
        *info = -7;
    } else if (lddc < max(1,m)) {
        *info = -10;
    }

    // TODO alloc after xerbla & quick return, else memory leak
    if (MAGMA_SUCCESS != magma_dmalloc( &dwork, n*nb )) {
        printf ("!!!! dorgqr_2stage magma_alloc failed for: dwork\n" );
        exit(-1);
    }

    if (*info != 0) {
        magma_xerbla( __func__, -(*info) );
        return *info;
    }

    /* Quick return if possible */
    if (m == 0 || n == 0 || k == 0) {
        return *info;
    }

    if ( (left && (! notran)) || ( (! left) && notran ) ) {
        i1 = 0;
        i2 = k;
        step = nb;
    } else {
        i1 = (k - 1) / nb * nb;
        i2 = 0;
        step = -nb;
    }

    // silence "uninitialized" warnings
    mi = 0;
    ni = 0;
    
    if (left) {
        ni = n;
        jc = 0;
    } else {
        mi = m;
        ic = 0;
    }

    for (i=i1; (step < 0 ? i >= i2 : i < i2); i += step) {
        ib = min(nb, k - i);
        if (left) {
            mi = m - i;
            ic = i;
        }
        else {
            ni = n - i;
            jc = i;
        }
        // TODO use info instead of ret, so info & return value always agree.
        ret = magma_dlarfb_gpu( MagmaLeft, trans, MagmaForward, MagmaColumnwise,
                               mi, ni, ib, dA(i,i), ldda, dT+i*nb, nb,
                               dC(ic,jc), lddc, dwork, nw);

        if ( ret != MAGMA_SUCCESS ) {
            magma_free(dwork);
            return ret;
        }
    }
    
    // TODO fix free(dwork) memory leak

    return MAGMA_SUCCESS;
} /* magma_dormqr_gpu_2stages */
