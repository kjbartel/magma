/*
    -- MAGMA (version 1.2.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       June 2012

       @author Stan Tomov

       @generated c Thu Jun 28 12:30:47 2012

*/
#include "common_magma.h"

extern "C" magma_int_t
magma_cunmqr_gpu(char side, char trans,
                 magma_int_t m, magma_int_t n, magma_int_t k,
                 cuFloatComplex *dA,    magma_int_t ldda, 
                 cuFloatComplex *tau,
                 cuFloatComplex *dC,    magma_int_t lddc,
                 cuFloatComplex *hwork, magma_int_t lwork,
                 cuFloatComplex *dT,    magma_int_t nb, 
                 magma_int_t *info)
{
/*  -- MAGMA (version 1.2.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       June 2012

    Purpose
    =======
    CUNMQR_GPU overwrites the general complex M-by-N matrix C with

                    SIDE = 'L'     SIDE = 'R'
    TRANS = 'N':      Q * C          C * Q
    TRANS = 'T':      Q**H * C       C * Q**H

    where Q is a complex orthogonal matrix defined as the product of k
    elementary reflectors

          Q = H(1) H(2) . . . H(k)

    as returned by CGEQRF. Q is of order M if SIDE = 'L' and of order N
    if SIDE = 'R'.

    Arguments
    =========
    SIDE    (input) CHARACTER*1
            = 'L': apply Q or Q**H from the Left;
            = 'R': apply Q or Q**H from the Right.

    TRANS   (input) CHARACTER*1
            = 'N':  No transpose, apply Q;
            = 'T':  Transpose, apply Q**H.

    M       (input) INTEGER
            The number of rows of the matrix C. M >= 0.

    N       (input) INTEGER
            The number of columns of the matrix C. N >= 0.

    K       (input) INTEGER
            The number of elementary reflectors whose product defines
            the matrix Q.
            If SIDE = 'L', M >= K >= 0;
            if SIDE = 'R', N >= K >= 0.

    DA      (input) COMPLEX array on the GPU, dimension (LDDA,K)
            The i-th column must contain the vector which defines the
            elementary reflector H(i), for i = 1,2,...,k, as returned by
            CGEQRF in the first k columns of its array argument DA.
            DA is modified by the routine but restored on exit.

    LDDA    (input) INTEGER
            The leading dimension of the array DA.
            If SIDE = 'L', LDDA >= max(1,M);
            if SIDE = 'R', LDDA >= max(1,N).

    TAU     (input) COMPLEX array, dimension (K)
            TAU(i) must contain the scalar factor of the elementary
            reflector H(i), as returned by CGEQRF.

    DC      (input/output) COMPLEX array on the GPU, dimension (LDDC,N)
            On entry, the M-by-N matrix C.
            On exit, C is overwritten by Q*C or Q**H * C or C * Q**H or C*Q.

    LDDC     (input) INTEGER
            The leading dimension of the array DC. LDDC >= max(1,M).

    HWORK    (workspace/output) COMPLEX array, dimension (MAX(1,LWORK))
            On exit, if INFO = 0, HWORK(1) returns the optimal LWORK.

    LWORK   (input) INTEGER
            The dimension of the array HWORK.
            LWORK >= (M-K+NB)*(N+2*NB) if SIDE = 'L', 
            and LWORK >= (N-K+NB)*(M+2*NB) if SIDE = 'R', where NB is the 
            optimal blocksize.

            If LWORK = -1, then a workspace query is assumed; the routine
            only calculates the optimal size of the HWORK array, returns
            this value as the first entry of the HWORK array, and no error
            message related to LWORK is issued by XERBLA.

    DT      (input) COMPLEX array on the GPU that is the output 
            (the 9th argument) of magma_cgeqrf_gpu.

    NB      (input) INTEGER
            This is the blocking size that was used in pre-computing DT, e.g.,
            the blocking size used in magma_cgeqrf_gpu.

    INFO    (output) INTEGER
            = 0:  successful exit
            < 0:  if INFO = -i, the i-th argument had an illegal value
    =====================================================================   */

    #define a_ref(a_1,a_2) (dA+(a_2)*(ldda) + (a_1))
    #define c_ref(a_1,a_2) (dC+(a_2)*(lddc) + (a_1))
    #define t_ref(a_1)     (dT+(a_1)*nb)

    cuFloatComplex c_one = MAGMA_C_ONE;

    char side_[2] = {side, 0};
    char trans_[2] = {trans, 0};

    cuFloatComplex *dwork;
    magma_int_t i, lddwork;

    magma_int_t i1, i2, i3, ib, ic, jc, mi, ni, nq, nw, ret;
    int left, notran, lquery;
    magma_int_t lwkopt;

    *info = 0;
    left   = lapackf77_lsame(side_, "L");
    notran = lapackf77_lsame(trans_, "N");
    lquery = (lwork == -1);

    if (!left || notran)
      printf("cunmqr_gpu called with arguments not yet supported\n");

    /* NQ is the order of Q and NW is the minimum dimension of WORK */
    if (left) {
        nq = m;
        nw = n;
    } else {
        nq = n;
        nw = m;
    }
    if ( (!left) && (!lapackf77_lsame(side_, "R")) ) {
        *info = -1;
    } else if ( (!notran) && (!lapackf77_lsame(trans_, MagmaConjTransStr)) ) {
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
    } else if (lwork < max(1,nw) && ! lquery) {
        *info = -12;
    }

    lwkopt = (m-k+nb)*(n+2*nb);
    hwork[0] = MAGMA_C_MAKE( lwkopt, 0 );

    if (*info != 0) {
        magma_xerbla( __func__, -(*info) );
        return *info;
    }
    else if (lquery) {
        return *info;
    }

    /* Quick return if possible */
    if (m == 0 || n == 0 || k == 0) {
        hwork[0] = c_one;
        return *info;
    }

    lddwork= k;
    dwork  = dT+2*lddwork*nb;

    if ( (left && (! notran)) || ( (!left) && notran ) ) {
        i1 = 0;
        i2 = k-nb;
        i3 = nb;
    } else {
        i1 = (k - 1 - nb) / nb * nb;
        i2 = 0;
        i3 = -nb;
    }

    if (left) {
        ni = n;
        jc = 0;
    } else {
        mi = m;
        ic = 0;
    }

    if (nb < k)
    {
        for (i=i1; i3<0 ? i>i2 : i<i2; i+=i3)
        {
            ib = min(nb, k - i);
            if (left){
                mi = m - i;
                ic = i;
            }
            else {
                ni = n - i;
                jc = i;
            }
            ret = magma_clarfb_gpu( MagmaLeft, MagmaConjTrans, MagmaForward, MagmaColumnwise,
                                    mi, ni, ib, 
                                    a_ref(i,  i ), ldda, t_ref(i), nb, 
                                    c_ref(ic, jc), lddc, dwork,    nw);
            if ( ret != MAGMA_SUCCESS )
              return ret;
        }
    }
    else
    {
        i = i1;
    }

    /* Use unblocked code to multiply the last or only block. */
    if (i < k) {
        ib   = k-i;
        if (left){
            mi = m - i;
            ic = i;
        }
        else {
            ni = n - i;
            jc = i;
        }

        magma_cgetmatrix( mi, ib, a_ref(i,  i ), ldda, hwork, mi );
        magma_cgetmatrix( mi, ni, c_ref(ic, jc), lddc, hwork+mi*ib, mi );

        magma_int_t lhwork = lwork - mi*(ib + ni);
        lapackf77_cunmqr( MagmaLeftStr, MagmaConjTransStr, 
                          &mi, &ni, &ib, 
                          hwork,       &mi, tau+i, 
                          hwork+mi*ib, &mi, 
                          hwork+mi*(ib+ni), &lhwork, info);

        // send the updated part of c back to the GPU
        magma_csetmatrix( mi, ni, hwork+mi*ib, mi, c_ref(ic, jc), lddc );
    }

    return *info;
    /* End of MAGMA_CUNMQR_GPU */
}
