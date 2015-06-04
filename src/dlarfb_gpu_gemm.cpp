/*
    -- MAGMA (version 1.5.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date September 2014

       @author Mark Gates
       @author Azzam Haidar
       @generated from zlarfb_gpu_gemm.cpp normal z -> d, Wed Sep 17 15:08:32 2014
*/
#include "common_magma.h"

/**
    Purpose
    -------
    DLARFB applies a real block reflector H or its transpose H^H to a
    DOUBLE_PRECISION m by n matrix C, from the left.
    
    __Note that this function assumes__ that the upper part of dV is 0
    because it is referenced. Same for upper/lower part of dT.

    Arguments
    ---------
    @param[in]
    side    magma_side_t
      -     = MagmaLeft:      apply H or H^H from the Left
      -     = MagmaRight:     apply H or H^H from the Right

    @param[in]
    trans   magma_trans_t
      -     = MagmaNoTrans:    apply H   (No transpose)
      -     = MagmaTrans: apply H^H (Conjugate transpose)

    @param[in]
    direct  magma_direct_t
            Indicates how H is formed from a product of elementary
            reflectors
      -     = MagmaForward:  H = H(1) H(2) . . . H(k) (Forward)
      -     = MagmaBackward: H = H(k) . . . H(2) H(1) (Backward)

    @param[in]
    storev  magma_storev_t
            Indicates how the vectors which define the elementary
            reflectors are stored:
      -     = MagmaColumnwise: Columnwise
      -     = MagmaRowwise:    Rowwise

    @param[in]
    m       INTEGER
            The number of rows of the matrix C.

    @param[in]
    n       INTEGER
            The number of columns of the matrix C.

    @param[in]
    k       INTEGER
            The order of the matrix T (= the number of elementary
            reflectors whose product defines the block reflector).

    @param[in]
    dV      DOUBLE_PRECISION array on the GPU, dimension
                (LDV,K) if STOREV = MagmaColumnwise
                (LDV,M) if STOREV = MagmaRowwise and SIDE = MagmaLeft
                (LDV,N) if STOREV = MagmaRowwise and SIDE = MagmaRight
            The matrix V. See further details.

    @param[in]
    ldv     INTEGER
            The leading dimension of the array V.
            If STOREV = MagmaColumnwise and SIDE = MagmaLeft, LDV >= max(1,M);
            if STOREV = MagmaColumnwise and SIDE = MagmaRight, LDV >= max(1,N);
            if STOREV = MagmaRowwise, LDV >= K.

    @param[in]
    dT      DOUBLE_PRECISION array on the GPU, dimension (LDT,K)
            The triangular k by k matrix T in the representation of the
            block reflector.

    @param[in]
    ldt     INTEGER
            The leading dimension of the array T. LDT >= K.

    @param[in,out]
    dC      DOUBLE_PRECISION array on the GPU, dimension (LDC,N)
            On entry, the m by n matrix C.
            On exit, C is overwritten by H*C, or H^H*C, or C*H, or C*H^H.

    @param[in]
    ldc     INTEGER
            The leading dimension of the array C. LDA >= max(1,M).

    @param
    dwork   (workspace) DOUBLE_PRECISION array, dimension (LDWORK,K)

    @param[in]
    ldwork  INTEGER
            The leading dimension of the array WORK.
            If SIDE = MagmaLeft,  LDWORK >= max(1,N);
            if SIDE = MagmaRight, LDWORK >= max(1,M);

    @param
    dworkvt (workspace) DOUBLE_PRECISION array, dimension (LDWORKT,K)

    @param[in]
    ldworkvt INTEGER
            The leading dimension of the array WORKVT.
            LDWORKVT >= max(1,min(M,N));

    Further Details
    ---------------
    The shape of the matrix V and the storage of the vectors which define
    the H(i) is best illustrated by the following example with n = 5 and
    k = 3.
    All elements including 0's and 1's are stored, unlike LAPACK.

        DIRECT = MagmaForward and         DIRECT = MagmaForward and
        STOREV = MagmaColumnwise:         STOREV = MagmaRowwise:

                 V = (  1  0  0 )                 V = (  1 v1 v1 v1 v1 )
                     ( v1  1  0 )                     (  0  1 v2 v2 v2 )
                     ( v1 v2  1 )                     (  0  0  1 v3 v3 )
                     ( v1 v2 v3 )
                     ( v1 v2 v3 )

        DIRECT = MagmaBackward and        DIRECT = MagmaBackward and 
        STOREV = MagmaColumnwise:         STOREV = MagmaRowwise:

                 V = ( v1 v2 v3 )                 V = ( v1 v1  1  0  0 )
                     ( v1 v2 v3 )                     ( v2 v2 v2  1  0 )
                     (  1 v2 v3 )                     ( v3 v3 v3 v3  1 )
                     (  0  1 v3 )
                     (  0  0  1 )

    @ingroup magma_daux3
    ********************************************************************/
extern "C" magma_int_t
magma_dlarfb_gpu_gemm( magma_side_t side, magma_trans_t trans, magma_direct_t direct, magma_storev_t storev,
                  magma_int_t m, magma_int_t n, magma_int_t k,
                  const double *dV,    magma_int_t ldv,
                  const double *dT,    magma_int_t ldt,
                  double *dC,          magma_int_t ldc,
                  double *dwork,       magma_int_t ldwork,
                  double *dworkvt,     magma_int_t ldworkvt)
{
    double c_zero    = MAGMA_D_ZERO;
    double c_one     = MAGMA_D_ONE;
    double c_neg_one = MAGMA_D_NEG_ONE;

    /* Function Body */
    magma_int_t info = 0;
    if (m <= 0 || n <= 0) {
        return info;
    }
    //internal variable
    magma_int_t ldwvt = m > n ?  k : m;
    magma_int_t ldw;
    if ( side == MagmaLeft ) {
        ldw = k;
    } else {
        ldw = m;
    }
    // opposite of trans
    magma_trans_t transt;
    if (trans == MagmaNoTrans)
        transt = MagmaTrans;
    else
        transt = MagmaNoTrans;
    
    // whether T is upper or lower triangular
    magma_uplo_t uplo;
    if (direct == MagmaForward)
        uplo = MagmaUpper;
    else
        uplo = MagmaLower;
    
    // whether V is stored transposed or not
    magma_trans_t notransV, transV;
    if (storev == MagmaColumnwise) {
        notransV = MagmaNoTrans;
        transV   = MagmaTrans;
    }
    else {
        notransV = MagmaTrans;
        transV   = MagmaNoTrans;
    }

    if ( side == MagmaLeft ) {
        // Form H C or H^H C
        // Comments assume H C.
        // When forming H^H C, T gets transposed via transt for m >= n or by trans for m < n.
        
        // W = V' C
        magma_dgemm( MagmaTrans,notransV,
                     k, n, m,
                     c_one,  dV,    ldv,
                             dC,    ldc,
                     c_zero, dwork, ldw);

        if (m <= n) {
            // W2 = V T
            magma_dgemm( notransV, trans,
                         m, k, k,
                         c_one,  dV, ldv,
                                 dT, ldt,
                         c_zero, dworkvt, ldwvt);
            // C = C - W2 W = C - V T V' C = (I - V T V') C = H C
            magma_dgemm( MagmaNoTrans, MagmaNoTrans,
                         m, n, k,
                         c_neg_one, dworkvt,  ldwvt,
                                    dwork,    ldw,
                         c_one,     dC,       ldc);
        } else {
            // W2 = T W  = T  V' C
            magma_dgemm( trans, MagmaNoTrans,
                         k, n, k,
                         c_one,  dT, ldt,
                                 dwork, ldw,
                         c_zero, dworkvt, ldwvt);
            // C = C - V W2 = C - V T V' C = (I - V T V') C = H C
            magma_dgemm( notransV, MagmaNoTrans,
                         m, n, k,
                         c_neg_one, dV,  ldv,
                                    dworkvt,  ldwvt,
                         c_one,     dC,       ldc);
        }
    }
    else {
        // Form C H or C H^H
        // Comments assume C H.
        // When forming C H^H, T gets transposed via trans.
        
        // W = C V
        magma_dgemm( MagmaNoTrans, notransV,
                     m, k, n,
                     c_one,  dC,    ldc,
                             dV,    ldv,
                     c_zero, dwork, ldw);
        if (m <= n) {
            // W2 = W T = C V T
            magma_dgemm( MagmaNoTrans, trans,
                         m, k, k,
                         c_one,  dwork, ldw,
                                 dT, ldt,
                         c_zero, dworkvt, ldwvt);
            // C = C - W2 V' = C - C V T V' = C (I - V T V') = C H
            magma_dgemm( MagmaNoTrans, transV,
                         m, n, k,
                         c_neg_one, dworkvt, ldwvt,
                                    dV,    ldv,
                         c_one,     dC,    ldc);
        } else {
            // W2 = T V'
            magma_dgemm( trans, transV,
                         k, n, k,
                         c_one,  dT, ldt,
                                 dV, ldv,
                         c_zero, dworkvt, ldwvt);
            // C = C - W W2 = C - C V T V' = C (I - V T V') = C H
            magma_dgemm( MagmaNoTrans, MagmaNoTrans,
                         m, n, k,
                         c_neg_one, dwork,   ldw,
                                    dworkvt, ldwvt,
                         c_one,     dC,      ldc);
        }
    }

    return info;
} /* magma_dlarfb */
