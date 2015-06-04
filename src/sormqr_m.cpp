/*
    -- MAGMA (version 1.3.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2012

       @author Stan Tomov

       @generated s Wed Nov 14 22:53:17 2012

*/
#include "common_magma.h"

#define  A(i, j) ( a+(j)*lda  + (i))
#define  C(i, j) ( c+(j)*ldc  + (i))

#define dC(gpui, i, j) (dw[gpui]+(j)*lddc + (i))
#define dA_c(gpui, ind, i, j) (dw[gpui] + n_l*lddc + (ind)*lddar*lddac + (i) + (j)*lddac)
#define dA_r(gpui, ind, i, j) (dw[gpui] + n_l*lddc + (ind)*lddar*lddac + (i) + (j)*lddar)
#define dt(gpui, ind)    (dw[gpui] + n_l*lddc + 2*lddac*lddar + (ind)*((nb+1)*nb))
#define dwork(gpui, ind) (dw[gpui] + n_l*lddc + 2*lddac*lddar + 2*((nb+1)*nb) + (ind)*(lddwork*nb))

extern"C"{
  void magmablas_ssetdiag1subdiag0_stream(char uplo, int k, int nb, float *A, int lda, cudaStream_t stream);
}

extern "C" magma_int_t
magma_sormqr_m(magma_int_t nrgpu, char side, char trans,
               magma_int_t m, magma_int_t n, magma_int_t k,
               float *a,    magma_int_t lda,
               float *tau,
               float *c,    magma_int_t ldc,
               float *work, magma_int_t lwork,
               magma_int_t *info)
{
/*  -- MAGMA (version 1.3.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2012

    Purpose
    =======
    SORMQR overwrites the general real M-by-N matrix C with

                    SIDE = 'L'     SIDE = 'R'
    TRANS = 'N':      Q * C          C * Q
    TRANS = 'T':      Q**T * C       C * Q**T

    where Q is a real orthogonal matrix defined as the product of k
    elementary reflectors

          Q = H(1) H(2) . . . H(k)

    as returned by SGEQRF. Q is of order M if SIDE = 'L' and of order N
    if SIDE = 'R'.

    Arguments
    =========
    SIDE    (input) CHARACTER*1
            = 'L': apply Q or Q**T from the Left;
            = 'R': apply Q or Q**T from the Right.

    TRANS   (input) CHARACTER*1
            = 'N':  No transpose, apply Q;
            = 'T':  Transpose, apply Q**T.

    M       (input) INTEGER
            The number of rows of the matrix C. M >= 0.

    N       (input) INTEGER
            The number of columns of the matrix C. N >= 0.

    K       (input) INTEGER
            The number of elementary reflectors whose product defines
            the matrix Q.
            If SIDE = 'L', M >= K >= 0;
            if SIDE = 'R', N >= K >= 0.

    A       (input) REAL array, dimension (LDA,K)
            The i-th column must contain the vector which defines the
            elementary reflector H(i), for i = 1,2,...,k, as returned by
            SGEQRF in the first k columns of its array argument A.

    LDA     (input) INTEGER
            The leading dimension of the array A.
            If SIDE = 'L', LDA >= max(1,M);
            if SIDE = 'R', LDA >= max(1,N).

    TAU     (input) REAL array, dimension (K)
            TAU(i) must contain the scalar factor of the elementary
            reflector H(i), as returned by SGEQRF.

    C       (input/output) REAL array, dimension (LDC,N)
            On entry, the M-by-N matrix C.
            On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q.

    LDC     (input) INTEGER
            The leading dimension of the array C. LDC >= max(1,M).

    WORK    (workspace/output) REAL array, dimension (MAX(1,LWORK))
            On exit, if INFO = 0, WORK(0) returns the optimal LWORK.

    LWORK   (input) INTEGER
            The dimension of the array WORK.
            If SIDE = 'L', LWORK >= max(1,N);
            if SIDE = 'R', LWORK >= max(1,M).
            For optimum performance LWORK >= N*NB if SIDE = 'L', and
            LWORK >= M*NB if SIDE = 'R', where NB is the optimal
            blocksize.

            If LWORK = -1, then a workspace query is assumed; the routine
            only calculates the optimal size of the WORK array, returns
            this value as the first entry of the WORK array, and no error
            message related to LWORK is issued by XERBLA.

    INFO    (output) INTEGER
            = 0:  successful exit
            < 0:  if INFO = -i, the i-th argument had an illegal value
    =====================================================================   */
    float c_one = MAGMA_S_ONE;

    char side_[2] = {side, 0};
    char trans_[2] = {trans, 0};

    float* dw[MagmaMaxGPUs];
    cudaStream_t stream [MagmaMaxGPUs][2];

    magma_int_t ind_c, kb;

    magma_int_t i__4;
    magma_int_t i;
    float t[4160];        /* was [65][64] */
    magma_int_t i1, i2, i3, ib, nb, nq, nw;
    magma_int_t left, notran, lquery;
    magma_int_t iinfo, lwkopt;

    magma_int_t igpu = 0;

    int gpu_b;
    magma_getdevice(&gpu_b);

    *info = 0;
    left = lapackf77_lsame(side_, "L");
    notran = lapackf77_lsame(trans_, "N");
    lquery = (lwork == -1);

    /* NQ is the order of Q and NW is the minimum dimension of WORK */
    if (left) {
        nq = m;
        nw = n;
    } else {
        nq = n;
        nw = m;
    }
    if (! left && ! lapackf77_lsame(side_, "R")) {
        *info = -1;
    } else if (! notran && ! lapackf77_lsame(trans_, "T")) {
        *info = -2;
    } else if (m < 0) {
        *info = -3;
    } else if (n < 0) {
        *info = -4;
    } else if (k < 0 || k > nq) {
        *info = -5;
    } else if (lda < max(1,nq)) {
        *info = -7;
    } else if (ldc < max(1,m)) {
        *info = -10;
    } else if (lwork < max(1,nw) && ! lquery) {
        *info = -12;
    }

    if (*info == 0)
    {
        /* Determine the block size.  NB may be at most NBMAX, where NBMAX
         is used to define the local array T.    */
        nb = 64;
        lwkopt = max(1,nw) * nb;
        MAGMA_S_SET2REAL( work[0], lwkopt );
    }

    if (*info != 0) {
        magma_xerbla( __func__, -(*info) );
        return *info;
    }
    else if (lquery) {
        return *info;
    }

    /* Quick return if possible */
    if (m == 0 || n == 0 || k == 0) {
        work[0] = c_one;
        return *info;
    }

    magma_int_t lddc = (m+63)/64*64;
    magma_int_t lddac = nq;
    magma_int_t lddar =nb;
    magma_int_t lddwork = nw;

    magma_int_t n_l = (n+nrgpu-1)/nrgpu; // local n
    n_l = ((n_l+63)/64)*64;

    if (n_l<256)
       n_l=256;

    nrgpu = min(nrgpu, (n+n_l-1)/n_l); // Don't use GPU that will not have data.

    for (igpu = 0; igpu < nrgpu; ++igpu){
        magma_setdevice(igpu);
        magmablasSetKernelStream(NULL);
        if (MAGMA_SUCCESS != magma_smalloc( &dw[igpu], (n_l*lddc + 2*lddac*lddar + 2*(nb + 1 + lddwork)*nb))) {
            magma_xerbla( __func__, -(*info) );
            *info = MAGMA_ERR_DEVICE_ALLOC;
            return *info;
        }
        magma_queue_create( &stream[igpu][0] );
        magma_queue_create( &stream[igpu][1] );
    }

    if (nb >= k)
    {
        /* Use CPU code */
        lapackf77_sormqr(side_, trans_, &m, &n, &k, a, &lda, tau,
                         c, &ldc, work, &lwork, &iinfo);
    }
    else
    {
        /* Use hybrid CPU-MGPU code */
        if (left) {

            //copy C to mgpus
            for (igpu = 0; igpu < nrgpu; ++igpu){
                magma_setdevice(igpu);
                kb = min(n_l, n-igpu*n_l);
                magma_ssetmatrix_async( m, kb,
                                        C(0, igpu*n_l), ldc,
                                        dC(igpu, 0, 0), lddc, stream[igpu][0] );
            }

            if ( !notran ) {
                i1 = 0;
                i2 = k;
                i3 = nb;
            } else {
                i1 = (k - 1) / nb * nb;
                i2 = 0;
                i3 = -nb;
            }

            kb = min(nb, k-i1);
            for (igpu = 0; igpu < nrgpu; ++igpu){
                magma_setdevice(igpu);
                magma_ssetmatrix_async( nq-i1, kb,
                                        A(i1, i1),            lda,
                                        dA_c(igpu, 0, i1, 0), lddac, stream[igpu][0] );
            }
            ind_c = 0;

            for (i = i1; i3 < 0 ? i >= i2 : i < i2; i += i3)
            {
                ib = min(nb, k - i);
                /* Form the triangular factor of the block reflector
                   H = H(i) H(i+1) . . . H(i+ib-1) */
                i__4 = nq - i;
                lapackf77_slarft("F", "C", &i__4, &ib, A(i, i), &lda,
                                 &tau[i], t, &ib);

                /* H or H' is applied to C(1:m,i:n) */

                /* Apply H or H'; First copy T to the GPU */
                for (igpu = 0; igpu < nrgpu; ++igpu){
                    magma_setdevice(igpu);
                    magma_ssetmatrix_async( ib, ib,
                                            t,               ib,
                                            dt(igpu, ind_c), ib, stream[igpu][ind_c] );

                    magma_queue_sync( stream[igpu][ind_c] ); // Makes sure that we can change t next iteration.
                }

                // start the copy of next A panel
                kb = min(nb, k - i - i3);
                if (kb > 0 && i+i3 >= 0){
                    for (igpu = 0; igpu < nrgpu; ++igpu){
                        magma_setdevice(igpu);
                        magma_ssetmatrix_async( nq-(i+i3), kb,
                                                A(i+i3, i+i3),                    lda,
                                                dA_c(igpu, (ind_c+1)%2, i+i3, 0), lddac, stream[igpu][(ind_c+1)%2] );
                    }
                }

                for (igpu = 0; igpu < nrgpu; ++igpu){
                    magma_setdevice(igpu);
                    // Put 0s in the upper triangular part of dA;
                    magmablas_ssetdiag1subdiag0_stream('L', ib, ib, dA_c(igpu, ind_c, i, 0), lddac, stream[igpu][ind_c]);

                    kb = min(n_l, n-igpu*n_l);
                    magmablasSetKernelStream(stream[igpu][ind_c]);
                    magma_slarfb_gpu( side, trans, MagmaForward, MagmaColumnwise,
                                     m-i, kb, ib,
                                     dA_c(igpu, ind_c, i, 0), lddac, dt(igpu, ind_c), ib,
                                     dC(igpu, i, 0), lddc,
                                     dwork(igpu, ind_c), lddwork);
                }

                ind_c = (ind_c+1)%2;
            }

            //copy C from mgpus
            for (igpu = 0; igpu < nrgpu; ++igpu){
                magma_setdevice(igpu);
                magma_queue_sync( stream[igpu][0] );
                magma_queue_sync( stream[igpu][1] );
                kb = min(n_l, n-igpu*n_l);
                //asynchronous copy gives problems sometimes...
//                magma_sgetmatrix_async( m, kb,
//                                        dC(igpu, 0, 0), lddc,
//                                        C(0, igpu*n_l), ldc, stream[igpu][0] );
                magma_sgetmatrix( m, kb,
                                  dC(igpu, 0, 0), lddc,
                                  C(0, igpu*n_l), ldc );
            }

        } else {

            fprintf(stderr, "The case (side == right) is not implemented\n");
            magma_xerbla( __func__, 1 );
            return *info;

            /*if ( notran ) {
                i1 = 0;
                i2 = k;
                i3 = nb;
            } else {
                i1 = (k - 1) / nb * nb;
                i2 = 0;
                i3 = -nb;
            }

            mi = m;
            ic = 0;

            for (i = i1; i3 < 0 ? i >= i2 : i < i2; i += i3)
            {
                ib = min(nb, k - i);

                // Form the triangular factor of the block reflector
                // H = H(i) H(i+1) . . . H(i+ib-1)
                i__4 = nq - i;
                lapackf77_slarft("F", "C", &i__4, &ib, A(i, i), &lda,
                                 &tau[i], t, &ib);

                // 1) copy the panel from A to the GPU, and
                // 2) Put 0s in the upper triangular part of dA;
                magma_ssetmatrix( i__4, ib, A(i, i), lda, dA(i, 0), ldda );
                magmablas_ssetdiag1subdiag0('L', ib, ib, dA(i, 0), ldda);


                // H or H' is applied to C(1:m,i:n)
                ni = n - i;
                jc = i;

                // Apply H or H'; First copy T to the GPU
                magma_ssetmatrix( ib, ib, t, ib, dt, ib );
                magma_slarfb_gpu( side, trans, MagmaForward, MagmaColumnwise,
                                 mi, ni, ib,
                                 dA(i, 0), ldda, dt, ib,
                                 dC(ic, jc), lddc,
                                 dwork, lddwork);
            }
            */
        }
    }
    MAGMA_S_SET2REAL( work[0], lwkopt );

    for (igpu = 0; igpu < nrgpu; ++igpu){
        magma_setdevice(igpu);
        magma_queue_sync( stream[igpu][0] );
        magmablasSetKernelStream(NULL);
        magma_queue_destroy( stream[igpu][0] );
        magma_queue_destroy( stream[igpu][1] );
        magma_free( dw[igpu] );
    }

    magma_setdevice(gpu_b);

    return *info;
} /* magma_sormqr */


