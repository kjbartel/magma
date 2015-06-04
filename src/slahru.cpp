/*
    -- MAGMA (version 1.3.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2012

       @generated s Wed Nov 14 22:53:33 2012

*/
#include "common_magma.h"

// === Define what BLAS to use ============================================
#define PRECISION_s
//#if (defined(PRECISION_s) || defined(PRECISION_d))
#if (defined(PRECISION_d))
 #define magma_sgemm magmablas_sgemm
#endif
// === End defining what BLAS to use =======================================

extern "C" magma_int_t
magma_slahru(magma_int_t n, magma_int_t ihi, magma_int_t k, magma_int_t nb, 
             float *a, magma_int_t lda,
             float *d_a, float *y,
             float *v, float *d_t, 
             float *d_work)
{
/*  -- MAGMA auxiliary routine (version 1.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2012

    Purpose
    =======
    SLAHRU is an auxiliary MAGMA routine that is used in SGEHRD to update
    the trailing sub-matrices after the reductions of the corresponding
    panels.
    See further details below.

    Arguments
    =========
    N       (input) INTEGER
            The order of the matrix A.  N >= 0.

    IHI     (input) INTEGER
            Last row to update. Same as IHI in sgehrd.

    K       (input) INTEGER
            Number of rows of the matrix M (see details below)

    NB      (input) INTEGER
            Block size

    A       (output) REAL array, dimension (LDA,N-K)
            On entry, the N-by-(N-K) general matrix to be updated. The
            computation is done on the GPU. After M is updated on the GPU
            only M(1:NB) is transferred to the CPU - to update the
            corresponding M matrix. See Further Details below.

    LDA     (input) INTEGER
            The leading dimension of the array A.  LDA >= max(1,N).

    D_A     (input/output) REAL array on the GPU, dimension
            (N,N-K). On entry, the N-by-(N-K) general matrix to be updated.
            On exit, the 1st K rows (matrix M) of A are updated by
            applying an orthogonal transformation from the right
            M = M (I-V T V'), and sub-matrix G is updated by
            G = (I - V T V') G (I - V T V(NB+1:)' )
            where Q = I - V T V' represent the orthogonal matrix
            (as a product of elementary reflectors V) used to reduce
            the current panel of A to upper Hessenberg form. After M
            is updated M(:,1:NB) is sent to the CPU.
            See Further Details below.

    Y       (input/workspace) REAL array on the GPU, dimension
            (N, NB). On entry the N-K-by-NB Y = A V. It is used internally
            as workspace, so its value is changed on exit.

    V       (input/workspace) REAL array onthe GPU, dimension
            (N, NB). On entry the N-K-by-NB matrix V of elementary reflectors
            used to reduce the current panel of A to upper Hessenberg form.
            The rest K-by-NB part is used as workspace. V is unchanged on
            exit.

    D_T     (input) REAL array on the GPU, dimension (NB, NB).
            On entry the NB-by-NB upper trinagular matrix defining the
            orthogonal Hessenberg reduction transformation matrix for
            the current panel. The lower triangular part are 0s.

    D_WORK  (workspace) REAL array on the GPU, dimension N*NB.

    Further Details
    ===============
    This implementation follows the algorithm and notations described in:

    S. Tomov and J. Dongarra, "Accelerating the reduction to upper Hessenberg
    form through hybrid GPU-based computing," University of Tennessee Computer
    Science Technical Report, UT-CS-09-642 (also LAPACK Working Note 219),
    May 24, 2009.

    The difference is that here M is computed on the GPU.
    =====================================================================    */

    float c_zero    = MAGMA_S_ZERO;
    float c_one     = MAGMA_S_ONE;
    float c_neg_one = MAGMA_S_NEG_ONE;

    magma_int_t ldda = lda;
    float *v0 = v + ihi - k;

    /* V0 = M V */
    magma_sgemm( MagmaNoTrans, MagmaNoTrans, k, nb, ihi-k,
                 c_one,  d_a, ldda,
                         v,   ldda,
                 c_zero, v0,  ldda);

    /* Update matrix M -= V0 T V' through
       1. d_work = T V'
       2. M -= V0 d_work                  */
    magma_sgemm( MagmaNoTrans, MagmaTrans, nb, ihi-k, nb,
                 c_one,  d_t,    nb,
                         v,      ldda,
                 c_zero, d_work, nb);

    magma_sgemm( MagmaNoTrans, MagmaNoTrans, k, ihi-k, nb,
                 c_neg_one, v0,     ldda,
                            d_work, nb,
                 c_one,     d_a,    ldda);
    magma_sgetmatrix( k, nb, d_a, ldda, a, lda );

    /* Update G -= Y T -= Y d_work */
    magma_sgemm( MagmaNoTrans, MagmaNoTrans, ihi-k, ihi-k-nb, nb,
                 c_neg_one, y,                ldda,
                            d_work+nb*nb,     nb,
                 c_one,     d_a   +nb*ldda+k, ldda);

    /* Update G2 = (I - V T V') G2 = (I - work' V') G2 through
       1. Y = V' G2
       2. G2 -= work' Y
       Note that G is A(k:ihi, nb+1:ihi-k)
       while    G2 is A(k:ihi, nb+1: n -k)   */
    magma_sgemm( MagmaTrans, MagmaNoTrans, nb, n-k-nb, ihi-k,
                 c_one,  v,               ldda,
                         d_a + nb*ldda+k, ldda,
                 c_zero, y,               nb);
    magma_sgemm( MagmaTrans, MagmaNoTrans, ihi-k, n-k-nb, nb,
                 c_neg_one, d_work,        nb,
                            y,             nb,
                 c_one,     d_a+nb*ldda+k, ldda);
    return 0;
}
