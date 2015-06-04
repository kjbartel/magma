/*
    -- MAGMA (version 1.3.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2012

       @generated s Wed Nov 14 22:53:33 2012

*/
#include "common_magma.h"

extern "C" magma_int_t
magma_sorghr(magma_int_t n, magma_int_t ilo, magma_int_t ihi, 
             float *a, magma_int_t lda, 
             float *tau,
             float *dT, magma_int_t nb,
             magma_int_t *info)
{
/*  -- MAGMA (version 1.3.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2012

    Purpose   
    =======   
    SORGHR generates a REAL unitary matrix Q which is defined as the   
    product of IHI-ILO elementary reflectors of order N, as returned by   
    SGEHRD:   

    Q = H(ilo) H(ilo+1) . . . H(ihi-1).   

    Arguments   
    =========   
    N       (input) INTEGER   
            The order of the matrix Q. N >= 0.   

    ILO     (input) INTEGER   
    IHI     (input) INTEGER   
            ILO and IHI must have the same values as in the previous call   
            of SGEHRD. Q is equal to the unit matrix except in the   
            submatrix Q(ilo+1:ihi,ilo+1:ihi).   
            1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0.   

    A       (input/output) REAL array, dimension (LDA,N)   
            On entry, the vectors which define the elementary reflectors,   
            as returned by SGEHRD.   
            On exit, the N-by-N unitary matrix Q.   

    LDA     (input) INTEGER   
            The leading dimension of the array A. LDA >= max(1,N).   

    TAU     (input) REAL array, dimension (N-1)   
            TAU(i) must contain the scalar factor of the elementary   
            reflector H(i), as returned by SGEHRD.   

    DT      (input) REAL array on the GPU device.
            DT contains the T matrices used in blocking the elementary
            reflectors H(i), e.g., this can be the 9th argument of
            magma_sgehrd.

    NB      (input) INTEGER
            This is the block size used in SGEHRD, and correspondingly
            the size of the T matrices, used in the factorization, and
            stored in DT.

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   
    ===================================================================== */

    #define a_ref(i,j) (a + (j)*lda+ (i))

    magma_int_t i, j, nh, iinfo;

    *info = 0;
    nh = ihi - ilo;
    if (n < 0)
      *info = -1;
    else if (ilo < 1 || ilo > max(1,n)) 
      *info = -2;
    else if (ihi < min(ilo,n) || ihi > n) 
      *info = -3;
    else if (lda < max(1,n)) 
        *info = -5;

    if (*info != 0) {
        magma_xerbla( __func__, -(*info) );
        return *info;
    }

    /* Quick return if possible */
    if (n == 0) 
      return *info;

    /* Shift the vectors which define the elementary reflectors one   
       column to the right, and set the first ilo and the last n-ihi   
       rows and columns to those of the unit matrix */
    for (j = ihi-1; j >= ilo; --j) {
      for (i = 0; i < j; ++i)
        *a_ref(i, j) = MAGMA_S_ZERO;
        
      for (i = j+1; i < ihi; ++i)
        *a_ref(i, j) = *a_ref(i, j - 1);
        
      for (i = ihi; i < n; ++i)
        *a_ref(i, j) = MAGMA_S_ZERO;
    }
    for (j = 0; j < ilo; ++j) {
      for (i = 0; i < n; ++i)
        *a_ref(i, j) = MAGMA_S_ZERO;
        
      *a_ref(j, j) = MAGMA_S_ONE;
    }
    for (j = ihi; j < n; ++j) {
      for (i = 0; i < n; ++i)
        *a_ref(i, j) = MAGMA_S_ZERO; 
        
      *a_ref(j, j) = MAGMA_S_ONE;
    }

    if (nh > 0)
      /* Generate Q(ilo+1:ihi,ilo+1:ihi) */
      magma_sorgqr(nh, nh, nh,
                   a_ref(ilo, ilo), lda,
                   tau+ilo-1, dT, nb, &iinfo);

    return *info;
} /* magma_sorghr */

#undef a_ref
