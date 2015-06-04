/*
    -- MAGMA (version 1.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2011

       @precisions normal z -> s d c

*/

#include "common_magma.h"

/* === Define what BLAS to use ============================================ */
#define PRECISION_z
#if (defined(PRECISION_s) || defined(PRECISION_d))
  #define cublasZgemm magmablas_zgemm
  #define cublasZtrsm magmablas_ztrsm
#endif
/* === End defining what BLAS to use ======================================= */


/* to appy pivoting from the previous big-panel: need some index-adjusting */
extern "C" void
magmablas_zpermute_long3( cuDoubleComplex *dAT, int lda, int *ipiv, int nb, int ind );


extern "C" magma_int_t
magma_zgetrf_ooc(magma_int_t m, magma_int_t n, cuDoubleComplex *a, magma_int_t lda, 
                 magma_int_t *ipiv, magma_int_t *info)
{
/*  -- MAGMA (version 1.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2011

    Purpose
    =======
    ZGETRF_OOC computes an LU factorization of a general M-by-N matrix A
    using partial pivoting with row interchanges.  This version does not
    require work space on the GPU passed as input. GPU memory is allocated
    in the routine. The matrix may not fit entirely in the GPU memory.

    The factorization has the form
       A = P * L * U
    where P is a permutation matrix, L is lower triangular with unit
    diagonal elements (lower trapezoidal if m > n), and U is upper
    triangular (upper trapezoidal if m < n).

    This is the right-looking Level 3 BLAS version of the algorithm.

    Arguments
    =========
    M       (input) INTEGER
            The number of rows of the matrix A.  M >= 0.

    N       (input) INTEGER
            The number of columns of the matrix A.  N >= 0.

    A       (input/output) COMPLEX_16 array, dimension (LDA,N)
            On entry, the M-by-N matrix to be factored.
            On exit, the factors L and U from the factorization
            A = P*L*U; the unit diagonal elements of L are not stored.

            Higher performance is achieved if A is in pinned memory, e.g.
            allocated using cudaMallocHost.

    LDA     (input) INTEGER
            The leading dimension of the array A.  LDA >= max(1,M).

    IPIV    (output) INTEGER array, dimension (min(M,N))
            The pivot indices; for 1 <= i <= min(M,N), row i of the
            matrix was interchanged with row IPIV(i).

    INFO    (output) INTEGER
            = 0:  successful exit
            < 0:  if INFO = -i, the i-th argument had an illegal value
                  if INFO = -7, the GPU memory allocation failed
            > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
                  has been completed, but the factor U is exactly
                  singular, and division by zero will occur if it is used
                  to solve a system of equations.

    =====================================================================    */

#define    A(i,j) (a   + (j)*lda + (i))
#define inAT(i,j) (dAT + (i)*nb*maxn + (j)*nb)
#define inPT(i,j) (dPT + (i)*nb*nb + (j)*nb)

    cuDoubleComplex        *dAT, *dA, *da, *dPT, *work;
    cuDoubleComplex        c_one     = MAGMA_Z_ONE;
    cuDoubleComplex        c_neg_one = MAGMA_Z_NEG_ONE;
    magma_int_t                iinfo, nb, maxm, maxn, maxdim;
    magma_int_t                N, M, NB, NBk, MB, I;
    magma_int_t                i, ii, jj, offset, ib, rows, cols, s, nb0, m0;
#if CUDA_VERSION > 3010
    size_t totalMem;
#else
    unsigned int totalMem;
#endif
    CUdevice dev;

    *info = 0;

    if (m < 0)
        *info = -1;
    else if (n < 0)
        *info = -2;
    else if (lda < max(1,m))
        *info = -4;

    if (*info != 0) {
        magma_xerbla( __func__, -(*info) );
        return MAGMA_ERR_ILLEGAL_VALUE;
    }

    /* Quick return if possible */
    if (m == 0 || n == 0)
        return MAGMA_SUCCESS;

    /* initialize nb */
    nb = magma_get_zgetrf_nb(m);

    /* figure out NB */
    cuDeviceGet( &dev, 0);
    cuDeviceTotalMem( &totalMem, dev );
    totalMem /= sizeof(cuDoubleComplex);
        
    /* printf( " max. matrix dimension (%d)\n",(int)sqrt((double)totalMem) ); */
    MB = m;                                      /* number of rows in the big panel    */
    NB = (magma_int_t)(0.8*totalMem/(2*m))-2*nb; /* number of columns in the big panel */
    if( NB >= n ) {
#ifdef CHECK_ZGETRF_OOC
      printf( "      * still fit in GPU memory.\n" );
#endif
      NB = n;
    } 
#ifdef CHECK_ZGETRF_OOC
    else {
      printf( "      * don't fit in GPU memory.\n" );
    }
#endif
    NB = (NB / nb) * nb;   /* making sure it's devisable by nb   */

#ifdef CHECK_ZGETRF_OOC
        if( NB != n ) printf( "      * running in out-core mode (n=%d, NB=%d, nb=%d).\n",n,NB,nb );
        else          printf( "      * running in in-core mode  (n=%d, NB=%d, nb=%d).\n",n,NB,nb );
    fflush(stdout);
#endif 

    if ( (nb <= 1) || (nb >= min(m,n)) ) {
        /* Use CPU code for scalar of one tile. */
            lapackf77_zgetrf(&m, &n, a, &lda, ipiv, info);
    } else {
        /* Use hybrid blocked code. */

        maxm = ((MB + 31)/32)*32;
        maxn = ((NB + 31)/32)*32;
        maxdim = max(maxm, maxn);

        /* allocate memory on GPU to store the big panel */
        if (CUBLAS_STATUS_SUCCESS != cublasAlloc((2*nb+maxn)*maxm, 
                                                 sizeof(cuDoubleComplex), (void**)&dA) ) {
          *info = -7; 
          return MAGMA_ERR_CUBLASALLOC;
        }
        da  = dA + 2*nb*maxm; /* for transposing the next panel to be sent to CPU */
        dPT = dA +   nb*maxm; /* for storing the previous panel from CPU          */

        /* allocate memory to store the transpose of A */
        if (CUBLAS_STATUS_SUCCESS != cublasAlloc(maxm*maxn, 
                                                 sizeof(cuDoubleComplex), (void**)&dAT) ) {
          cublasFree(dA);
          *info = -7; 
          return MAGMA_ERR_CUBLASALLOC;
        }

        for( I=0; I<n; I+=NB ) {
                  M = MB;
                  N = min( NB, n-I );       /* number of columns in this big panel             */
                  s = min(max(m-I,0),N)/nb; /* number of small block-columns in this big panel */

                  /* upload the next big panel into GPU, transpose (A->A'), and pivot it */
                  cublasSetMatrix( M, N, sizeof(cuDoubleComplex), A(0, I), lda, da, maxm);
                  magmablas_ztranspose2( dAT, maxn, da, maxm, M, N );

                  /* == --------------------------------------------------------------- == */
                  /* == loop around the previous big-panels to update the new big-panel == */
                  for( offset = 0; offset<min(m,I); offset+=NB ) {

                        /* applying the pivot from the big-panel */
                        NBk    = min( m-offset, NB );
                        magmablas_zpermute_long3( dAT, maxn, ipiv, NBk, offset );

                        /* == going through each block-column of this big-panel == */
                        for( jj=0, ib=offset/nb; jj<NBk; jj+=nb, ib++ ) 
                          {
                            nb0  = min(NBk-jj,nb);
                            ii   = offset+jj;
                            rows = maxm - ii;

                            /* upload the previous block-column to GPU */
                            cublasSetMatrix( M-ii, nb, sizeof(cuDoubleComplex), A(ii, ii), lda, dA, rows);
                            magmablas_ztranspose2( dPT, nb, dA, rows, M-ii, nb0);

                            /* update with the block column */
                            cuCtxSynchronize();
                            cublasZtrsm( MagmaRight, MagmaUpper, MagmaNoTrans, MagmaUnit, 
                                         N, nb0, c_one, inPT(0,0), nb, inAT(ib,0), maxn );
                            if( M > ii+nb0 ) {
                              cublasZgemm( MagmaNoTrans, MagmaNoTrans, 
                                           N, M-(ii+nb0), nb0, c_neg_one, inAT(ib,0), maxn, 
                                           inPT(1,0), nb, c_one, inAT(ib+1,0), maxn );
                            }
                            
                          } /* end of for each block-columns in a big-panel */
                        
                  } /* end of for each previous big-panels */

                  nb0 = min( nb, n-I );
                  m0  = M-I;
                  work = &a[I*lda];   /* using the first nb0 columns as the workspace */
                  if( m0 > 0 ) {      /* if more rows to be factorized */

                    /* download the first block-column in this big-panel to CPU */
                    if( I > 0 ) {
                      cols = maxm - I;    /* the number of columns in At */

                      cuCtxSynchronize();
                      magmablas_ztranspose2( dA, cols, inAT(I/nb,0), maxn, nb0, cols );
                      cublasGetMatrix( M-I, nb0, sizeof(cuDoubleComplex), dA, cols, work, lda);
                    }

                    /* factorize the first diagonal block of this big panel; ipiv is 1-base */
                    lapackf77_zgetrf( &m0, &nb0, work, &lda, ipiv+I, &iinfo);
                    if( iinfo != 0 ) {
                          *info = iinfo;
                          break;
                    }

                    /* for each small block-columns in this big panel */
                    for( ii = 0; ii < s; ii++ ) {

                          i = I/nb+ii;         /* row-index of the current diagonal block in global A */
                  cols = maxm - i*nb;  /* the number of columns in At                         */
            
                  if (ii>0) {

                    /* download i-th panel to CPU (into work)                                     */
                            /* dtranspose makes the assumption of the matrix size being a multiple of 32. */
                        magmablas_ztranspose( dA, cols, inAT(i,ii), maxn, nb, cols );
                        cublasGetMatrix( m-i*nb, nb, sizeof(cuDoubleComplex), dA, cols, work, lda);
                
                        /* make sure that gpu queue is empty */
                        cuCtxSynchronize();
                
                            /* update the remaining matrix with (i-1)-th panel */
                        cublasZtrsm( MagmaRight, MagmaUpper, MagmaNoTrans, MagmaUnit, 
                                         N - (ii+1)*nb, nb, 
                                         c_one, inAT(i-1,ii-1), maxn, 
                                         inAT(i-1,ii+1), maxn );
                        cublasZgemm( MagmaNoTrans, MagmaNoTrans, 
                                         N-(ii+1)*nb, M-i*nb, nb, 
                                         c_neg_one, inAT(i-1,ii+1), maxn, 
                                         inAT(i,  ii-1), maxn, 
                                         c_one, inAT(i,  ii+1), maxn );

                        /* do the cpu part; i.e., factorize the i-th panel  */
                        rows = m - i*nb;
                        lapackf77_zgetrf( &rows, &nb, work, &lda, ipiv+i*nb, &iinfo);
                  }
                  if (*info == 0 && iinfo > 0)
                        *info = iinfo + i*nb;

                          /* apply the pivoting from the i-th panel   */
                          /* to the columns in the current big panel  */
                  magmablas_zpermute_long2( dAT, maxn, ipiv, nb, i*nb );

                  /* upload i-th panel to GPU, and transpose it */
                  cublasSetMatrix( m-i*nb, nb, sizeof(cuDoubleComplex), work, lda, dA, cols);
                  magmablas_ztranspose( inAT(i,ii), maxn, dA, cols, cols, nb);

                  /* do the small non-parallel computations;              */
                          /* i.e., update the (i+1)-th column with the i-th panel */
                  if (s > (ii+1)) {
                            cublasZtrsm( MagmaRight, MagmaUpper, MagmaNoTrans, MagmaUnit, nb, nb, 
                                         c_one, inAT(i, ii  ), maxn,   /* diagonal of i-th panel         */
                                    inAT(i, ii+1), maxn);  /* upper-block in (i+1)-th column */
                        cublasZgemm( MagmaNoTrans, MagmaNoTrans, nb, M-(i+1)*nb, nb, 
                                         c_neg_one, inAT(i,   ii+1), maxn,    /* upper-block of (i+1)-th column      */
                                        inAT(i+1, ii  ), maxn,    /* off-diagonal blocks from i-th panel */
                                         c_one,     inAT(i+1, ii+1), maxn );  /* blocks to be updated                */
                  } else {
                            cublasZtrsm( MagmaRight, MagmaUpper, MagmaNoTrans, MagmaUnit, N-s*nb, nb,
                             c_one, inAT(i, ii  ), maxn,
                                    inAT(i, ii+1), maxn);
                        cublasZgemm( MagmaNoTrans, MagmaNoTrans, N-s*nb, M-(i+1)*nb, nb,
                             c_neg_one, inAT(i,   ii+1), maxn,
                                        inAT(i+1, ii  ), maxn, 
                             c_one,     inAT(i+1, ii+1), maxn );
                  }
                } /* end of for i=0,..,s-1 */

                    /* the last off-set */
            i    = I/nb+s;
                    nb0  = min(M - i*nb, N - s*nb);
            rows = M    - i*nb;
            cols = maxm - i*nb;

                    if( nb0 > 0 ) {
                      /* download the last columns to CPU */
              magmablas_ztranspose2( dA, cols, inAT(i,s), maxn, nb0, rows);
              cublasGetMatrix(rows, nb0, sizeof(cuDoubleComplex), dA, cols, work, lda);

              /* make sure that gpu queue is empty */
              cuCtxSynchronize();

              /* do the cpu part; factorize the last column  */
              lapackf77_zgetrf( &rows, &nb0, work, &lda, ipiv+i*nb, &iinfo);
              if (*info == 0 && iinfo > 0)
                *info = iinfo + s*nb;

                      /* apply the pivoting from the last columns to those in GPU */
              magmablas_zpermute_long2( dAT, maxn, ipiv, nb0, i*nb );

                      /* upload the last panel to GPU, and transpose it */
              cublasSetMatrix(rows, nb0, sizeof(cuDoubleComplex), work, lda, dA, cols);
              magmablas_ztranspose2( inAT(i,s), maxn, dA, cols, rows, nb0);

                      /* update with the last (in case the matrix is wide; i.e., n > m). */
              cublasZtrsm( MagmaRight, MagmaUpper, MagmaNoTrans, MagmaUnit, 
                           N-s*nb-nb0, nb0,
                           c_one, inAT(i, s),     maxn, 
                                  inAT(i, s)+nb0, maxn);
                    } /* end of big-panel factorization */
                  } /* end if more row to be factorized */

                  /* download the current big panel to CPU */
          magmablas_ztranspose2( da, maxm, dAT, maxn, N, M );
          cublasGetMatrix( M, N, sizeof(cuDoubleComplex), da, maxm, A(0, I), lda);

            } /* end of for */

        cublasFree(dAT); 
        cublasFree(dA); 
    }
    
    return MAGMA_SUCCESS;
} /* magma_zgetrf_ooc */




extern "C" magma_int_t
magma_zgetrf_piv(magma_int_t m, magma_int_t n, cuDoubleComplex *a, magma_int_t lda, 
                     magma_int_t *ipiv, magma_int_t *info)
{
    magma_int_t nb;
    magma_int_t NB, MB, I, k1, k2, incx, minmn;

    /* Function Body */
    *info = 0;

    if (m < 0)
        *info = -1;
    else if (n < 0)
        *info = -2;
    else if (lda < max(1,m))
        *info = -4;

    if (*info != 0)
        return MAGMA_ERR_ILLEGAL_VALUE;

    /* Quick return if possible */
    if (m == 0 || n == 0)
        return MAGMA_SUCCESS;

        /* initialize nb */
    nb = magma_get_zgetrf_nb(m);

        /* figure out NB */
#if CUDA_VERSION > 3010
    size_t totalMem;
#else
    unsigned int totalMem;
#endif
    CUdevice dev;
    cuDeviceGet( &dev, 0);
    cuDeviceTotalMem( &totalMem, dev );
    totalMem /= sizeof(cuDoubleComplex);
        MB = m;                                             /* number of rows in the big panel    */
    NB = (magma_int_t)min((0.8*totalMem/(2*m))-2*nb,n); /* number of columns in the big panel */
        NB = (NB / nb) * nb;   /* making sure it's devisable by nb   */
    minmn = min(m,n);

        for( I=0; I<minmn-NB; I+=NB ) {
                k1 = 1+I+NB;
                k2 = minmn;
                incx = 1;
                lapackf77_zlaswp(&NB, &a[I*lda], &lda, &k1, &k2, ipiv, &incx);
        }

    return MAGMA_SUCCESS;
} /* magma_zgetrf_piv */

#undef inAT
