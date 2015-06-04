/*
  -- MAGMA (version 1.2.1) --
     Univ. of Tennessee, Knoxville
     Univ. of California, Berkeley
     Univ. of Colorado, Denver
     June 2012

     @generated s Thu Jun 28 12:30:47 2012

*/

#include "common_magma.h"
#include <cblas.h>


extern "C" magma_int_t 
magma_sgeqp3(magma_int_t m_, magma_int_t n_, float *a, 
             magma_int_t lda_, magma_int_t *jpvt, float *tau, float *work, 
             magma_int_t lwork_, magma_int_t *info)
{
/*  -- MAGMA (version 1.2.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       June 2012
       
    Purpose   
    =======   
    SGEQP3 computes a QR factorization with column pivoting of a   
    matrix A:  A*P = Q*R  using Level 3 BLAS.   

    Arguments   
    =========   
    M       (input) INTEGER
            The number of rows of the matrix A. M >= 0.   

    N       (input) INTEGER   
            The number of columns of the matrix A.  N >= 0.   

    A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)   
            On entry, the M-by-N matrix A.   
            On exit, the upper triangle of the array contains the   
            min(M,N)-by-N upper trapezoidal matrix R; the elements below   
            the diagonal, together with the array TAU, represent the   
            orthogonal matrix Q as a product of min(M,N) elementary   
            reflectors.   

    LDA     (input) INTEGER   
            The leading dimension of the array A. LDA >= max(1,M).   

    JPVT    (input/output) INTEGER array, dimension (N)   
            On entry, if JPVT(J).ne.0, the J-th column of A is permuted   
            to the front of A*P (a leading column); if JPVT(J)=0,   
            the J-th column of A is a free column.   
            On exit, if JPVT(J)=K, then the J-th column of A*P was the   
            the K-th column of A.   

    TAU     (output) DOUBLE PRECISION array, dimension (min(M,N))   
            The scalar factors of the elementary reflectors.   

    WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))   
            On exit, if INFO=0, WORK(1) returns the optimal LWORK.   

    LWORK   (input) INTEGER   
            The dimension of the array WORK, LWORK >= 2*N + ( N+1 )*NB, 
            where NB is the optimal blocksize.   

            If LWORK = -1, then a workspace query is assumed; the routine   
            only calculates the optimal size of the WORK array, returns   
            this value as the first entry of the WORK array, and no error   
            message related to LWORK is issued by XERBLA.   

    INFO    (output) INTEGER   
            = 0: successful exit.   
            < 0: if INFO = -i, the i-th argument had an illegal value.   

    Further Details   
    ===============   
    The matrix Q is represented as a product of elementary reflectors   

      Q = H(1) H(2) . . . H(k), where k = min(m,n).   

    Each H(i) has the form   

      H(i) = I - tau * v * v'   

    where tau is a real scalar, and v is a real vector   
    with v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in   
    A(i+1:m,i), and tau in TAU(i).   
    =====================================================================   */

#define  A(i, j) (a    +(j)*(*lda) + (i))
#define dA(i, j) (dwork+(j)* ldda  + (i))

    magma_int_t *m     = &m_;
    magma_int_t *n     = &n_;
    magma_int_t *lda   = &lda_;
    magma_int_t *lwork = &lwork_;

    float   *dwork, *df;

    magma_int_t c__1 = 1;
    
    magma_int_t i__1, i__2, ldda;
    magma_int_t j, jb, na, nb, sm, sn, nx, fjb, iws, nfxd, nbmin, minmn, minws;

    magma_int_t topbmn, sminmn, lwkopt, lquery;
    
    a -= 1 + *lda;
    --jpvt;
    --tau;
    --work;
    
    *info = 0;
    lquery = *lwork == -1;
    if (*m < 0) {
       *info = -1;
    } else if (*n < 0) {
       *info = -2;
    } else if (*lda < max(1,*m)) {
       *info = -4;
    }
    
    if (*info == 0) {
        minmn = min(*m,*n);
        if (minmn == 0) {
            iws = 1;
            lwkopt = 1;
        } else {
            iws = *n + 1;
            nb = magma_get_sgeqp3_nb(min(*m, *n));
            lwkopt = 2*(*n) + (*n + 1) * nb;
        }
        MAGMA_S_SET2REAL(work[1],(float)lwkopt);

        if (*lwork < iws && ! lquery) {
            *info = -8;
        }
    }

    if (*info != 0) {
        i__1 = -(*info);
        magma_xerbla( __func__, -(*info) );
        return *info;
    } else if (lquery) {
        return *info;
    }

    if (minmn == 0)
        return *info;

    ldda = ((*m+31)/32)*32;
    if (MAGMA_SUCCESS != magma_smalloc( &dwork, (*n)*ldda +2*(*n)+(*n+1)*nb)){
        *info = MAGMA_ERR_DEVICE_ALLOC;
        return *info;
    }
    df     = dwork + (*n)*ldda + 2*(*n); 
    dwork -= 1 + ldda;
    

    /* Move initial columns up front. */
    nfxd = 1;
    for (j = 1; j <= *n; ++j) {
        if (jpvt[j] != 0) {
           if (j != nfxd) {
               blasf77_sswap(m, A(1, j), &c__1, A(1, nfxd), &c__1);
               jpvt[j] = jpvt[nfxd];
               jpvt[nfxd] = j;
           } else {
               jpvt[j] = j;
           }
           ++nfxd;
        } else
            jpvt[j] = j;
    }
    --nfxd;

    /*     Factorize fixed columns   
           =======================   
           Compute the QR factorization of fixed columns and update   
           remaining columns. */
    if (nfxd > 0) {
        na = min(*m,nfxd);
        lapackf77_sgeqrf(m, &na, A(1, 1), lda, &tau[1], &work[1], lwork, info);

        iws = max(iws, (magma_int_t)work[1]);
        if (na < *n) {
            i__1 = *n - na;
            lapackf77_sormqr(MagmaLeftStr, MagmaTransStr, m, &i__1, &na, A(1,1),
                             lda, &tau[1], A(1, na + 1), lda, &work[1], 
                             lwork, info);
            iws = max(iws, (magma_int_t) work[1]);
        }
    }
    
    /*  Factorize free columns */
    if (nfxd < minmn) {
       sm = *m - nfxd;
       sn = *n - nfxd;
       sminmn = minmn - nfxd;

       /* Determine the block size. */
       nbmin = 2;
       nx = 0;

       if (nb > 1 && nb < sminmn) {
           /* Determine when to cross over from blocked to unblocked code */
           nx = nb;

           if (nx < sminmn) {               
               /* Determine if workspace is large enough for blocked code. */
               minws = 2*sn + (sn + 1) * nb;
               iws = max(iws,minws);
               if (*lwork < minws) {
       
                   /* Not enough workspace to use optimal NB: Reduce NB and   
                      determine the minimum value of NB. */
                   nb = (*lwork - 2*sn)/ (sn + 1);
                   nbmin = max(i__1,i__2);
               }
           }
       }

       /*       Initialize partial column norms. The first N elements of work   
                store the exact column norms. */
       for (j = nfxd + 1; j <= *n; ++j) {
           work[j] = cblas_snrm2(sm, A(nfxd + 1, j), c__1);
           work[*n + j] = work[j];
       }
       
       if (nb >= nbmin && nb < sminmn && nx < sminmn) {
           /* Use blocked code initially. */   
           j = nfxd + 1;
           
           // Set the original matrix to the GPU
           magma_ssetmatrix( *m, sn,
                             A(1,j),  *lda,
                             dA(1,j), ldda );
           
           /* Compute factorization: while loop. */
           topbmn = minmn - nx;
       L30:
           if (j <= topbmn) {
               jb = min(nb, topbmn - j + 1);
              
               /* Factorize JB columns among columns J:N. */               
               i__1 = *n - j + 1;
               i__2 = j - 1;
               
               // Get panel to the CPU
               magma_sgetmatrix( *m-j+1, jb,
                                 dA(j,j), ldda,
                                 A (j,j), *lda );
 
               // Get the rows 
               magma_sgetmatrix( jb, i__1 - jb, 
                                 dA(j,j + jb), ldda,
                                 A (j,j + jb), *lda );
               
               magma_slaqps(m, &i__1, &i__2, &jb, &fjb, 
                            A(1, j), lda,
                            dA(1, j), &ldda,  
                            &jpvt[j], &tau[j], &work[j], &work[*n + j], 
                            &work[2*(*n)+1], &work[2*(*n)+1+jb], &i__1,
                            &df[jb], &i__1);

               j += fjb;
               goto L30;
           }
       } else {
           j = nfxd + 1;
       }
        
       /* Use unblocked code to factor the last or only block. */
       if (j <= minmn) {
           i__1 = *n - j + 1;
           i__2 = j - 1;
           
           if (j > nfxd + 1)
               magma_sgetmatrix( *m-j+1, i__1,
                                 dA(j,j), ldda,
                                 A (j,j), *lda );
           
           lapackf77_slaqp2(m, &i__1, &i__2, A(1, j), lda, &jpvt[j], 
                            &tau[j], &work[j], &work[*n + j], &work[2*(*n) + 1]);
       }
    }

    work[1] = MAGMA_S_MAKE((float) iws, 0.);
    magma_free( dA(1,1) );

    return *info;
} /* sgeqp3_ */
