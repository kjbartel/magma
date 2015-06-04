/*
    -- MAGMA (version 1.2.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       June 2012

       @author Azzam Haidar
       @author Stan Tomov

       @generated c Thu Jun 28 12:30:59 2012

*/
#include "common_magma.h"

// === Define what BLAS to use ============================================
#define PRECISION_c

#if (defined(PRECISION_s))
     #define magma_ssyr2k magmablas_ssyr2k
#endif
// === End defining what BLAS to use ======================================


extern "C" magma_int_t
magma_chebbd(char uplo, magma_int_t n, magma_int_t nb,
              cuFloatComplex *a, magma_int_t lda, 
              cuFloatComplex *tau,
              cuFloatComplex *work, magma_int_t lwork,
              cuFloatComplex *dT,
              magma_int_t *info)
{
/*  -- MAGMA (version 1.2.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       June 2012

    Purpose   
    =======   
    CHEBBD reduces a complex Hermitian matrix A to real symmetric   
    band-diagonal form T by an orthogonal similarity transformation:   
    Q**H * A * Q = T.   
    This version stores the triangular matrices T used in the accumulated
    Householder transformations (I - V T V').

    Arguments   
    =========   
    UPLO    (input) CHARACTER*1   
            = 'U':  Upper triangle of A is stored;   
            = 'L':  Lower triangle of A is stored.   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    A       (input/output) COMPLEX array, dimension (LDA,N)   
            On entry, the Hermitian matrix A.  If UPLO = 'U', the leading   
            N-by-N upper triangular part of A contains the upper   
            triangular part of the matrix A, and the strictly lower   
            triangular part of A is not referenced.  If UPLO = 'L', the   
            leading N-by-N lower triangular part of A contains the lower   
            triangular part of the matrix A, and the strictly upper   
            triangular part of A is not referenced.   
            On exit, if UPLO = 'U', the Upper band-diagonal of A is 
            overwritten by the corresponding elements of the   
            band-diagonal matrix T, and the elements above the band   
            diagonal, with the array TAU, represent the orthogonal   
            matrix Q as a product of elementary reflectors; if UPLO   
            = 'L', the the Lower band-diagonal of A is overwritten by 
            the corresponding elements of the band-diagonal   
            matrix T, and the elements below the band-diagonal, with   
            the array TAU, represent the orthogonal matrix Q as a product   
            of elementary reflectors. See Further Details.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= max(1,N).   

    TAU     (output) COMPLEX array, dimension (N-1)   
            The scalar factors of the elementary reflectors (see Further   
            Details).   

    WORK    (workspace/output) COMPLEX array, dimension (MAX(1,LWORK))   
            On exit, if INFO = 0, WORK(1) returns the optimal LWORK.   

    LWORK   (input) INTEGER   
            The dimension of the array WORK.  LWORK >= 1.   
            For optimum performance LWORK >= N*NB, where NB is the   
            optimal blocksize.   

            If LWORK = -1, then a workspace query is assumed; the routine   
            only calculates the optimal size of the WORK array, returns   
            this value as the first entry of the WORK array, and no error   
            message related to LWORK is issued by XERBLA.   

    dT      (output) COMPLEX array on the GPU, dimension N*NB, 
            where NB is the optimal blocksize.
            On exit dT holds the upper triangular matrices T from the 
            accumulated Householder transformations (I - V T V') used
            in the factorization. The nb x nb matrices T are ordered 
            consecutively in memory one after another.

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   

    Further Details   
    ===============   
    If UPLO = 'U', the matrix Q is represented as a product of elementary   
    reflectors   

       Q = H(n-1) . . . H(2) H(1).   

    Each H(i) has the form   

       H(i) = I - tau * v * v'

    where tau is a complex scalar, and v is a complex vector with   
    v(i+1:n) = 0 and v(i) = 1; v(1:i-1) is stored on exit in   
    A(1:i-1,i+1), and tau in TAU(i).   

    If UPLO = 'L', the matrix Q is represented as a product of elementary   
    reflectors   

       Q = H(1) H(2) . . . H(n-1).   

    Each H(i) has the form   

       H(i) = I - tau * v * v'   

    where tau is a complex scalar, and v is a complex vector with   
    v(1:i) = 0 and v(i+1) = 1; v(i+2:n) is stored on exit in A(i+2:n,i),   
    and tau in TAU(i).

    The contents of A on exit are illustrated by the following examples   
    with n = 5:   

    if UPLO = 'U':                       if UPLO = 'L':   

      (  d   e   v2  v3  v4 )              (  d                  )   
      (      d   e   v3  v4 )              (  e   d              )   
      (          d   e   v4 )              (  v1  e   d          )   
      (              d   e  )              (  v1  v2  e   d      )   
      (                  d  )              (  v1  v2  v3  e   d  )   

    where d and e denote diagonal and off-diagonal elements of T, and vi   
    denotes an element of the vector defining H(i).   
    =====================================================================    */

    #define a_ref(a_1,a_2) ( a+((a_2)-1)*( lda) + (a_1)-1)
    #define da_ref(a_1,a_2) (da+((a_2)-1)*(ldda) + (a_1)-1)
    #define tau_ref(a_1)    (tau + (a_1)-1)
    #define t_ref(a_1)      (dT  + ((a_1)-1)*(lddt))

    char uplo_[2] = {uplo, 0};

    int ldda = ((n+31)/32)*32; // magma_get_chebbd_nb(n); 
    int lddt = nb;
   
    cuFloatComplex c_neg_one  = MAGMA_C_NEG_ONE, c_neg_half;
    MAGMA_C_SSCALE(c_neg_half, c_neg_one, 2.);
    cuFloatComplex c_one  = MAGMA_C_ONE ;
    cuFloatComplex c_zero = MAGMA_C_ZERO;
    float  d_one = MAGMA_D_ONE;

    magma_int_t pm, pn, indi, indj, pk;
    magma_int_t pm_old=0, pn_old=0, indi_old=0, indj_old=0;

    int i;
    int lwkopt;
    int lquery;

    *info = 0;
    int upper = lapackf77_lsame(uplo_, "U");
    lquery = lwork == -1;
    if (! upper && ! lapackf77_lsame(uplo_, "L")) {
        *info = -1;
    } else if (n < 0) {
        *info = -2;
    } else if (lda < max(1,n)) {
        *info = -4;
    } else if (lwork < 1 && ! lquery) {
        *info = -9;
    }

    if (*info == 0) {
      /* Determine the block size. */
      lwkopt = n * nb;
      MAGMA_C_SET2REAL( work[0], lwkopt );
    }

    if (*info != 0)
      return *info;
    else if (lquery)
      return *info;

    /* Quick return if possible */
    if (n == 0) {
        work[0] = c_one;
        return *info;
    }

    cuFloatComplex *da;
    if (MAGMA_SUCCESS != magma_cmalloc( &da, (n + 2*nb)*ldda )) {
        *info = MAGMA_ERR_DEVICE_ALLOC;
        return *info;
    }

    /* Use the first panel of da as work space */
    cuFloatComplex *dwork = da+n*ldda;
    cuFloatComplex *dW    = dwork + nb*ldda;

    cudaStream_t stream[2];
    magma_queue_create( &stream[0] );
    magma_queue_create( &stream[1] );

    cuFloatComplex *hT = work + lwork - nb*nb;
    lwork -= nb*nb;
    memset( hT, 0, nb*nb*sizeof(cuFloatComplex));


    if (upper) {

      printf("CHEBBD is not yet implemented for upper matrix storage. Exit.\n");
      exit(1);

    }else {
        /* Copy the matrix to the GPU */
        if (1 <= n-nb){
            magma_csetmatrix_async( (n-nb), (n-nb),
                                    a_ref(nb+1, nb+1),  lda,
                                    da_ref(nb+1, nb+1), ldda, stream[0] );
        }

        /* Reduce the lower triangle of A */
        for (i = 1; i <= n-nb; i += nb) 
        {
             indi = i+nb;
             indj = i;
             pm   = n - i - nb + 1;
             //pn   = min(i+nb-1, n-nb) -i + 1;
             pn   = nb;
             
             /*   Get the current panel (no need for the 1st iteration) */
             if (i > 1 ){
                 // cpanel_to_q copy the upper oof diagonal part of 
                 // the matrix to work to be restored later. acctually
                 //  the zero's and one's putted are not used this is only
                 //   because we don't have a function that copy only the
                 //    upper part of A to be restored after copying the 
                 //    lookahead panel that has been computted from GPU to CPU. 
                 cpanel_to_q(MagmaUpper, pn-1, a_ref(i, i+1), lda, work);
                 magma_cgetmatrix_async( (pm+pn), pn,
                                         da_ref( i, i), ldda,
                                         a_ref ( i, i), lda, stream[1] );

                 magma_cher2k(MagmaLower, MagmaNoTrans, pm_old-pn_old, pn_old, c_neg_one,
                      da_ref(indi_old+pn_old, indj_old), ldda,
                      dW + pn_old           , pm_old, d_one,
                      da_ref(indi_old+pn_old, indi_old+pn_old), ldda);

                 magma_queue_sync( stream[1] );
                 cq_to_panel(MagmaUpper, pn-1, a_ref(i, i+1), lda, work);
             }

             /* ==========================================================
                QR factorization on a panel starting nb off of the diagonal.
                Prepare the V and T matrices. 
                ==========================================================  */
             lapackf77_cgeqrf(&pm, &pn, a_ref(indi, indj), &lda, 
                        tau_ref(i), work, &lwork, info);
             
             /* Form the matrix T */
                         pk=min(pm,pn);
             lapackf77_clarft( MagmaForwardStr, MagmaColumnwiseStr,
                           &pm, &pk, a_ref(indi, indj), &lda,
                           tau_ref(i), hT, &nb);

             /* Prepare V - put 0s in the upper triangular part of the panel
                (and 1s on the diagonal), temporaly storing the original in work */
             cpanel_to_q(MagmaUpper, pk, a_ref(indi, indj), lda, work);

             /* Send V from the CPU to the GPU */
             magma_csetmatrix_async( pm, pk,
                                     a_ref(indi, indj),  lda,
                                     da_ref(indi, indj), ldda, stream[0] );

             /* Send the triangular factor T to the GPU */
             magma_csetmatrix_async( pk, pk,
                                     hT,       nb,
                                     t_ref(i), lddt, stream[0] );
             /* ==========================================================
                Compute W:
                1. X = A (V T)
                2. W = X - 0.5* V * (T' * (V' * X)) 
                ==========================================================  */
             /* dwork = V T */
             magma_queue_sync( stream[0] );
             magma_cgemm(MagmaNoTrans, MagmaNoTrans, pm, pk, pk,
                         c_one, da_ref(indi, indj), ldda, 
                         t_ref(i), lddt,
                         c_zero, dwork, pm);
             /* W = X = A*V*T = A dwork */ 
             magma_chemm(MagmaLeft, uplo, pm, pk,
                         c_one, da_ref(indi, indi), ldda,
                         dwork, pm,
                         c_zero, dW, pm);
             /* restore the panel */
             cq_to_panel(MagmaUpper, pk, a_ref(indi, indj), lda, work);
             
             /* dwork = V*T already ==> dwork' =V'*T'
              * compute T'*V'*X ==> dwork'*W ==>
              * dwork + pm*nb = ((T' * V') * X) = dwork' * X = dwork' * W */
             magma_cgemm(MagmaConjTrans, MagmaNoTrans, pk, pk, pm,
                         c_one, dwork, pm, 
                         dW, pm,
                         c_zero, dwork + pm*nb, nb);
             /* W = X - 0.5 * V * T'*V'*X
              *   = X - 0.5 * V * (dwork + pm*nb) = W - 0.5 * V * (dwork + pm*nb) */
             magma_cgemm(MagmaNoTrans, MagmaNoTrans, pm, pk, pk,
                         c_neg_half, da_ref(indi, indj), ldda,
                         dwork + pm*nb, nb, 
                         c_one,     dW, pm);
             
             /* ==========================================================
                Update the unreduced submatrix A(i+ib:n,i+ib:n), using   
                an update of the form:  A := A - V*W' - W*V' 
                ==========================================================  */
             if (i + nb <= n-nb){
                 /* There would be next iteration;
                    do lookahead - update the next panel */
                 magma_cgemm(MagmaNoTrans, MagmaConjTrans, pm, pn, pn, c_neg_one,
                             da_ref(indi, indj), ldda,
                             dW                , pm, c_one,
                             da_ref(indi, indi), ldda);
                 magma_cgemm(MagmaNoTrans, MagmaConjTrans, pm, pn, pn, c_neg_one,
                             dW                , pm,
                             da_ref(indi, indj), ldda, c_one,
                             da_ref(indi, indi), ldda);
             }
             else {
                 /* no look-ahead as this is last iteration */
                 magma_cher2k(MagmaLower, MagmaNoTrans, pk, pk, c_neg_one,
                              da_ref(indi, indj), ldda,
                              dW                , pm, d_one,
                              da_ref(indi, indi), ldda);
             }
             
             indi_old = indi;
             indj_old = indj;
             pm_old   = pm;
             pn_old   = pn;
        }  // end loop for(i)

        /* Send the last block to the CPU */
        pk = min(pm,pn);
        if (1 <= n-nb){
            cpanel_to_q(MagmaUpper, pk-1, a_ref(n-pk+1, n-pk+2), lda, work);
            magma_cgetmatrix( pk, pk,
                              da_ref(n-pk+1, n-pk+1), ldda,
                              a_ref(n-pk+1, n-pk+1),  lda );
            cq_to_panel(MagmaUpper, pk-1, a_ref(n-pk+1, n-pk+2), lda, work);
        }
    }// end of LOWER

    magma_queue_destroy( stream[0] );
    magma_queue_destroy( stream[1] );
    magma_free( da );
    MAGMA_C_SET2REAL( work[0], lwkopt );
    return *info;
} /* chebbd_ */