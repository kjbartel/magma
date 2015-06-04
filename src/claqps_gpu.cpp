/*
    -- MAGMA (version 1.4.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       December 2013

       @generated c Tue Dec 17 13:18:36 2013

*/
#include "common_magma.h"
#include <cblas.h>

#define PRECISION_c

extern "C" magma_int_t
magma_claqps_gpu(magma_int_t m, magma_int_t n, magma_int_t offset,
             magma_int_t nb, magma_int_t *kb,
             magmaFloatComplex *A,  magma_int_t lda,
             magma_int_t *jpvt, magmaFloatComplex *tau,
             float *vn1, float *vn2,
             magmaFloatComplex *auxv,
             magmaFloatComplex *F,  magma_int_t ldf)
{
/*  -- MAGMA (version 1.4.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       December 2013

    Purpose
    =======
    CLAQPS computes a step of QR factorization with column pivoting
    of a complex M-by-N matrix A by using Blas-3.  It tries to factorize
    NB columns from A starting from the row OFFSET+1, and updates all
    of the matrix with Blas-3 xGEMM.

    In some cases, due to catastrophic cancellations, it cannot
    factorize NB columns.  Hence, the actual number of factorized
    columns is returned in KB.

    Block A(1:OFFSET,1:N) is accordingly pivoted, but not factorized.

    Arguments
    =========
    M       (input) INTEGER
            The number of rows of the matrix A. M >= 0.

    N       (input) INTEGER
            The number of columns of the matrix A. N >= 0

    OFFSET  (input) INTEGER
            The number of rows of A that have been factorized in
            previous steps.

    NB      (input) INTEGER
            The number of columns to factorize.

    KB      (output) INTEGER
            The number of columns actually factorized.

    A       (input/output) COMPLEX array, dimension (LDA,N)
            On entry, the M-by-N matrix A.
            On exit, block A(OFFSET+1:M,1:KB) is the triangular
            factor obtained and block A(1:OFFSET,1:N) has been
            accordingly pivoted, but no factorized.
            The rest of the matrix, block A(OFFSET+1:M,KB+1:N) has
            been updated.

    LDA     (input) INTEGER
            The leading dimension of the array A. LDA >= max(1,M).

    JPVT    (input/output) INTEGER array, dimension (N)
            JPVT(I) = K <==> Column K of the full matrix A has been
            permuted into position I in AP.

    TAU     (output) COMPLEX array, dimension (KB)
            The scalar factors of the elementary reflectors.

    VN1     (input/output) REAL array, dimension (N)
            The vector with the partial column norms.

    VN2     (input/output) REAL array, dimension (N)
            The vector with the exact column norms.

    AUXV    (input/output) COMPLEX array, dimension (NB)
            Auxiliar vector.

    F       (input/output) COMPLEX array, dimension (LDF,NB)
            Matrix F' = L*Y'*A.

    LDF     (input) INTEGER
            The leading dimension of the array F. LDF >= max(1,N).

    =====================================================================    */
    
#define  A(i, j) (A  + (i) + (j)*(lda ))
#define  F(i, j) (F  + (i) + (j)*(ldf ))

    magmaFloatComplex c_zero    = MAGMA_C_MAKE( 0.,0.);
    magmaFloatComplex c_one     = MAGMA_C_MAKE( 1.,0.);
    magmaFloatComplex c_neg_one = MAGMA_C_MAKE(-1.,0.);
    magma_int_t ione = 1;
    
    magma_int_t i__1, i__2;
    //float d__1;
    magmaFloatComplex z__1;
    
    //magma_int_t j;
    magma_int_t k, rk;
    //magmaFloatComplex Akk;
    magmaFloatComplex *Aks;
    magmaFloatComplex tauk = MAGMA_C_ZERO;
    magma_int_t pvt;
    //float temp, temp2;
    float tol3z;
    magma_int_t itemp;

    float lsticc, *lsticcs;
    magma_int_t lastrk;
    magma_smalloc( &lsticcs, 1+256*(n+255)/256 );

    lastrk = min( m, n + offset );
    tol3z = magma_ssqrt( lapackf77_slamch("Epsilon"));

    lsticc = 0;
    k = 0;
    magma_cmalloc( &Aks, nb );

    while( k < nb && lsticc == 0 ) {
        rk = offset + k;
        
        /* Determine ith pivot column and swap if necessary */
        // Fortran: pvt, k, isamax are all 1-based; subtract 1 from k.
        // C:       pvt, k, isamax are all 0-based; don't subtract 1.
        pvt = k - 1 + magma_isamax( n-k, &vn1[k], ione );
        
        if (pvt != k) {

            /*if (pvt >= nb) {
                // 1. Start copy from GPU
                magma_cgetmatrix_async( m - offset - nb, 1,
                                        dA(offset + nb, pvt), ldda,
                                        A (offset + nb, pvt), lda, stream );
            }*/

            /* F gets swapped so F must be sent at the end to GPU   */
            i__1 = k;
            /*if (pvt < nb){
                // no need of transfer if pivot is within the panel
                blasf77_cswap( &m, A(0, pvt), &ione, A(0, k), &ione );
            }
            else {
                // 1. Finish copy from GPU
                magma_queue_sync( stream );

                // 2. Swap as usual on CPU
                blasf77_cswap(&m, A(0, pvt), &ione, A(0, k), &ione);

                // 3. Restore the GPU
                magma_csetmatrix_async( m - offset - nb, 1,
                                        A (offset + nb, pvt), lda,
                                        dA(offset + nb, pvt), ldda, stream);
            }*/
            magmablas_cswap( m, A(0, pvt), ione, A(0, k), ione );

            //blasf77_cswap( &i__1, F(pvt,0), &ldf, F(k,0), &ldf );
            magmablas_cswap( i__1, F(pvt, 0), ldf, F(k, 0), ldf);
            itemp     = jpvt[pvt];
            jpvt[pvt] = jpvt[k];
            jpvt[k]   = itemp;
            //vn1[pvt] = vn1[k];
            //vn2[pvt] = vn2[k];
            #if defined(PRECISION_d) || defined(PRECISION_z)
                //magma_dswap( 1, &vn1[pvt], 1, &vn1[k], 1 );
                //magma_dswap( 1, &vn2[pvt], 1, &vn2[k], 1 );
                magma_dswap( 2, &vn1[pvt], n+offset, &vn1[k], n+offset );
            #else
                //magma_sswap( 1, &vn1[pvt], 1, &vn1[k], 1 );
                //magma_sswap( 1, &vn2[pvt], 1, &vn2[k], 1 );
                magma_sswap(2, &vn1[pvt], n+offset, &vn1[k], n+offset);
            #endif

        }

        /* Apply previous Householder reflectors to column K:
           A(RK:M,K) := A(RK:M,K) - A(RK:M,1:K-1)*F(K,1:K-1)'.
           Optimization: multiply with beta=0; wait for vector and subtract */
        if (k > 0) {
            /*#if (defined(PRECISION_c) || defined(PRECISION_z))
            for (j = 0; j < k; ++j){
                *F(k,j) = MAGMA_C_CNJG( *F(k,j) );
            }
            #endif*/

//#define RIGHT_UPDATE
#ifdef RIGHT_UPDATE
            i__1 = m - offset - nb;
            i__2 = k;
            magma_cgemv( MagmaNoTrans, i__1, i__2,
                         c_neg_one, A(offset+nb, 0), lda,
                                    F(k,         0), ldf,
                         c_one,     A(offset+nb, k), ione );
#else
            i__1 = m - rk;
            i__2 = k;
            /*blasf77_cgemv( MagmaNoTransStr, &i__1, &i__2,
                           &c_neg_one, A(rk, 0), &lda,
                                       F(k,  0), &ldf,
                           &c_one,     A(rk, k), &ione );*/
            magma_cgemv( MagmaNoTrans, i__1, i__2,
                         c_neg_one, A(rk, 0), lda,
                                    F(k,  0), ldf,
                         c_one,     A(rk, k), ione );
#endif

            /*#if (defined(PRECISION_c) || defined(PRECISION_z))
            for (j = 0; j < k; ++j) {
                *F(k,j) = MAGMA_C_CNJG( *F(k,j) );
            }
            #endif*/
        }
        
        /*  Generate elementary reflector H(k). */
        magma_clarfg_gpu(m-rk, A(rk, k), A(rk + 1, k), &tau[k], &vn1[k], &Aks[k]);

        //Akk = *A(rk, k);
        //*A(rk, k) = c_one;
        //magma_cgetvector( 1, &Aks[k],  1, &Akk,     1 );

        /* needed to avoid the race condition */
        if (k == 0) magma_csetvector(  1,    &c_one,       1, A(rk, k), 1 );
        else        magma_ccopymatrix( 1, 1, A(offset, 0), 1, A(rk, k), 1 );

        /* Compute Kth column of F:
           Compute  F(K+1:N,K) := tau(K)*A(RK:M,K+1:N)'*A(RK:M,K) on the GPU */
        if (k < n-1 || k > 0) magma_cgetvector( 1, &tau[k], 1, &tauk, 1 );
        if (k < n-1) {
            i__1 = m - rk;
            i__2 = n - k - 1;

            /* Send the vector to the GPU */
            //magma_csetmatrix( i__1, 1, A(rk, k), lda, dA(rk,k), ldda );

            /* Multiply on GPU */
            // was CALL CGEMV( 'Conjugate transpose', M-RK+1, N-K,
            //                 TAU( K ), A( RK,  K+1 ), LDA,
            //                           A( RK,  K   ), 1,
            //                 CZERO,    F( K+1, K   ), 1 )
            //magma_cgetvector( 1, &tau[k], 1, &tauk, 1 );
            magma_cgemv( MagmaConjTrans, m-rk, n-k-1,
                         tauk,   A( rk,  k+1 ), lda,
                                 A( rk,  k   ), 1,
                         c_zero, F( k+1, k   ), 1 );
            //magma_cscal( m-rk, tau[k], F( k+1, k), 1 );
            //magma_int_t i__3 = nb-k-1;
            //magma_int_t i__4 = i__2 - i__3;
            //magma_int_t i__5 = nb-k;
            //magma_cgemv( MagmaConjTrans, i__1 - i__5, i__2 - i__3,
            //             tau[k], dA(rk +i__5, k+1+i__3), ldda,
            //                     dA(rk +i__5, k       ), ione,
            //             c_zero, dF(k+1+i__3, k       ), ione );
            
            //magma_cgetmatrix_async( i__2-i__3, 1,
            //                        dF(k + 1 +i__3, k), i__2,
            //                        F (k + 1 +i__3, k), i__2, stream );
            
            //blasf77_cgemv( MagmaConjTransStr, &i__1, &i__3,
            //               &tau[k], A(rk,  k+1), &lda,
            //                        A(rk,  k  ), &ione,
            //               &c_zero, F(k+1, k  ), &ione );
            
            //magma_queue_sync( stream );
            //blasf77_cgemv( MagmaConjTransStr, &i__5, &i__4,
            //               &tau[k], A(rk, k+1+i__3), &lda,
            //                        A(rk, k       ), &ione,
            //               &c_one,  F(k+1+i__3, k ), &ione );
        }
        
        /* Padding F(1:K,K) with zeros.
        for (j = 0; j <= k; ++j) {
            magma_csetvector( 1, &c_zero, 1, F(j, k), 1 );
        }*/
        
        /* Incremental updating of F:
           F(1:N,K) := F(1:N,K)                        - tau(K)*F(1:N,1:K-1)*A(RK:M,1:K-1)'*A(RK:M,K).
           F(1:N,K) := tau(K)*A(RK:M,K+1:N)'*A(RK:M,K) - tau(K)*F(1:N,1:K-1)*A(RK:M,1:K-1)'*A(RK:M,K)
                    := tau(K)(A(RK:M,K+1:N)' - F(1:N,1:K-1)*A(RK:M,1:K-1)') A(RK:M,K)
           so, F is (updated A)*V */
        //if (k > 0 && k<n-1) {
        if (k > 0) {
            //magma_cgetvector( 1, &tau[k], 1, &tauk, 1 );
            z__1 = MAGMA_C_NEGATE( tauk );
#ifdef RIGHT_UPDATE
            i__1 = m - offset - nb;
            i__2 = k;
            magma_cgemv( MagmaConjTrans, i__1, i__2,
                         z__1,   A(offset+nb, 0), lda,
                                 A(offset+nb, k), ione,
                         c_zero, auxv, ione );
            
            i__1 = k;
            magma_cgemv( MagmaNoTrans, n-k-1, i__1,
                         c_one, F(k+1,0), ldf,
                                auxv,     ione,
                         c_one, F(k+1,k), ione );
#else
            i__1 = m - rk;
            i__2 = k;
            //blasf77_cgemv( MagmaConjTransStr, &i__1, &i__2,
            //               &z__1,   A(rk, 0), &lda,
            //                        A(rk, k), &ione,
            //               &c_zero, auxv, &ione );

            magma_cgemv( MagmaConjTrans, i__1, i__2,
                         z__1,   A(rk, 0), lda,
                                 A(rk, k), ione,
                         c_zero, auxv, ione );
            
            //i__1 = k;
            //blasf77_cgemv( MagmaNoTransStr, &n, &i__1,
            //               &c_one, F(0,0), &ldf,
            //                       auxv,   &ione,
            //               &c_one, F(0,k), &ione );
            /*magma_cgemv( MagmaNoTrans, n, i__1,
                           c_one, F(0,0), ldf,
                                  auxv,   ione,
                           c_one, F(0,k), ione );*/
            /* I think we only need stricly lower-triangular part :) */
            magma_cgemv( MagmaNoTrans, n-k-1, i__2,
                         c_one, F(k+1,0), ldf,
                                auxv,     ione,
                         c_one, F(k+1,k), ione );
#endif
        }
        
        /* Optimization: On the last iteration start sending F back to the GPU */
        
        /* Update the current row of A:
           A(RK,K+1:N) := A(RK,K+1:N) - A(RK,1:K)*F(K+1:N,1:K)'.               */
        if (k < n-1) {
            i__1 = n - k - 1;
            i__2 = k + 1;
            //blasf77_cgemm( MagmaNoTransStr, MagmaConjTransStr, &ione, &i__1, &i__2,
            //               &c_neg_one, A(rk, 0  ), &lda,
            //                           F(k+1,0  ), &ldf,
            //               &c_one,     A(rk, k+1), &lda );
#ifdef RIGHT_UPDATE
            /* right-looking update of rows,                     */
            magma_cgemm( MagmaNoTrans, MagmaConjTrans, nb-k, i__1, ione,
                         c_neg_one, A(rk,  k  ), lda,
                                    F(k+1, k  ), ldf,
                         c_one,     A(rk,  k+1), lda );
#else
            /* left-looking update of rows,                     *
             * since F=A'v with original A, so no right-looking */
            magma_cgemm( MagmaNoTrans, MagmaConjTrans, ione, i__1, i__2,
                         c_neg_one, A(rk, 0  ), lda,
                                    F(k+1,0  ), ldf,
                         c_one,     A(rk, k+1), lda );
#endif
        }
        
        /* Update partial column norms. */
        if (rk < min(m, n+offset)-1 ){
            magmablas_scnrm2_row_check_adjust(n-k-1, tol3z, &vn1[k+1], &vn2[k+1], A(rk,k+1), lda, lsticcs);

            magma_device_sync();
            #if defined(PRECISION_d) || defined(PRECISION_z)
            magma_dgetvector( 1, &lsticcs[0], 1, &lsticc, 1 );
            #else
            magma_sgetvector( 1, &lsticcs[0], 1, &lsticc, 1 );
            #endif
        }


        /*if (rk < lastrk) {
            for (j = k + 1; j < n; ++j) {
                if (vn1[j] != 0.) {
                    // NOTE: The following 4 lines follow from the analysis in
                    //   Lapack Working Note 176.
                    temp = MAGMA_C_ABS( *A(rk,j) ) / vn1[j];
                    temp = max( 0., ((1. + temp) * (1. - temp)) );

                    d__1 = vn1[j] / vn2[j];
                    temp2 = temp * (d__1 * d__1);

                    if (temp2 <= tol3z) {
                        vn2[j] = (float) lsticc;
                        lsticc = j;
                    } else {
                        vn1[j] *= magma_ssqrt(temp);
                    }
                }
            }
        }*/
        
        //*A(rk, k) = Akk;
        //magma_csetvector( 1, &Akk, 1, A(rk, k), 1 );
        //magma_cswap( 1, &Aks[k], 1, A(rk, k), 1 );
        
        ++k;
    }
    magma_ccopymatrix( 1, k, Aks, 1, A(offset, 0), lda+1 );

    // leave k as the last column done
    --k;
    *kb = k + 1;
    rk = offset + *kb - 1;

    /* Apply the block reflector to the rest of the matrix:
       A(OFFSET+KB+1:M,KB+1:N) := A(OFFSET+KB+1:M,KB+1:N) - A(OFFSET+KB+1:M,1:KB)*F(KB+1:N,1:KB)'  */
    if (*kb < min(n, m - offset)) {
        i__1 = m - rk - 1;
        i__2 = n - *kb;
        
        /* Send F to the GPU
        magma_csetmatrix( i__2, *kb,
                          F (*kb, 0), ldf,
                          dF(*kb, 0), i__2 );*/

        magma_cgemm( MagmaNoTrans, MagmaConjTrans, i__1, i__2, *kb,
                     c_neg_one, A(rk+1, 0  ), lda,
                                F(*kb,  0  ), ldf,
                     c_one,     A(rk+1, *kb), lda );
    }
    /* Recomputation of difficult columns. */
    if( lsticc > 0 ) {
        printf( " -- recompute dnorms --\n" );
        magmablas_scnrm2_check(m-rk-1, n-*kb, A(rk+1,*kb), lda,
                               &vn1[*kb], lsticcs);
#if defined(PRECISION_d) || defined(PRECISION_z)
        magma_dcopymatrix( n-*kb, 1, &vn1[*kb], *kb, &vn2[*kb], *kb);
#else
        magma_scopymatrix( n-*kb, 1, &vn1[*kb], *kb, &vn2[*kb], *kb);
#endif
    /*while( lsticc > 0 ) {
        itemp = (magma_int_t)(vn2[lsticc] >= 0. ? floor(vn2[lsticc] + .5) : -floor(.5 - vn2[lsticc]));
        i__1 = m - rk - 1;
        if (lsticc <= nb)
            vn1[lsticc] = cblas_scnrm2(i__1, A(rk + 1, lsticc), ione);
        else {
            // Where is the data, CPU or GPU ?
            float r1, r2;
            
            r1 = cblas_scnrm2(nb-k, A(rk + 1, lsticc), ione);
            r2 = magma_scnrm2(m-offset-nb, dA(offset + nb + 1, lsticc), ione);
            
            vn1[lsticc] = magma_ssqrt(r1*r1+r2*r2);
        }
        
        // NOTE: The computation of VN1( LSTICC ) relies on the fact that
        //   SNRM2 does not fail on vectors with norm below the value of SQRT(SLAMCH('S'))
        vn2[lsticc] = vn1[lsticc];
        lsticc = itemp;*/
    }
    magma_free(Aks);
    magma_free(lsticcs);

    return MAGMA_SUCCESS;
} /* magma_claqps */
