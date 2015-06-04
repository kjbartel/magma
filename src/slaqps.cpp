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

#define PRECISION_s

extern "C" magma_int_t 
magma_slaqps(magma_int_t *m, magma_int_t *n, magma_int_t *offset, 
             magma_int_t *nb, magma_int_t *kb, 
             float *a, magma_int_t *lda,
             float *da, magma_int_t *ldda,
             magma_int_t *jpvt, float *tau, float *vn1, float *vn2, 
             float *auxv, 
             float *f, magma_int_t *ldf,
             float *df, magma_int_t *lddf)
{
/*
    -- MAGMA (version 1.2.1) --
    Univ. of Tennessee, Knoxville
    Univ. of California, Berkeley
    Univ. of Colorado, Denver
    June 2012

    Purpose   
    =======   
    SLAQPS computes a step of QR factorization with column pivoting   
    of a real M-by-N matrix A by using Blas-3.  It tries to factorize   
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

    A       (input/output) COMPLEX*16 array, dimension (LDA,N)   
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

    TAU     (output) COMPLEX*16 array, dimension (KB)   
            The scalar factors of the elementary reflectors.   

    VN1     (input/output) DOUBLE PRECISION array, dimension (N)   
            The vector with the partial column norms.   

    VN2     (input/output) DOUBLE PRECISION array, dimension (N)   
            The vector with the exact column norms.   

    AUXV    (input/output) COMPLEX*16 array, dimension (NB)   
            Auxiliar vector.   

    F       (input/output) COMPLEX*16 array, dimension (LDF,NB)   
            Matrix F' = L*Y'*A.   

    LDF     (input) INTEGER   
            The leading dimension of the array F. LDF >= max(1,N).   

    =====================================================================    */
    
#define  A(i, j) (a    +(j)*(*lda)  + (i))
#define dA(i, j) (da   +(j)*(*ldda) + (i))

    float zero = MAGMA_S_MAKE(0.,0.);
    float one  = MAGMA_S_MAKE(1.,0.);
    float mone = MAGMA_S_MAKE(-1.,0.);
    magma_int_t c__1 = 1;
    
    magma_int_t f_dim1, i__1, i__2;
    float d__1;
    float z__1;
    
    magma_int_t j, k, rk;
    float akk;
    magma_int_t pvt;
    float temp, temp2, tol3z;
    magma_int_t itemp;

    magma_int_t lsticc;
    magma_int_t lastrk;    
    
    a -= 1 + *lda;
    --jpvt;
    --tau;
    --vn1;
    --vn2;
    --auxv;
    f_dim1 = *ldf;
    f -= 1 + f_dim1;

    i__1 = *m, i__2 = *n + *offset;
    lastrk = min(i__1,i__2);
    lsticc = 0;
    k = 0;
    tol3z = magma_ssqrt(slamch_("Epsilon"));

    da -= 1 + *ldda;
    df -= 1 + *lddf;

    cudaStream_t stream;
    magma_queue_create( &stream );

 L10:
    if (k < *nb && lsticc == 0) {
        ++k;
        rk = *offset + k;
        
        /* Determine ith pivot column and swap if necessary */
        i__1 = *n - k + 1;
        
        /* Comment:
           Fortran BLAS does not have to add 1
           C       BLAS must add one to cblas_isamax */
        //pvt = k - 1 + isamax_(&i__1, &vn1[k], &c__1);
        pvt = k - 1 + cblas_isamax(i__1, &vn1[k], c__1) +1;
        
        if (pvt != k) {
            if (pvt <= *nb){
                /* no need of transfer if pivot is within the panel */ 
                blasf77_sswap(m, A(1, pvt), &c__1, A(1, k), &c__1);
            }
            else {
                /* 1. Start copy from GPU                           */
                magma_sgetmatrix( *m - *offset - *nb, 1,
                                  dA(*offset + *nb+1,pvt), *ldda,
                                  A (*offset + *nb+1,pvt), *lda );
                
                /* 2. Swap as usual on CPU                          */
                blasf77_sswap(m, A(1, pvt), &c__1, A(1, k), &c__1);

                /* 3. Restore the GPU                               */
                magma_ssetmatrix_async( *m - *offset - *nb, 1,
                                        A (*offset + *nb+1,pvt), *lda,
                                        dA(*offset + *nb+1,pvt), *ldda, stream);
            }
            
            i__1 = k - 1;

            /* F gets swapped so F must be sent at the end to GPU   */
            blasf77_sswap(&i__1, &f[pvt + f_dim1], ldf, &f[k + f_dim1], ldf);
            itemp = jpvt[pvt];
            jpvt[pvt] = jpvt[k];
            jpvt[k] = itemp;
            vn1[pvt] = vn1[k];
            vn2[pvt] = vn2[k];
        }

        /* Apply previous Householder reflectors to column K:   
           A(RK:M,K) := A(RK:M,K) - A(RK:M,1:K-1)*F(K,1:K-1)'. 
           Optimization: multiply with beta=0; wait for vector and subtract */
        if (k > 1) {
            #if (defined(PRECISION_c) || defined(PRECISION_z))
            for (j = 1; j < k; ++j){
                f[ k + j * f_dim1 ] = MAGMA_S_CNJG(f[k + j * f_dim1]);
            }
            #endif

            i__1 = *m - rk + 1;
            i__2 = k - 1;
            blasf77_sgemv(MagmaNoTransStr, &i__1, &i__2, &mone, A(rk,1), lda, 
                          &f[k + f_dim1], ldf, &one, A(rk, k), &c__1);

            #if (defined(PRECISION_c) || defined(PRECISION_z))
            for (j = 1; j < k; ++j) {
                f[k + j * f_dim1] =MAGMA_S_CNJG( f[k + j * f_dim1] );
            }
            #endif
        }
        
        /*  Generate elementary reflector H(k). */
        if (rk < *m) {
            i__1 = *m - rk + 1;
            lapackf77_slarfg(&i__1, A(rk, k), A(rk + 1, k), &c__1, &tau[k]);
        } else {
            lapackf77_slarfg(&c__1, A(rk, k), A(rk, k), &c__1, &tau[k]);
        }
        
        akk = *A(rk, k);
        *A(rk, k) = one;

       /* Compute Kth column of F:   
          Compute  F(K+1:N,K) := tau(K)*A(RK:M,K+1:N)'*A(RK:M,K) on the GPU */
       if (k < *n) {
           i__1 = *m - rk + 1;
           i__2 = *n - k;

           /* Send the vector to the GPU                                    */
           magma_ssetmatrix( i__1, 1, A(rk, k), *lda, dA(rk,k), *ldda );

           /* Multiply on GPU                                               */
           magma_int_t i__3 = *nb-k;
           magma_int_t i__4 = i__2 -i__3, i__5 = *nb-k+1;
           magma_sgemv(MagmaTrans, i__1 - i__5, i__2 - i__3, tau[k],
                       dA(rk+i__5, k+1+i__3), *ldda, dA(rk+i__5, k),
                       c__1, zero, &df[k + 1 +i__3 + k * f_dim1], c__1);

           magma_sgetmatrix_async(i__2-i__3, 1,
                                  &df[k + 1 +i__3+k * f_dim1], i__2,
                                  & f[k + 1 +i__3+k * f_dim1], i__2, stream);

           blasf77_sgemv(MagmaTransStr, &i__1, &i__3, &tau[k],
                         A(rk, k + 1), lda, A(rk, k),
                         &c__1, &zero, &f[k + 1 + k * f_dim1], &c__1);

           magma_queue_sync( stream );
           blasf77_sgemv(MagmaTransStr, &i__5, &i__4, &tau[k],
                         A(rk, k+1 +i__3), lda, A(rk, k),
                         &c__1, &one, &f[k + 1 + i__3 + k * f_dim1], &c__1);
       }
       
       /* Padding F(1:K,K) with zeros. */
       for (j = 1; j <= k; ++j)
           f[j + k * f_dim1] = zero;
       
       /* Incremental updating of F:   
          F(1:N,K) := F(1:N,K) - tau(K)*F(1:N,1:K-1)*A(RK:M,1:K-1)'*A(RK:M,K). */
       if (k > 1) {
           i__1 = *m - rk + 1;
           i__2 = k - 1;
           float temporary1 = MAGMA_S_REAL(tau[k]);
           float temporary2 = MAGMA_S_IMAG(tau[k]);
           z__1 = MAGMA_S_MAKE(-temporary1, -temporary2);
           blasf77_sgemv(MagmaTransStr, &i__1, &i__2, &z__1, A(rk, 1),
                         lda, A(rk,k), &c__1, &zero, &auxv[1], &c__1);
           
           i__1 = k - 1;
           blasf77_sgemv(MagmaNoTransStr, n, &i__1, &one, &f[f_dim1 + 1], ldf, 
                         &auxv[1], &c__1, &one, &f[k * f_dim1 + 1], &c__1);
       }
       
       /* Optimization: On the last iteration start sending F back to the GPU */

       /* Update the current row of A:   
          A(RK,K+1:N) := A(RK,K+1:N) - A(RK,1:K)*F(K+1:N,1:K)'.               */
       if (k < *n) {
           i__1 = *n - k;
           blasf77_sgemm(MagmaNoTransStr, MagmaTransStr, &c__1, &i__1, &k,
                         &mone, A(rk, 1), lda, &f[k + 1 + f_dim1], ldf,
                         &one, A(rk, k + 1), lda);
       }
       
       /* Update partial column norms. */
       if (rk < lastrk) {
           for (j = k + 1; j <= *n; ++j) {
               if (vn1[j] != 0.) {                   
                   /* NOTE: The following 4 lines follow from the analysis in   
                      Lapack Working Note 176. */
                   temp = MAGMA_S_ABS( *A(rk,j) ) / vn1[j];
                   temp = max(0., ((temp + 1.) * (1. - temp)));

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
       }
       
       *A(rk, k) = akk;
       
       /* End of while loop. */
       goto L10;
    }
    *kb = k;
    rk = *offset + *kb;

    /* Apply the block reflector to the rest of the matrix:   
       A(OFFSET+KB+1:M,KB+1:N) := A(OFFSET+KB+1:M,KB+1:N) - A(OFFSET+KB+1:M,1:KB)*F(KB+1:N,1:KB)'  */
    if (*kb < min(*n, *m - *offset)) {
       i__1 = *m - rk;
       i__2 = *n - *kb;
       
       /* Send f to the GPU                             */
       magma_ssetmatrix( i__2, *kb,
                         &f [*kb + 1 + f_dim1], *ldf,
                         &df[*kb + 1 + f_dim1], i__2);

       magma_sgemm(MagmaNoTrans, MagmaTrans, i__1, i__2, *kb, mone,
                   dA(rk+1,1), *ldda, &df[*kb + 1 + f_dim1], i__2, one,
                   dA(rk+1, *kb + 1), *ldda);
    }
    
    /* Recomputation of difficult columns. */
 L60:
    if (lsticc > 0) {
       itemp = (magma_int_t)(vn2[lsticc] >= 0. ? floor(vn2[lsticc] + .5) : -floor(.5 - vn2[lsticc]));  
       i__1 = *m - rk;
       if (lsticc <= *nb)
           vn1[lsticc] = cblas_snrm2(i__1, A(rk + 1, lsticc), c__1);
       else {
           /* Where is the data, CPU or GPU ? */
           float r1, r2;

           r1 = cblas_snrm2(*nb-k, A(rk + 1, lsticc), c__1);
           r2 = cublasSnrm2(*m-*offset-*nb, dA(*offset + *nb + 1, lsticc), c__1);

           //vn1[lsticc] = cublasSnrm2(i__1, dA(rk + 1, lsticc), c__1);
           vn1[lsticc] = magma_ssqrt(r1*r1+r2*r2);
       }
   
       /* NOTE: The computation of VN1( LSTICC ) relies on the fact that   
          SNRM2 does not fail on vectors with norm below the value of SQRT(SLAMCH('S')) */       
       vn2[lsticc] = vn1[lsticc];
       lsticc = itemp;
       goto L60;
    }
    
    magma_queue_destroy( stream );

    return MAGMA_SUCCESS;
} /* magma_slaqps */
