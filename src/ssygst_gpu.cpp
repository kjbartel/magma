/*
    -- MAGMA (version 1.4.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       December 2013

       @author Raffaele Solca
       @author Azzam Haidar

       @generated s Tue Dec 17 13:18:36 2013
*/

#include "common_magma.h"

#define A(i, j) (w + (j)*lda + (i))
#define B(i, j) (w+nb*lda + (j)*ldb + (i))

#define dA(i, j) (da + (j)*ldda + (i))
#define dB(i, j) (db + (j)*lddb + (i))

extern "C" magma_int_t
magma_ssygst_gpu(magma_int_t itype, char uplo, magma_int_t n,
                 float *da, magma_int_t ldda,
                 float *db, magma_int_t lddb, magma_int_t *info)
{
/*  -- MAGMA (version 1.4.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       December 2013

    Purpose
    =======
    SSYGST_GPU reduces a real symmetric-definite generalized
    eigenproblem to standard form.
    
    If ITYPE = 1, the problem is A*x = lambda*B*x,
    and A is overwritten by inv(U**T)*A*inv(U) or inv(L)*A*inv(L**T)
    
    If ITYPE = 2 or 3, the problem is A*B*x = lambda*x or
    B*A*x = lambda*x, and A is overwritten by U*A*U**T or L**T*A*L.
    
    B must have been previously factorized as U**T*U or L*L**T by SPOTRF.
    
    Arguments
    =========
    ITYPE   (input) INTEGER
            = 1: compute inv(U**T)*A*inv(U) or inv(L)*A*inv(L**T);
            = 2 or 3: compute U*A*U**T or L**T*A*L.
    
    UPLO    (input) CHARACTER*1
            = 'U':  Upper triangle of A is stored and B is factored as
                    U**T*U;
            = 'L':  Lower triangle of A is stored and B is factored as
                    L*L**T.
    
    N       (input) INTEGER
            The order of the matrices A and B.  N >= 0.
    
    DA      (device input/output) REAL array, dimension (LDA,N)
            On entry, the symmetric matrix A.  If UPLO = 'U', the leading
            N-by-N upper triangular part of A contains the upper
            triangular part of the matrix A, and the strictly lower
            triangular part of A is not referenced.  If UPLO = 'L', the
            leading N-by-N lower triangular part of A contains the lower
            triangular part of the matrix A, and the strictly upper
            triangular part of A is not referenced.
    
            On exit, if INFO = 0, the transformed matrix, stored in the
            same format as A.
    
    LDDA    (input) INTEGER
            The leading dimension of the array A.  LDA >= max(1,N).
    
    DB      (device input) REAL array, dimension (LDB,N)
            The triangular factor from the Cholesky factorization of B,
            as returned by SPOTRF.
    
    LDDB    (input) INTEGER
            The leading dimension of the array B.  LDB >= max(1,N).
    
    INFO    (output) INTEGER
            = 0:  successful exit
            < 0:  if INFO = -i, the i-th argument had an illegal value
    =====================================================================*/
    
    char uplo_[2] = {uplo, 0};
    magma_int_t        nb;
    magma_int_t        k, kb, kb2;
    float    c_one      = MAGMA_S_ONE;
    float    c_neg_one  = MAGMA_S_NEG_ONE;
    float    c_half     = MAGMA_S_HALF;
    float    c_neg_half = MAGMA_S_NEG_HALF;
    float   *w;
    magma_int_t        lda;
    magma_int_t        ldb;
    float             d_one = 1.0;
    int upper = lapackf77_lsame(uplo_, "U");
    
    /* Test the input parameters. */
    *info = 0;
    if (itype<1 || itype>3){
        *info = -1;
    }else if ((! upper) && (! lapackf77_lsame(uplo_, "L"))) {
        *info = -2;
    } else if (n < 0) {
        *info = -3;
    } else if (ldda < max(1,n)) {
        *info = -5;
    }else if (lddb < max(1,n)) {
        *info = -7;
    }
    if (*info != 0) {
        magma_xerbla( __func__, -(*info) );
        return *info;
    }
    
    /* Quick return */
    if ( n == 0 )
        return *info;
    
    nb = magma_get_ssygst_nb(n);
    
    lda = nb;
    ldb = nb;
    
    if (MAGMA_SUCCESS != magma_smalloc_pinned( &w, 2*nb*nb )) {
        *info = MAGMA_ERR_DEVICE_ALLOC;
        return *info;
    }
    
    magma_queue_t stream[3];
    magma_queue_create( &stream[0] );
    magma_queue_create( &stream[1] );
    magma_queue_create( &stream[2] );
    
    /* Use hybrid blocked code */
    if (itype==1) {
        if (upper) {
            kb = min(n,nb);
            
            /* Compute inv(U')*A*inv(U) */
            magma_sgetmatrix_async( kb, kb,
                                    dB(0, 0), lddb,
                                    B(0, 0),  nb, stream[2] );
            magma_sgetmatrix_async( kb, kb,
                                    dA(0, 0), ldda,
                                    A(0, 0),  nb, stream[1] );
            
            for(k = 0; k<n; k+=nb){
                kb = min(n-k,nb);
                kb2= min(n-k-nb,nb);
                
                /* Update the upper triangle of A(k:n,k:n) */
                
                magma_queue_sync( stream[2] );
                magma_queue_sync( stream[1] );
                
                lapackf77_ssygst( &itype, uplo_, &kb, A(0,0), &lda, B(0,0), &ldb, info);
                
                magma_ssetmatrix_async( kb, kb,
                                        A(0, 0),  lda,
                                        dA(k, k), ldda, stream[0] );
                
                if(k+kb<n){
                    // Start copying the new B block
                    magma_sgetmatrix_async( kb2, kb2,
                                            dB(k+kb, k+kb), lddb,
                                            B(0, 0),        nb, stream[2] );
                    
                    magma_strsm(MagmaLeft, MagmaUpper, MagmaTrans, MagmaNonUnit,
                                kb, n-k-kb,
                                c_one, dB(k,k), lddb,
                                dA(k,k+kb), ldda);
                    
                    magma_queue_sync( stream[0] );
                    
                    magma_ssymm(MagmaLeft, MagmaUpper,
                                kb, n-k-kb,
                                c_neg_half, dA(k,k), ldda,
                                dB(k,k+kb), lddb,
                                c_one, dA(k, k+kb), ldda);
                    
                    magma_ssyr2k(MagmaUpper, MagmaTrans,
                                 n-k-kb, kb,
                                 c_neg_one, dA(k,k+kb), ldda,
                                 dB(k,k+kb), lddb,
                                 d_one, dA(k+kb,k+kb), ldda);
                    
                    magma_sgetmatrix_async( kb2, kb2,
                                            dA(k+kb, k+kb), ldda,
                                            A(0, 0),        lda, stream[1] );
                    
                    magma_ssymm(MagmaLeft, MagmaUpper,
                                kb, n-k-kb,
                                c_neg_half, dA(k,k), ldda,
                                dB(k,k+kb), lddb,
                                c_one, dA(k, k+kb), ldda);
                    
                    magma_strsm(MagmaRight, MagmaUpper, MagmaNoTrans, MagmaNonUnit,
                                kb, n-k-kb,
                                c_one ,dB(k+kb,k+kb), lddb,
                                dA(k,k+kb), ldda);
                }
            }
            
            magma_queue_sync( stream[0] );
        }
        else {
            kb = min(n,nb);
            
            /* Compute inv(L)*A*inv(L') */
            magma_sgetmatrix_async( kb, kb,
                                    dB(0, 0), lddb,
                                    B(0, 0),  nb, stream[2] );
            magma_sgetmatrix_async( kb, kb,
                                    dA(0, 0), ldda,
                                    A(0, 0),  nb, stream[1] );
            
            for(k = 0; k<n; k+=nb){
                kb= min(n-k,nb);
                kb2= min(n-k-nb,nb);
                
                /* Update the lower triangle of A(k:n,k:n) */
                
                magma_queue_sync( stream[2] );
                magma_queue_sync( stream[1] );
                
                lapackf77_ssygst( &itype, uplo_, &kb, A(0, 0), &lda, B(0, 0), &ldb, info);
                
                magma_ssetmatrix_async( kb, kb,
                                        A(0, 0),  lda,
                                        dA(k, k), ldda, stream[0] );
                
                if(k+kb<n){
                    // Start copying the new B block
                    magma_sgetmatrix_async( kb2, kb2,
                                            dB(k+kb, k+kb), lddb,
                                            B(0, 0),        nb, stream[2] );
                    
                    magma_strsm(MagmaRight, MagmaLower, MagmaTrans, MagmaNonUnit,
                                n-k-kb, kb,
                                c_one, dB(k,k), lddb,
                                dA(k+kb,k), ldda);
                    
                    magma_queue_sync( stream[0] );
                    
                    magma_ssymm(MagmaRight, MagmaLower,
                                n-k-kb, kb,
                                c_neg_half, dA(k,k), ldda,
                                dB(k+kb,k), lddb,
                                c_one, dA(k+kb, k), ldda);
                    
                    magma_ssyr2k(MagmaLower, MagmaNoTrans,
                                 n-k-kb, kb,
                                 c_neg_one, dA(k+kb,k), ldda,
                                 dB(k+kb,k), lddb,
                                 d_one, dA(k+kb,k+kb), ldda);
                    
                    magma_sgetmatrix_async( kb2, kb2,
                                            dA(k+kb, k+kb), ldda,
                                            A(0, 0),        lda, stream[1] );
                    
                    magma_ssymm(MagmaRight, MagmaLower,
                                n-k-kb, kb,
                                c_neg_half, dA(k,k), ldda,
                                dB(k+kb,k), lddb,
                                c_one, dA(k+kb, k), ldda);
                    
                    magma_strsm(MagmaLeft, MagmaLower, MagmaNoTrans, MagmaNonUnit,
                                n-k-kb, kb,
                                c_one, dB(k+kb,k+kb), lddb,
                                dA(k+kb,k), ldda);
                }
            }
        }
        
        magma_queue_sync( stream[0] );
    }
    else {
        if (upper) {
            /* Compute U*A*U' */
            for(k = 0; k<n; k+=nb){
                kb= min(n-k,nb);
                
                magma_sgetmatrix_async( kb, kb,
                                        dB(k, k), lddb,
                                        B(0, 0),  nb, stream[2] );
                
                /* Update the upper triangle of A(1:k+kb-1,1:k+kb-1) */
                if(k>0){
                    magma_strmm(MagmaLeft, MagmaUpper, MagmaNoTrans, MagmaNonUnit,
                                k, kb,
                                c_one ,dB(0,0), lddb,
                                dA(0,k), ldda);
                    
                    magma_ssymm(MagmaRight, MagmaUpper,
                                k, kb,
                                c_half, dA(k,k), ldda,
                                dB(0,k), lddb,
                                c_one, dA(0, k), ldda);
                    
                    magma_queue_sync( stream[1] );
                }
                
                magma_sgetmatrix_async( kb, kb,
                                        dA(k, k), ldda,
                                        A(0, 0),  lda, stream[0] );
                
                if(k>0){
                    magma_ssyr2k(MagmaUpper, MagmaNoTrans,
                                 k, kb,
                                 c_one, dA(0,k), ldda,
                                 dB(0,k), lddb,
                                 d_one, dA(0,0), ldda);
                    
                    magma_ssymm(MagmaRight, MagmaUpper,
                                k, kb,
                                c_half, dA(k,k), ldda,
                                dB(0,k), lddb,
                                c_one, dA(0, k), ldda);
                    
                    magma_strmm(MagmaRight, MagmaUpper, MagmaTrans, MagmaNonUnit,
                                k, kb,
                                c_one, dB(k,k), lddb,
                                dA(0,k), ldda);
                }
                
                magma_queue_sync( stream[2] );
                magma_queue_sync( stream[0] );
                
                lapackf77_ssygst( &itype, uplo_, &kb, A(0, 0), &lda, B(0, 0), &ldb, info);
                
                magma_ssetmatrix_async( kb, kb,
                                        A(0, 0),  lda,
                                        dA(k, k), ldda, stream[1] );
            }
            
            magma_queue_sync( stream[1] );
        }
        else {
            /* Compute L'*A*L */
            for(k = 0; k<n; k+=nb){
                kb= min(n-k,nb);
                
                magma_sgetmatrix_async( kb, kb,
                                        dB(k, k), lddb,
                                        B(0, 0),  nb, stream[2] );
                
                /* Update the lower triangle of A(1:k+kb-1,1:k+kb-1) */
                if(k>0){
                    magma_strmm(MagmaRight, MagmaLower, MagmaNoTrans, MagmaNonUnit,
                                kb, k,
                                c_one ,dB(0,0), lddb,
                                dA(k,0), ldda);
                    
                    magma_ssymm(MagmaLeft, MagmaLower,
                                kb, k,
                                c_half, dA(k,k), ldda,
                                dB(k,0), lddb,
                                c_one, dA(k, 0), ldda);
                    
                    magma_queue_sync( stream[1] );
                }
                
                magma_sgetmatrix_async( kb, kb,
                                        dA(k, k), ldda,
                                        A(0, 0),  lda, stream[0] );
                
                if(k>0){
                    magma_ssyr2k(MagmaLower, MagmaTrans,
                                 k, kb,
                                 c_one, dA(k,0), ldda,
                                 dB(k,0), lddb,
                                 d_one, dA(0,0), ldda);
                    
                    magma_ssymm(MagmaLeft, MagmaLower,
                                kb, k,
                                c_half, dA(k,k), ldda,
                                dB(k,0), lddb,
                                c_one, dA(k, 0), ldda);
                    
                    magma_strmm(MagmaLeft, MagmaLower, MagmaTrans, MagmaNonUnit,
                                kb, k,
                                c_one, dB(k,k), lddb,
                                dA(k,0), ldda);
                }
                
                magma_queue_sync( stream[2] );
                magma_queue_sync( stream[0] );
                
                lapackf77_ssygst( &itype, uplo_, &kb, A(0, 0), &lda, B(0, 0), &ldb, info);
                
                magma_ssetmatrix_async( kb, kb,
                                        A(0, 0),  lda,
                                        dA(k, k), ldda, stream[1] );
            }
            
            magma_queue_sync( stream[1] );
        }
    }
    magma_queue_destroy( stream[0] );
    magma_queue_destroy( stream[1] );
    magma_queue_destroy( stream[2] );
    
    magma_free_pinned( w );
    
    return *info;
} /* magma_ssygst_gpu */
#undef A
#undef B
#undef dA
#undef dB

