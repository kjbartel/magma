/*
    -- MAGMA (version 1.3.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2012

       @author Raffaele Solca

       @generated s Wed Nov 14 22:53:30 2012
*/
#include "common_magma.h"
#include <cblas.h>

extern "C"
magma_int_t magma_get_ssygst_m_nb() { return 256;}

#define A(i, j) (a+(j)*nb*lda + (i)*nb)
#define B(i, j) (b+(j)*nb*ldb + (i)*nb)

#define dA(gpui, i, j) (dw[gpui] + (j)*nb*ldda + (i)*nb)
#define dB_c(gpui, i, j) (dw[gpui] + dima*ldda + (i)*nb + (j)*nb*lddbc)
#define dB_r(gpui, i, j) (dw[gpui] + dima*ldda + (i)*nb + (j)*nb*lddbr)

extern "C" magma_int_t
magma_ssygst_m(magma_int_t nrgpu, magma_int_t itype, char uplo, magma_int_t n,
               float *a, magma_int_t lda,
               float *b, magma_int_t ldb, magma_int_t *info)
{
/*
  -- MAGMA (version 1.3.0) --
     Univ. of Tennessee, Knoxville
     Univ. of California, Berkeley
     Univ. of Colorado, Denver
     November 2012


   Purpose
   =======

   SSYGST reduces a real symmetric-definite generalized
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

   A       (input/output) COMPLEX*16 array, dimension (LDA,N)
           On entry, the symmetric matrix A.  If UPLO = 'U', the leading
           N-by-N upper triangular part of A contains the upper
           triangular part of the matrix A, and the strictly lower
           triangular part of A is not referenced.  If UPLO = 'L', the
           leading N-by-N lower triangular part of A contains the lower
           triangular part of the matrix A, and the strictly upper
           triangular part of A is not referenced.

           On exit, if INFO = 0, the transformed matrix, stored in the
           same format as A.

   LDA     (input) INTEGER
           The leading dimension of the array A.  LDA >= max(1,N).

   B       (input) COMPLEX*16 array, dimension (LDB,N)
           The triangular factor from the Cholesky factorization of B,
           as returned by SPOTRF.

   LDB     (input) INTEGER
           The leading dimension of the array B.  LDB >= max(1,N).

   INFO    (output) INTEGER
           = 0:  successful exit
           < 0:  if INFO = -i, the i-th argument had an illegal value

   =====================================================================*/

    char uplo_[2] = {uplo, 0};
    magma_int_t        k, kb, j, jb, kb2;
    magma_int_t        ldda, dima, lddbr, lddbc;
    float    c_one      = MAGMA_S_ONE;
    float    c_neg_one  = MAGMA_S_NEG_ONE;
    float    c_half     = MAGMA_S_HALF;
    float    c_neg_half = MAGMA_S_NEG_HALF;
    float* dw[MagmaMaxGPUs];
    cudaStream_t stream [MagmaMaxGPUs][3];
    magma_int_t igpu = 0;

    int gpu_b;
    magma_getdevice(&gpu_b);

    float             d_one = 1.0;
    int upper = lapackf77_lsame(uplo_, "U");

    magma_int_t nb = magma_get_ssygst_m_nb();

    /* Test the input parameters. */
    *info = 0;
    if (itype<1 || itype>3){
        *info = -1;
    }else if ((! upper) && (! lapackf77_lsame(uplo_, "L"))) {
        *info = -2;
    } else if (n < 0) {
        *info = -3;
    } else if (lda < max(1,n)) {
        *info = -5;
    }else if (ldb < max(1,n)) {
        *info = -7;
    }
    if (*info != 0) {
        magma_xerbla( __func__, -(*info) );
        return *info;
    }

    /* Quick return */
    if ( n == 0 )
        return *info;

    magma_int_t nbl = (n-1)/nb+1; // number of blocks

    if ( (itype==1 && upper) || (itype!=1 && !upper) ){
        ldda = ((nbl-1)/nrgpu+1)*nb;
        dima = n;
    } else {
        ldda = n;
        dima = ((nbl-1)/nrgpu+1)*nb;
    }
    lddbr = 2 * nb;
    lddbc = n;
    for (igpu = 0; igpu < nrgpu; ++igpu){
        magma_setdevice(igpu);
        if (MAGMA_SUCCESS != magma_smalloc( &dw[igpu], (dima*ldda + lddbc*lddbr) )) {
            *info = MAGMA_ERR_DEVICE_ALLOC;
            return *info;
        }
        magma_queue_create( &stream[igpu][0] );
        magma_queue_create( &stream[igpu][1] );
        magma_queue_create( &stream[igpu][2] );
    }

    /* Use hybrid blocked code */

    if (itype==1) {
        if (upper) {

            /* Compute inv(U')*A*inv(U) */

            //copy A to mgpus
            for (k = 0; k < nbl; ++k){
                igpu = k%nrgpu;
                magma_setdevice(igpu);
                kb = min(nb, n-k*nb);
                magma_ssetmatrix_async(kb, n-k*nb,
                                       A(k, k),              lda,
                                       dA(igpu, k/nrgpu, k), ldda, stream[igpu][0] );
            }
            kb= min(n,nb);
            igpu = 0;
            magma_setdevice(igpu);
            // dB_r(0,0) is used to store B(k,k)
            magma_ssetmatrix_async(kb, kb,
                                   B(0, 0),          ldb,
                                   dB_r(igpu, 0, 0), lddbr, stream[igpu][1] );

            for(k = 0; k<nbl; ++k){
                kb= min(n-k*nb,nb);
                kb2= min(n-(k+1)*nb,nb);

                if(k+1<nbl){
                    for (igpu = 0; igpu < nrgpu; ++igpu){
                        magma_setdevice(igpu);
                        magma_queue_sync( stream[igpu][0] );
                        magma_ssetmatrix_async(kb, n-(k+1)*nb,
                                               B(k, k+1),          ldb,
                                               dB_r(igpu, 0, k+1), lddbr, stream[igpu][0] );
                    }
                }

                igpu = k%nrgpu;
                magma_setdevice(igpu);

                magma_queue_sync( stream[igpu][1] ); // Needed, otherwise conflicts reading B(k,k) between sygs2 and cudaMemcpy2D
                magma_queue_sync( stream[igpu][2] );

                if(k+1<nbl){
                    magmablasSetKernelStream(stream[igpu][1]);
                    // dB_r(0,0) stores B(k,k)
                    magma_strsm(MagmaLeft, uplo, MagmaTrans, MagmaNonUnit,
                                kb, n-(k+1)*nb,
                                c_one, dB_r(igpu, 0, 0), lddbr,
                                dA(igpu, k/nrgpu, k+1), ldda);
                }

                lapackf77_ssygs2( &itype, uplo_, &kb, A(k,k), &lda, B(k,k), &ldb, info);

                if (k+1<nbl) {
                    magma_ssetmatrix_async(kb, kb,
                                           A(k, k),              lda,
                                           dA(igpu, k/nrgpu, k), ldda, stream[igpu][0] );

                    magma_queue_sync( stream[igpu][1] );
                    magmablasSetKernelStream(stream[igpu][0]);

                    magma_ssymm(MagmaLeft, uplo,
                                kb, n-(k+1)*nb,
                                c_neg_half, dA(igpu, k/nrgpu, k), ldda,
                                dB_r(igpu, 0, k+1), lddbr,
                                c_one, dA(igpu, k/nrgpu, k+1), ldda);

                    magma_queue_sync( stream[igpu][0] );

                    magma_sgetmatrix(kb, n-(k+1)*nb,
                                     dA(igpu, k/nrgpu, k+1), ldda,
                                     A(k, k+1),              lda );

                    // send the partially updated panel of dA to each gpu in the second dB block
                    // to overlap symm computation

                    for (igpu = 0; igpu < nrgpu; ++igpu){
                        magma_setdevice(igpu);
                        magma_ssetmatrix_async(kb, n-(k+1)*nb,
                                               A(k, k+1),          lda,
                                               dB_r(igpu, 1, k+1), lddbr, stream[igpu][0] );
                    }

                    igpu = k%nrgpu;
                    magma_setdevice(igpu);
                    magmablasSetKernelStream(stream[igpu][1]);

                    magma_ssymm(MagmaLeft, uplo,
                                kb, n-(k+1)*nb,
                                c_neg_half, dA(igpu, k/nrgpu, k), ldda,
                                dB_r(igpu, 0, k+1), lddbr,
                                c_one, dA(igpu, k/nrgpu, k+1), ldda);

                    for (igpu = 0; igpu < nrgpu; ++igpu){
                        magma_queue_sync( stream[igpu][0] );
                    }

                    for (j = k+1; j < nbl; ++j){
                        jb = min(nb, n-j*nb);
                        igpu = j%nrgpu;
                        magma_setdevice(igpu);
                        magmablasSetKernelStream(stream[igpu][(j/nrgpu)%3]);
                        magma_ssyr2k(uplo, MagmaTrans,
                                     jb, nb,
                                     c_neg_one, dB_r(igpu, 1, j), lddbr,
                                     dB_r(igpu, 0, j), lddbr,
                                     d_one, dA(igpu, j/nrgpu, j), ldda);
                        magma_queue_sync( stream[igpu][((j)/nrgpu)%3] ); // Needed for correctness. Why?
                        if (j == k+1){
                            magma_queue_sync( stream[igpu][(j/nrgpu)%3] );
                            magma_sgetmatrix_async(kb2, kb2,
                                                   dA(igpu, (k+1)/nrgpu, k+1), ldda,
                                                   A(k+1, k+1),                lda, stream[igpu][2] );
                            // dB_r(0,0) is used to store B(k,k)
                            magma_ssetmatrix_async(kb2, kb2,
                                                   B(k+1, k+1),      ldb,
                                                   dB_r(igpu, 0, 0), lddbr, stream[igpu][1] );
                        }
                    }
                    for (j = k+1; j < nbl-1; ++j){
                        igpu = j%nrgpu;
                        magma_setdevice(igpu);
                        magmablasSetKernelStream(stream[igpu][0]);
                        magma_sgemm(MagmaTrans, MagmaNoTrans, nb, n-(j+1)*nb, nb, c_neg_one, dB_r(igpu, 0, j), lddbr,
                                    dB_r(igpu, 1, j+1), lddbr, c_one, dA(igpu, j/nrgpu, j+1), ldda );

                        magma_sgemm(MagmaTrans, MagmaNoTrans, nb, n-(j+1)*nb, nb, c_neg_one, dB_r(igpu, 1, j), lddbr,
                                    dB_r(igpu, 0, j+1), lddbr, c_one, dA(igpu, j/nrgpu, j+1), ldda );
                    }
                }
            }

            for (igpu = 0; igpu < nrgpu; ++igpu){
                magma_queue_sync( stream[igpu][0] );
                magma_queue_sync( stream[igpu][1] );
            }

            if (n > nb){

                magma_int_t nloc[MagmaMaxGPUs];

                jb = min(nb, n-nb);
                for (igpu = 0; igpu < nrgpu; ++igpu){
                    nloc[igpu]=0;
                    magma_setdevice(igpu);
                    magma_ssetmatrix_async(jb, n-nb,
                                           B(1, 1),          ldb,
                                           dB_r(igpu, 1, 1), lddbr, stream[igpu][1] );
                }
                for (j = 1; j < nbl; ++j){
                    if ((j+1)*nb < n){
                        jb = min(nb, n-(j+1)*nb);
                        for (igpu = 0; igpu < nrgpu; ++igpu){
                            magma_setdevice(igpu);
                            magma_ssetmatrix_async(jb, n-(j+1)*nb,
                                                   B(j+1, j+1),              ldb,
                                                   dB_r(igpu, (j+1)%2, j+1), lddbr, stream[igpu][(j+1)%2] );
                        }
                    }
                    jb = min(nb, n-j*nb);
                    nloc[(j-1)%nrgpu] += nb;

                    for (igpu = 0; igpu < nrgpu; ++igpu){
                        magma_setdevice(igpu);
                        magmablasSetKernelStream(stream[igpu][j%2]);
                        magma_strsm(MagmaRight, uplo, MagmaNoTrans, MagmaNonUnit, nloc[igpu], jb, c_one, dB_r(igpu, j%2, j), lddbr,
                                    dA(igpu, 0, j), ldda );
                    }

                    if ( j < nbl-1 ){

                        for (igpu = 0; igpu < nrgpu; ++igpu){
                            magma_setdevice(igpu);
                            magmablasSetKernelStream(stream[igpu][j%2]);
                            magma_sgemm(MagmaNoTrans, MagmaNoTrans, nloc[igpu], n-(j+1)*nb, nb, c_neg_one, dA(igpu, 0, j), ldda,
                                        dB_r(igpu, j%2, j+1), lddbr, c_one, dA(igpu, 0, j+1), ldda );
                        }
                    }

                    for (igpu = 0; igpu < nrgpu; ++igpu){
                        magma_queue_sync( stream[igpu][j%2] );
                    }

                    for (k = 0; k < j; ++k){
                        igpu = k%nrgpu;
                        magma_setdevice(igpu);
                        kb = min(nb, n-k*nb);
                        magma_sgetmatrix_async(kb, jb,
                                               dA(igpu, k/nrgpu, j), ldda,
                                               A(k, j),              lda, stream[igpu][2] );
                    }
                }
            }


        } else {
            /* Compute inv(L)*A*inv(L') */
            //copy A to mgpus
            for (k = 0; k < nbl; ++k){
                igpu = k%nrgpu;
                magma_setdevice(igpu);
                kb = min(nb, n-k*nb);
                magma_ssetmatrix_async((n-k*nb), kb,
                                       A(k, k),              lda,
                                       dA(igpu, k, k/nrgpu), ldda, stream[igpu][0] );
            }
            kb= min(n,nb);
            igpu = 0;
            magma_setdevice(igpu);
            // dB_c(0,0) is used to store B(k,k)
            magma_ssetmatrix_async(kb, kb,
                                   B(0, 0),          ldb,
                                   dB_c(igpu, 0, 0), lddbc, stream[igpu][1] );

            for(k = 0; k<nbl; ++k){
                kb= min(n-k*nb,nb);
                kb2= min(n-(k+1)*nb,nb);

                if(k+1<nbl){
                    for (igpu = 0; igpu < nrgpu; ++igpu){
                        magma_setdevice(igpu);
                        magma_queue_sync( stream[igpu][0] );
                        magma_ssetmatrix_async((n-(k+1)*nb), kb,
                                               B(k+1, k),          ldb,
                                               dB_c(igpu, k+1, 0), lddbc, stream[igpu][0] );
                    }
                }

                igpu = k%nrgpu;
                magma_setdevice(igpu);

                magma_queue_sync( stream[igpu][1] ); // Needed, otherwise conflicts reading B(k,k) between sygs2 and cudaMemcpy2D
                magma_queue_sync( stream[igpu][2] );

                if(k+1<nbl){
                    magmablasSetKernelStream(stream[igpu][1]);
                    // dB_c(0,0) stores B(k,k)
                    magma_strsm(MagmaRight, uplo, MagmaTrans, MagmaNonUnit,
                                n-(k+1)*nb, kb,
                                c_one, dB_c(igpu, 0, 0), lddbc,
                                dA(igpu, k+1, k/nrgpu), ldda);
                }

                lapackf77_ssygs2( &itype, uplo_, &kb, A(k,k), &lda, B(k,k), &ldb, info);

                if (k+1<nbl) {
                    magma_ssetmatrix_async(kb, kb,
                                           A(k, k),               lda,
                                           dA(igpu, k , k/nrgpu), ldda, stream[igpu][0] );

                    magma_queue_sync( stream[igpu][1] );
                    magmablasSetKernelStream(stream[igpu][0]);

                    magma_ssymm(MagmaRight, uplo,
                                n-(k+1)*nb, kb,
                                c_neg_half, dA(igpu, k, k/nrgpu), ldda,
                                dB_c(igpu, k+1, 0), lddbc,
                                c_one, dA(igpu, k+1, k/nrgpu), ldda);

                    magma_queue_sync( stream[igpu][0] );

                    magma_sgetmatrix(n-(k+1)*nb, kb,
                                     dA(igpu, k+1, k/nrgpu), ldda,
                                     A(k+1, k),              lda );

                    // send the partially updated panel of dA to each gpu in the second dB block
                    // to overlap symm computation

                    for (igpu = 0; igpu < nrgpu; ++igpu){
                        magma_setdevice(igpu);
                        magma_ssetmatrix_async((n-(k+1)*nb), kb,
                                               A(k+1, k),          lda,
                                               dB_c(igpu, k+1, 1), lddbc, stream[igpu][0] );
                    }

                    igpu = k%nrgpu;
                    magma_setdevice(igpu);
                    magmablasSetKernelStream(stream[igpu][1]);

                    magma_ssymm(MagmaRight, uplo,
                                n-(k+1)*nb, kb,
                                c_neg_half, dA(igpu, k, k/nrgpu), ldda,
                                dB_c(igpu, k+1, 0), lddbc,
                                c_one, dA(igpu, k+1, k/nrgpu), ldda);

                    for (igpu = 0; igpu < nrgpu; ++igpu){
                        magma_queue_sync( stream[igpu][0] );
                    }

                    for (j = k+1; j < nbl; ++j){
                        jb = min(nb, n-j*nb);
                        igpu = j%nrgpu;
                        magma_setdevice(igpu);
                        magmablasSetKernelStream(stream[igpu][(j/nrgpu)%3]);
                        magma_ssyr2k(uplo, MagmaNoTrans,
                                     jb, nb,
                                     c_neg_one, dB_c(igpu, j, 1), lddbc,
                                     dB_c(igpu, j, 0), lddbc,
                                     d_one, dA(igpu, j, j/nrgpu), ldda);
                        magma_queue_sync( stream[igpu][((j)/nrgpu)%3] ); // Needed for correctness. Why?
                        if (j == k+1){
                            magma_queue_sync( stream[igpu][(j/nrgpu)%3] );
                            magma_sgetmatrix_async(kb2, kb2,
                                                   dA(igpu, k+1, (k+1)/nrgpu), ldda,
                                                   A(k+1, k+1),                lda, stream[igpu][2] );
                            // dB_c(0,0) is used to store B(k,k)
                            magma_ssetmatrix_async(kb2, kb2,
                                                   B(k+1, k+1),      ldb,
                                                   dB_c(igpu, 0, 0), lddbc, stream[igpu][1] );
                        }
                    }
                    for (j = k+1; j < nbl-1; ++j){
                        igpu = j%nrgpu;
                        magma_setdevice(igpu);
                        magmablasSetKernelStream(stream[igpu][0]);
                        magma_sgemm(MagmaNoTrans, MagmaTrans, n-(j+1)*nb, nb, nb, c_neg_one, dB_c(igpu, j+1, 1), lddbc,
                                    dB_c(igpu, j, 0), lddbc, c_one, dA(igpu, j+1, j/nrgpu), ldda );

                        magma_sgemm(MagmaNoTrans, MagmaTrans, n-(j+1)*nb, nb, nb, c_neg_one, dB_c(igpu, j+1, 0), lddbc,
                                    dB_c(igpu, j, 1), lddbc, c_one, dA(igpu, j+1, j/nrgpu), ldda );
                    }
                }
            }

            for (igpu = 0; igpu < nrgpu; ++igpu){
                magma_queue_sync( stream[igpu][0] );
                magma_queue_sync( stream[igpu][1] );
            }

            if (n > nb){

                magma_int_t nloc[MagmaMaxGPUs];

                jb = min(nb, n-nb);
                for (igpu = 0; igpu < nrgpu; ++igpu){
                    nloc[igpu]=0;
                    magma_setdevice(igpu);
                    magma_ssetmatrix_async((n-nb), jb,
                                           B(1, 1),          ldb,
                                           dB_c(igpu, 1, 1), lddbc, stream[igpu][1] );
                }
                for (j = 1; j < nbl; ++j){
                    if ((j+1)*nb < n){
                        jb = min(nb, n-(j+1)*nb);
                        for (igpu = 0; igpu < nrgpu; ++igpu){
                            magma_setdevice(igpu);
                            magma_ssetmatrix_async((n-(j+1)*nb), jb,
                                                   B(j+1, j+1),              ldb,
                                                   dB_c(igpu, j+1, (j+1)%2), lddbc, stream[igpu][(j+1)%2] );
                        }
                    }
                    jb = min(nb, n-j*nb);
                    nloc[(j-1)%nrgpu] += nb;

                    for (igpu = 0; igpu < nrgpu; ++igpu){
                        magma_setdevice(igpu);
                        magmablasSetKernelStream(stream[igpu][j%2]);
                        magma_strsm(MagmaLeft, uplo, MagmaNoTrans, MagmaNonUnit, jb, nloc[igpu], c_one, dB_c(igpu, j, j%2), lddbc,
                                    dA(igpu, j, 0), ldda );
                    }

                    if ( j < nbl-1 ){

                        for (igpu = 0; igpu < nrgpu; ++igpu){
                            magma_setdevice(igpu);
                            magmablasSetKernelStream(stream[igpu][j%2]);
                            magma_sgemm(MagmaNoTrans, MagmaNoTrans, n-(j+1)*nb, nloc[igpu], nb, c_neg_one, dB_c(igpu, j+1, j%2), lddbc,
                                        dA(igpu, j, 0), ldda, c_one, dA(igpu, j+1, 0), ldda );
                        }
                    }

                    for (igpu = 0; igpu < nrgpu; ++igpu){
                        magma_queue_sync( stream[igpu][j%2] );
                    }

                    for (k = 0; k < j; ++k){
                        igpu = k%nrgpu;
                        magma_setdevice(igpu);
                        kb = min(nb, n-k*nb);
                        magma_sgetmatrix_async(jb, kb,
                                               dA(igpu, j, k/nrgpu), ldda,
                                               A(j, k),              lda, stream[igpu][2] );
                    }
                }
            }

        }

    } else {

        if (upper) {

            /* Compute U*A*U' */

            if (n > nb){
                magma_int_t nloc[MagmaMaxGPUs];
                magma_int_t iloc[MagmaMaxGPUs];
                for(igpu = 0; igpu < nrgpu; ++igpu){
                    nloc[igpu] = 0;
                    iloc[igpu] = 0;
                }

                kb = min(nb, n);
                for (j = 0; j < nbl; ++j){
                    igpu = j%nrgpu;
                    magma_setdevice(igpu);
                    jb = min(nb, n-j*nb);
                    nloc[igpu] += jb;
                    magma_ssetmatrix_async(kb, jb,
                                           A(0, j),              lda,
                                           dA(igpu, 0, j/nrgpu), ldda, stream[igpu][0] );
                }
                for (igpu = 0; igpu < nrgpu; ++igpu){
                    magma_setdevice(igpu);
                    magma_ssetmatrix_async(kb, kb,
                                           B(0, 0),          ldb,
                                           dB_c(igpu, 0, 0), lddbc, stream[igpu][0] );
                }
                for (k = 0; k < nbl-1; ++k){
                    ++iloc[k%nrgpu];
                    kb = min(nb, n-(k+1)*nb);
                    for (j = k+1; j < nbl; ++j){
                        igpu = j%nrgpu;
                        magma_setdevice(igpu);
                        jb = min(nb, n-j*nb);
                        magma_ssetmatrix_async(kb, jb,
                                               A(k+1, j),              lda,
                                               dA(igpu, k+1, j/nrgpu), ldda, stream[igpu][(k+1)%2] );
                    }
                    for (igpu = 0; igpu < nrgpu; ++igpu){
                        magma_setdevice(igpu);
                        magma_ssetmatrix_async((k+1)*nb + kb, kb,
                                               B(0, k+1),              ldb,
                                               dB_c(igpu, 0, (k+1)%2), lddbc, stream[igpu][(k+1)%2] );
                    }

                    kb = min(nb, n-k*nb);

                    if (k > 0){
                        for (igpu = 0; igpu < nrgpu; ++igpu){
                            magma_setdevice(igpu);
                            magmablasSetKernelStream(stream[igpu][k%2]);
                            magma_sgemm(MagmaNoTrans, MagmaNoTrans, k*nb, nloc[igpu]-nb*iloc[igpu], kb, c_one, dB_c(igpu, 0, k%2), lddbc,
                                        dA(igpu, k, iloc[igpu]), ldda, c_one, dA(igpu, 0, iloc[igpu]), ldda );
                        }
                    }
                    for (igpu = 0; igpu < nrgpu; ++igpu){
                        magma_setdevice(igpu);
                        magmablasSetKernelStream(stream[igpu][k%2]);
                        magma_strmm(MagmaLeft, uplo, MagmaNoTrans, MagmaNonUnit, kb, nloc[igpu]-nb*iloc[igpu], c_one, dB_c(igpu, k, k%2), lddbc,
                                    dA(igpu, k, iloc[igpu]), ldda );
                    }

                    for (igpu = 0; igpu < nrgpu; ++igpu){
                        magma_queue_sync( stream[igpu][k%2] );
                    }
                }
            }

            for(k = 0; k<nbl; ++k){
                kb= min(n-k*nb,nb);
                kb2= min(n-(k+1)*nb,nb);

                if(k>0){

                    igpu = k%nrgpu;
                    magma_setdevice(igpu);
                    magmablasSetKernelStream(stream[igpu][0]);
                    magma_ssymm(MagmaRight, uplo,
                                k*nb, kb,
                                c_half, dA(igpu, k, k/nrgpu), ldda,
                                dB_c(igpu, 0, 0), lddbc,
                                c_one, dA(igpu, 0, k/nrgpu), ldda);

                    magma_sgetmatrix_async(k*nb, kb,
                                           dA(igpu, 0, k/nrgpu), ldda,
                                           A(0, k),              lda, stream[igpu][0] );
                }

                if(k>0){
                    magma_queue_sync( stream[igpu][0] );
                    // send the partially updated panel of dA to each gpu in the second dB block
                    // to overlap symm computation

                    for (igpu = 0; igpu < nrgpu; ++igpu){
                        magma_setdevice(igpu);
                        magma_ssetmatrix_async(k*nb, kb,
                                               A(0, k),          lda,
                                               dB_c(igpu, 0, 1), lddbc, stream[igpu][0] );
                    }

                    igpu = k%nrgpu;
                    magma_setdevice(igpu);
                    magmablasSetKernelStream(stream[igpu][2]);
                    magma_ssymm(MagmaRight, uplo,
                                k*nb, kb,
                                c_half, dA(igpu, k, k/nrgpu), ldda,
                                dB_c(igpu, 0, 0), lddbc,
                                c_one, dA(igpu, 0, k/nrgpu), ldda);

                    magma_strmm(MagmaRight, uplo, MagmaTrans, MagmaNonUnit,
                                k*nb, kb,
                                c_one, dB_c(igpu, k, 0), lddbc,
                                dA(igpu, 0, k/nrgpu), ldda);

                    for (igpu = 0; igpu < nrgpu; ++igpu){
                        magma_queue_sync( stream[igpu][0] );
                    }

                    for (j = 0; j < k; ++j){
                        igpu = j%nrgpu;
                        magma_setdevice(igpu);
                        magmablasSetKernelStream(stream[igpu][(j/nrgpu)%3]);
                        magma_ssyr2k(uplo, MagmaNoTrans,
                                     nb, kb,
                                     c_one, dB_c(igpu, j, 1), lddbc,
                                     dB_c(igpu, j, 0), lddbc,
                                     d_one, dA(igpu, j, j/nrgpu), ldda);
                        magma_queue_sync( stream[igpu][((j)/nrgpu)%3] ); // Needed for correctness. Why?
                    }

                    for (j = 1; j < k; ++j){
                        igpu = j%nrgpu;
                        magma_setdevice(igpu);
                        magmablasSetKernelStream(stream[igpu][0]);
                        magma_sgemm(MagmaNoTrans, MagmaTrans, j*nb, nb, kb, c_one, dB_c(igpu, 0, 0), lddbc,
                                    dB_c(igpu, j, 1), lddbc, c_one, dA(igpu, 0, j/nrgpu), ldda );

                        magma_sgemm(MagmaNoTrans, MagmaTrans, j*nb, nb, kb, c_one, dB_c(igpu, 0, 1), lddbc,
                                    dB_c(igpu, j, 0), lddbc, c_one, dA(igpu, 0, j/nrgpu), ldda );
                    }

                }
                if (k < nbl-1){
                    for (igpu = 0; igpu < nrgpu; ++igpu){
                        magma_setdevice(igpu);
                        magma_queue_sync( stream[igpu][0] );
                        magma_ssetmatrix_async((k+1)*nb+kb2, kb2,
                                               B(0, k+1),        ldb,
                                               dB_c(igpu, 0, 0), lddbc, stream[igpu][0] );
                    }
                }
                lapackf77_ssygs2( &itype, uplo_, &kb, A(k,k), &lda, B(k,k), &ldb, info);

                igpu = k%nrgpu;
                magma_setdevice(igpu);
                magma_ssetmatrix_async(kb, kb,
                                       A(k, k),              lda,
                                       dA(igpu, k, k/nrgpu), ldda, stream[igpu][0] );
            }
            for (igpu = 0; igpu < nrgpu; ++igpu){
                magma_queue_sync( stream[igpu][0] );
                magma_queue_sync( stream[igpu][1] );
                magma_queue_sync( stream[igpu][2] );
            }

            //copy A from mgpus
            for (j = 0; j < nbl; ++j){
                igpu = j%nrgpu;
                magma_setdevice(igpu);
                jb = min(nb, n-j*nb);
                magma_sgetmatrix_async(j*nb+jb, jb,
                                       dA(igpu, 0, j/nrgpu), ldda,
                                       A(0, j),              lda, stream[igpu][0] );
            }

        } else {
            /* Compute L'*A*L */

            if (n > nb){

                magma_int_t nloc[MagmaMaxGPUs];
                magma_int_t iloc[MagmaMaxGPUs];
                for(igpu = 0; igpu < nrgpu; ++igpu){
                    nloc[igpu] = 0;
                    iloc[igpu] = 0;
                }

                kb = min(nb, n);
                for (j = 0; j < nbl; ++j){
                    igpu = j%nrgpu;
                    magma_setdevice(igpu);
                    jb = min(nb, n-j*nb);
                    nloc[igpu] += jb;
                    magma_ssetmatrix_async(jb, kb,
                                           A(j, 0),              lda,
                                           dA(igpu, j/nrgpu, 0), ldda, stream[igpu][0] );
                }
                for (igpu = 0; igpu < nrgpu; ++igpu){
                    magma_setdevice(igpu);
                    magma_ssetmatrix_async(kb, kb,
                                           B(0, 0),          ldb,
                                           dB_r(igpu, 0, 0), lddbr, stream[igpu][0] );
                }
                for (k = 0; k < nbl-1; ++k){
                    ++iloc[k%nrgpu];
                    kb = min(nb, n-(k+1)*nb);
                    for (j = k+1; j < nbl; ++j){
                        igpu = j%nrgpu;
                        magma_setdevice(igpu);
                        jb = min(nb, n-j*nb);
                        magma_ssetmatrix_async(jb, kb,
                                               A(j, k+1),              lda,
                                               dA(igpu, j/nrgpu, k+1), ldda, stream[igpu][(k+1)%2] );
                    }
                    for (igpu = 0; igpu < nrgpu; ++igpu){
                        magma_setdevice(igpu);
                        magma_ssetmatrix_async(kb, (k+1)*nb + kb,
                                               B(k+1, 0),              ldb,
                                               dB_r(igpu, (k+1)%2, 0), lddbr, stream[igpu][(k+1)%2] );
                    }

                    kb = min(nb, n-k*nb);

                    if (k > 0){
                        for (igpu = 0; igpu < nrgpu; ++igpu){
                            magma_setdevice(igpu);
                            magmablasSetKernelStream(stream[igpu][k%2]);
                            magma_sgemm(MagmaNoTrans, MagmaNoTrans, nloc[igpu]-nb*iloc[igpu], k*nb, kb, c_one, dA(igpu, iloc[igpu], k), ldda,
                                        dB_r(igpu, k%2, 0), lddbr, c_one, dA(igpu, iloc[igpu], 0), ldda );
                        }
                    }
                    for (igpu = 0; igpu < nrgpu; ++igpu){
                        magma_setdevice(igpu);
                        magmablasSetKernelStream(stream[igpu][k%2]);
                        magma_strmm(MagmaRight, uplo, MagmaNoTrans, MagmaNonUnit, nloc[igpu]-nb*iloc[igpu], kb, c_one, dB_r(igpu, k%2, k), lddbr,
                                    dA(igpu, iloc[igpu], k), ldda );
                    }

                    for (igpu = 0; igpu < nrgpu; ++igpu){
                        magma_queue_sync( stream[igpu][k%2] );
                    }
                }
            }

            for(k = 0; k<nbl; ++k){
                kb= min(n-k*nb,nb);
                kb2= min(n-(k+1)*nb,nb);

                if(k>0){

                    igpu = k%nrgpu;
                    magma_setdevice(igpu);
                    magmablasSetKernelStream(stream[igpu][0]);
                    magma_ssymm(MagmaLeft, uplo,
                                kb, k*nb,
                                c_half, dA(igpu, k/nrgpu, k), ldda,
                                dB_r(igpu, 0, 0), lddbr,
                                c_one, dA(igpu, k/nrgpu, 0), ldda);

                    magma_sgetmatrix_async(kb, k*nb,
                                           dA(igpu, k/nrgpu, 0), ldda,
                                           A(k, 0),              lda, stream[igpu][0] );
                }

                if(k>0){
                    magma_queue_sync( stream[igpu][0] );
                    // send the partially updated panel of dA to each gpu in the second dB block
                    // to overlap symm computation

                    for (igpu = 0; igpu < nrgpu; ++igpu){
                        magma_setdevice(igpu);
                        magma_ssetmatrix_async(kb, k*nb,
                                               A(k, 0),          lda,
                                               dB_r(igpu, 1, 0), lddbr, stream[igpu][0] );
                    }

                    igpu = k%nrgpu;
                    magma_setdevice(igpu);
                    magmablasSetKernelStream(stream[igpu][2]);
                    magma_ssymm(MagmaLeft, uplo,
                                kb, k*nb,
                                c_half, dA(igpu, k/nrgpu, k), ldda,
                                dB_r(igpu, 0, 0), lddbr,
                                c_one, dA(igpu, k/nrgpu, 0), ldda);

                    magma_strmm(MagmaLeft, uplo, MagmaTrans, MagmaNonUnit,
                                kb, k*nb,
                                c_one, dB_r(igpu, 0, k), lddbr,
                                dA(igpu, k/nrgpu, 0), ldda);

                    for (igpu = 0; igpu < nrgpu; ++igpu){
                        magma_queue_sync( stream[igpu][0] );
                    }

                    for (j = 0; j < k; ++j){
                        igpu = j%nrgpu;
                        magma_setdevice(igpu);
                        magmablasSetKernelStream(stream[igpu][(j/nrgpu)%3]);
                        magma_ssyr2k(uplo, MagmaTrans,
                                     nb, kb,
                                     c_one, dB_r(igpu, 1, j), lddbr,
                                     dB_r(igpu, 0, j), lddbr,
                                     d_one, dA(igpu, j/nrgpu, j), ldda);
                        magma_queue_sync( stream[igpu][((j)/nrgpu)%3] ); // Needed for correctness. Why?
                    }

                    for (j = 1; j < k; ++j){
                        igpu = j%nrgpu;
                        magma_setdevice(igpu);
                        magmablasSetKernelStream(stream[igpu][0]);
                        magma_sgemm(MagmaTrans, MagmaNoTrans, nb, j*nb, kb, c_one, dB_r(igpu, 0, j), lddbr,
                                    dB_r(igpu, 1, 0), lddbr, c_one, dA(igpu, j/nrgpu, 0), ldda );

                        magma_sgemm(MagmaTrans, MagmaNoTrans, nb, j*nb, kb, c_one, dB_r(igpu, 1, j), lddbr,
                                    dB_r(igpu, 0, 0), lddbr, c_one, dA(igpu, j/nrgpu, 0), ldda );
                    }

                }
                if (k < nbl-1){
                    for (igpu = 0; igpu < nrgpu; ++igpu){
                        magma_setdevice(igpu);
                        magma_queue_sync( stream[igpu][0] );
                        magma_ssetmatrix_async(kb2, (k+1)*nb+kb2,
                                               B(k+1, 0),        ldb,
                                               dB_r(igpu, 0, 0), lddbr, stream[igpu][0] );
                    }
                }
                lapackf77_ssygs2( &itype, uplo_, &kb, A(k,k), &lda, B(k,k), &ldb, info);

                igpu = k%nrgpu;
                magma_setdevice(igpu);
                magma_ssetmatrix_async(kb, kb,
                                       A(k, k),              lda,
                                       dA(igpu, k/nrgpu, k), ldda, stream[igpu][0] );
            }
            for (igpu = 0; igpu < nrgpu; ++igpu){
                magma_queue_sync( stream[igpu][0] );
                magma_queue_sync( stream[igpu][1] );
                magma_queue_sync( stream[igpu][2] );
            }

            //copy A from mgpus
            for (j = 0; j < nbl; ++j){
                igpu = j%nrgpu;
                magma_setdevice(igpu);
                jb = min(nb, n-j*nb);
                magma_sgetmatrix_async(jb, j*nb+jb,
                                       dA(igpu, j/nrgpu, 0), ldda,
                                       A(j, 0),              lda, stream[igpu][0] );
            }

        }
    }

    for (igpu = 0; igpu < nrgpu; ++igpu){
        magma_setdevice(igpu);
        magmablasSetKernelStream(NULL);
        magma_queue_sync( stream[igpu][2] );
        magma_queue_destroy( stream[igpu][0] );
        magma_queue_destroy( stream[igpu][1] );
        magma_queue_destroy( stream[igpu][2] );
        magma_free( dw[igpu] );
    }

    magma_setdevice(gpu_b);

    return *info;
} /* magma_ssygst_gpu */
