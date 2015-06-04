/*
    -- MAGMA (version 1.2.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       May 2012

       @author Raffaele Solca

       @generated d Tue May 15 18:17:49 2012
*/
#define N_MAX_GPU 8

#include "common_magma.h"
#include "cblas.h"

extern "C"
magma_int_t magma_get_dsygst_m_nb() { return 256;}

#define A(i, j) (a+(j)*nb*lda + (i)*nb)
#define B(i, j) (b+(j)*nb*ldb + (i)*nb)

#define dA(gpui, i, j) (dw[gpui] + (j)*nb*ldda + (i)*nb)
#define dB_c(gpui, i, j) (dw[gpui] + dima*ldda + (i)*nb + (j)*nb*lddbc)
#define dB_r(gpui, i, j) (dw[gpui] + dima*ldda + (i)*nb + (j)*nb*lddbr)

extern "C" magma_int_t
magma_dsygst_m(magma_int_t nrgpu, magma_int_t itype, char uplo, magma_int_t n,
               double *a, magma_int_t lda,
               double *b, magma_int_t ldb, magma_int_t *info)
{
/*
  -- MAGMA (version 1.2.0) --
     Univ. of Tennessee, Knoxville
     Univ. of California, Berkeley
     Univ. of Colorado, Denver
     May 2012


   Purpose
   =======

   DSYGST reduces a real symmetric-definite generalized
   eigenproblem to standard form.

   If ITYPE = 1, the problem is A*x = lambda*B*x,
   and A is overwritten by inv(U**T)*A*inv(U) or inv(L)*A*inv(L**T)

   If ITYPE = 2 or 3, the problem is A*B*x = lambda*x or
   B*A*x = lambda*x, and A is overwritten by U*A*U**T or L**T*A*L.

   B must have been previously factorized as U**T*U or L*L**T by DPOTRF.

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
           as returned by DPOTRF.

   LDB     (input) INTEGER
           The leading dimension of the array B.  LDB >= max(1,N).

   INFO    (output) INTEGER
           = 0:  successful exit
           < 0:  if INFO = -i, the i-th argument had an illegal value

   =====================================================================*/

    char uplo_[2] = {uplo, 0};
    magma_int_t        k, kb, j, jb, kb2;
    magma_int_t        ldda, dima, lddbr, lddbc;
    double    c_one      = MAGMA_D_ONE;
    double    c_neg_one  = MAGMA_D_NEG_ONE;
    double    c_half     = MAGMA_D_HALF;
    double    c_neg_half = MAGMA_D_NEG_HALF;
    double* dw[N_MAX_GPU];
    cudaStream_t stream [N_MAX_GPU][3];
    magma_int_t igpu = 0;

    int gpu_b;
    magma_getdevice(&gpu_b);

    double             d_one = 1.0;
    long int           upper = lapackf77_lsame(uplo_, "U");

    magma_int_t nb = magma_get_dsygst_m_nb();

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
        if (MAGMA_SUCCESS != magma_dmalloc( &dw[igpu], (dima*ldda + lddbc*lddbr) )) {
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
                magma_dsetmatrix_async( kb, n-k*nb,
                                        A(k, k),              lda,
                                        dA(igpu, k/nrgpu, k), ldda, stream[igpu][0] );
            }
            kb= min(n,nb);
            igpu = 0;
            magma_setdevice(igpu);
            // dB_r(0,0) is used to store B(k,k)
            magma_dsetmatrix_async( kb, kb,
                                    B(0, 0),          ldb,
                                    dB_r(igpu, 0, 0), lddbr, stream[igpu][1] );

            for(k = 0; k<nbl; ++k){
                kb= min(n-k*nb,nb);
                kb2= min(n-(k+1)*nb,nb);

                if(k+1<nbl){
                    for (igpu = 0; igpu < nrgpu; ++igpu){
                        magma_setdevice(igpu);
                        magma_queue_sync( stream[igpu][0] );
                        magma_dsetmatrix_async( kb, n-(k+1)*nb,
                                                B(k, k+1),          ldb,
                                                dB_r(igpu, 0, k+1), lddbr, stream[igpu][0] );
                    }
                }

                igpu = k%nrgpu;
                magma_setdevice(igpu);

                magma_queue_sync( stream[igpu][1] ); // Needed, otherwise conflicts reading B(k,k) between hegs2 and cudaMemcpy2D
                magma_queue_sync( stream[igpu][2] );

                if(k+1<nbl){
                    magmablasSetKernelStream(stream[igpu][1]);
                    // dB_r(0,0) stores B(k,k)
                    magma_dtrsm(MagmaLeft, MagmaUpper, MagmaTrans, MagmaNonUnit,
                                kb, n-(k+1)*nb,
                                c_one, dB_r(igpu, 0, 0), lddbr,
                                dA(igpu, k/nrgpu, k+1), ldda);
                }

                lapackf77_dhegs2( &itype, uplo_, &kb, A(k,k), &lda, B(k,k), &ldb, info);
printf("hegs2%d\n", k);
                if (k+1<nbl) {
                    magma_dsetmatrix_async( kb, kb,
                                            A(k, k),              lda,
                                            dA(igpu, k/nrgpu, k), ldda, stream[igpu][0] );

                    magma_queue_sync( stream[igpu][1] );
                    magmablasSetKernelStream(stream[igpu][0]);

                    magma_dsymm(MagmaLeft, MagmaUpper,
                                kb, n-(k+1)*nb,
                                c_neg_half, dA(igpu, k/nrgpu, k), ldda,
                                dB_r(igpu, 0, k+1), lddbr,
                                c_one, dA(igpu, k/nrgpu, k+1), ldda);

                    magma_queue_sync( stream[igpu][0] );

                    magma_dgetmatrix( kb, n-(k+1)*nb,
                                      dA(igpu, k/nrgpu, k+1), ldda,
                                      A(k, k+1),              lda );

                    // send the partially updated panel of dA to each gpu in the second dB block
                    // to overlap symm computation

                    for (igpu = 0; igpu < nrgpu; ++igpu){
                        magma_setdevice(igpu);
                        magma_dsetmatrix_async( kb, n-(k+1)*nb,
                                                A(k, k+1),          lda,
                                                dB_r(igpu, 1, k+1), lddbr, stream[igpu][0] );
                    }

                    igpu = k%nrgpu;
                    magma_setdevice(igpu);
                    magmablasSetKernelStream(stream[igpu][1]);

                    magma_dsymm(MagmaLeft, MagmaUpper,
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
                        magma_dsyr2k(MagmaUpper, MagmaTrans,
                                     jb, nb,
                                     c_neg_one, dB_r(igpu, 1, j), lddbr,
                                     dB_r(igpu, 0, j), lddbr,
                                     d_one, dA(igpu, j/nrgpu, j), ldda);
                        magma_queue_sync( stream[igpu][((j)/nrgpu)%3] ); // Needed for correctness. Why?
                        if (j == k+1){
                            magma_queue_sync( stream[igpu][(j/nrgpu)%3] );
                            magma_dgetmatrix_async( kb2, kb2,
                                                    dA(igpu, (k+1)/nrgpu, k+1), ldda,
                                                    A(k+1, k+1),                lda, stream[igpu][2] );
                            // dB_r(0,0) is used to store B(k,k)
                            magma_dsetmatrix_async( kb2, kb2,
                                                    B(k+1, k+1),      ldb,
                                                    dB_r(igpu, 0, 0), lddbr, stream[igpu][1] );
                        }
                    }
                    for (j = k+1; j < nbl-1; ++j){
                        igpu = j%nrgpu;
                        magma_setdevice(igpu);
                        magmablasSetKernelStream(stream[igpu][0]);
                        magma_dgemm(MagmaTrans, MagmaNoTrans, nb, n-(j+1)*nb, nb, c_neg_one, dB_r(igpu, 0, j), lddbr,
                                    dB_r(igpu, 1, j+1), lddbr, c_one, dA(igpu, j/nrgpu, j+1), ldda );

                        magma_dgemm(MagmaTrans, MagmaNoTrans, nb, n-(j+1)*nb, nb, c_neg_one, dB_r(igpu, 1, j), lddbr,
                                    dB_r(igpu, 0, j+1), lddbr, c_one, dA(igpu, j/nrgpu, j+1), ldda );
                    }
                }
            }

            for (igpu = 0; igpu < nrgpu; ++igpu){
                magma_queue_sync( stream[igpu][0] );
                magma_queue_sync( stream[igpu][1] );
            }

            if (n > nb){

                magma_int_t nloc[N_MAX_GPU];

                jb = min(nb, n-nb);
                for (igpu = 0; igpu < nrgpu; ++igpu){
                    nloc[igpu]=0;
                    magma_setdevice(igpu);
                    magma_dsetmatrix_async( jb, n-nb,
                                            B(1, 1),          ldb,
                                            dB_r(igpu, 1, 1), lddbr, stream[igpu][1] );
                }
                for (j = 1; j < nbl; ++j){
                    if ((j+1)*nb < n){
                        jb = min(nb, n-(j+1)*nb);
                        for (igpu = 0; igpu < nrgpu; ++igpu){
                            magma_setdevice(igpu);
                            magma_dsetmatrix_async( jb, n-(j+1)*nb,
                                                    B(j+1, j+1),              ldb,
                                                    dB_r(igpu, (j+1)%2, j+1), lddbr, stream[igpu][(j+1)%2] );
                        }
                    }
                    jb = min(nb, n-j*nb);
                    nloc[(j-1)%nrgpu] += nb;

                    for (igpu = 0; igpu < nrgpu; ++igpu){
                        magma_setdevice(igpu);
                        magmablasSetKernelStream(stream[igpu][j%2]);
                        magma_dtrsm(MagmaRight, uplo, MagmaNoTrans, MagmaNonUnit, nloc[igpu], jb, c_one, dB_r(igpu, j%2, j), lddbr,
                                    dA(igpu, 0, j), ldda );
                    }

                    if ( j < nbl-1 ){

                        for (igpu = 0; igpu < nrgpu; ++igpu){
                            magma_setdevice(igpu);
                            magmablasSetKernelStream(stream[igpu][j%2]);
                            magma_dgemm(MagmaNoTrans, MagmaNoTrans, nloc[igpu], n-(j+1)*nb, nb, c_neg_one, dA(igpu, 0, j), ldda,
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
                        magma_dgetmatrix_async( kb, jb,
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
                magma_dsetmatrix_async( (n-k*nb), kb,
                                        A(k, k),              lda,
                                        dA(igpu, k, k/nrgpu), ldda, stream[igpu][0] );
            }
            kb= min(n,nb);
            igpu = 0;
            magma_setdevice(igpu);
            // dB_c(0,0) is used to store B(k,k)
            magma_dsetmatrix_async( kb, kb,
                                    B(0, 0),          ldb,
                                    dB_c(igpu, 0, 0), lddbc, stream[igpu][1] );

            for(k = 0; k<nbl; ++k){
                kb= min(n-k*nb,nb);
                kb2= min(n-(k+1)*nb,nb);

                if(k+1<nbl){
                    for (igpu = 0; igpu < nrgpu; ++igpu){
                        magma_setdevice(igpu);
                        magma_queue_sync( stream[igpu][0] );
                        magma_dsetmatrix_async( (n-(k+1)*nb), kb,
                                                B(k+1, k),          ldb,
                                                dB_c(igpu, k+1, 0), lddbc, stream[igpu][0] );
                    }
                }

                igpu = k%nrgpu;
                magma_setdevice(igpu);

                magma_queue_sync( stream[igpu][1] ); // Needed, otherwise conflicts reading B(k,k) between hegs2 and cudaMemcpy2D
                magma_queue_sync( stream[igpu][2] );

                if(k+1<nbl){
                    magmablasSetKernelStream(stream[igpu][1]);
                    // dB_c(0,0) stores B(k,k)
                    magma_dtrsm(MagmaRight, MagmaLower, MagmaTrans, MagmaNonUnit,
                                n-(k+1)*nb, kb,
                                c_one, dB_c(igpu, 0, 0), lddbc,
                                dA(igpu, k+1, k/nrgpu), ldda);
                }

                lapackf77_dhegs2( &itype, uplo_, &kb, A(k,k), &lda, B(k,k), &ldb, info);

                if (k+1<nbl) {
                    magma_dsetmatrix_async( kb, kb,
                                            A(k, k),               lda,
                                            dA(igpu, k , k/nrgpu), ldda, stream[igpu][0] );

                    magma_queue_sync( stream[igpu][1] );
                    magmablasSetKernelStream(stream[igpu][0]);

                    magma_dsymm(MagmaRight, MagmaLower,
                                n-(k+1)*nb, kb,
                                c_neg_half, dA(igpu, k, k/nrgpu), ldda,
                                dB_c(igpu, k+1, 0), lddbc,
                                c_one, dA(igpu, k+1, k/nrgpu), ldda);

                    magma_queue_sync( stream[igpu][0] );

                    magma_dgetmatrix( n-(k+1)*nb, kb,
                                      dA(igpu, k+1, k/nrgpu), ldda,
                                      A(k+1, k),              lda );

                    // send the partially updated panel of dA to each gpu in the second dB block
                    // to overlap symm computation

                    for (igpu = 0; igpu < nrgpu; ++igpu){
                        magma_setdevice(igpu);
                        magma_dsetmatrix_async( (n-(k+1)*nb), kb,
                                                A(k+1, k),          lda,
                                                dB_c(igpu, k+1, 1), lddbc, stream[igpu][0] );
                    }

                    igpu = k%nrgpu;
                    magma_setdevice(igpu);
                    magmablasSetKernelStream(stream[igpu][1]);

                    magma_dsymm(MagmaRight, MagmaLower,
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
                        magma_dsyr2k(MagmaLower, MagmaNoTrans,
                                     jb, nb,
                                     c_neg_one, dB_c(igpu, j, 1), lddbc,
                                     dB_c(igpu, j, 0), lddbc,
                                     d_one, dA(igpu, j, j/nrgpu), ldda);
                        magma_queue_sync( stream[igpu][((j)/nrgpu)%3] ); // Needed for correctness. Why?
                        if (j == k+1){
                            magma_queue_sync( stream[igpu][(j/nrgpu)%3] );
                            magma_dgetmatrix_async( kb2, kb2,
                                                    dA(igpu, k+1, (k+1)/nrgpu), ldda,
                                                    A(k+1, k+1),                lda, stream[igpu][2] );
                            // dB_c(0,0) is used to store B(k,k)
                            magma_dsetmatrix_async( kb2, kb2,
                                                    B(k+1, k+1),      ldb,
                                                    dB_c(igpu, 0, 0), lddbc, stream[igpu][1] );
                        }
                    }
                    for (j = k+1; j < nbl-1; ++j){
                        igpu = j%nrgpu;
                        magma_setdevice(igpu);
                        magmablasSetKernelStream(stream[igpu][0]);
                        magma_dgemm(MagmaNoTrans, MagmaTrans, n-(j+1)*nb, nb, nb, c_neg_one, dB_c(igpu, j+1, 1), lddbc,
                                    dB_c(igpu, j, 0), lddbc, c_one, dA(igpu, j+1, j/nrgpu), ldda );

                        magma_dgemm(MagmaNoTrans, MagmaTrans, n-(j+1)*nb, nb, nb, c_neg_one, dB_c(igpu, j+1, 0), lddbc,
                                    dB_c(igpu, j, 1), lddbc, c_one, dA(igpu, j+1, j/nrgpu), ldda );
                    }
                }
            }

            for (igpu = 0; igpu < nrgpu; ++igpu){
                magma_queue_sync( stream[igpu][0] );
                magma_queue_sync( stream[igpu][1] );
            }

            if (n > nb){

                magma_int_t nloc[N_MAX_GPU];

                jb = min(nb, n-nb);
                for (igpu = 0; igpu < nrgpu; ++igpu){
                    nloc[igpu]=0;
                    magma_setdevice(igpu);
                    magma_dsetmatrix_async( (n-nb), jb,
                                            B(1, 1),          ldb,
                                            dB_c(igpu, 1, 1), lddbc, stream[igpu][1] );
                }
                for (j = 1; j < nbl; ++j){
                    if ((j+1)*nb < n){
                        jb = min(nb, n-(j+1)*nb);
                        for (igpu = 0; igpu < nrgpu; ++igpu){
                            magma_setdevice(igpu);
                            magma_dsetmatrix_async( (n-(j+1)*nb), jb,
                                                    B(j+1, j+1),              ldb,
                                                    dB_c(igpu, j+1, (j+1)%2), lddbc, stream[igpu][(j+1)%2] );
                        }
                    }
                    jb = min(nb, n-j*nb);
                    nloc[(j-1)%nrgpu] += nb;

                    for (igpu = 0; igpu < nrgpu; ++igpu){
                        magma_setdevice(igpu);
                        magmablasSetKernelStream(stream[igpu][j%2]);
                        magma_dtrsm(MagmaLeft, uplo, MagmaNoTrans, MagmaNonUnit, jb, nloc[igpu], c_one, dB_c(igpu, j, j%2), lddbc,
                                    dA(igpu, j, 0), ldda );
                    }

                    if ( j < nbl-1 ){

                        for (igpu = 0; igpu < nrgpu; ++igpu){
                            magma_setdevice(igpu);
                            magmablasSetKernelStream(stream[igpu][j%2]);
                            magma_dgemm(MagmaNoTrans, MagmaNoTrans, n-(j+1)*nb, nloc[igpu], nb, c_neg_one, dB_c(igpu, j+1, j%2), lddbc,
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
                        magma_dgetmatrix_async( jb, kb,
                                                dA(igpu, j, k/nrgpu), ldda,
                                                A(j, k),              lda, stream[igpu][2] );
                    }
                }
            }

        }

    } else {

        if (upper) {

            printf("dsygst_m: type2 upper not implemented\n");
            exit(-1);

            /* Compute U*A*U' */

/*            for(k = 0; k<n; k+=nb){
                kb= min(n-k,nb);

                magma_dgetmatrix_async( kb, kb,
                                        dA(k, k), ldda,
                                        A(k, k),  lda, stream[0] );

                // Update the upper triangle of A(1:k+kb-1,1:k+kb-1)
                if(k>0){

                    magma_dtrmm(MagmaLeft, MagmaUpper, MagmaNoTrans, MagmaNonUnit,
                                k, kb,
                                c_one ,dB(0,0), lddb,
                                dA(0,k), ldda);

                    magma_dsymm(MagmaRight, MagmaUpper,
                                k, kb,
                                c_half, dA(k,k), ldda,
                                dB(0,k), lddb,
                                c_one, dA(0, k), ldda);

                    magma_queue_sync( stream[1] );

                    magma_dsyr2k(MagmaUpper, MagmaNoTrans,
                                 k, kb,
                                 c_one, dA(0,k), ldda,
                                 dB(0,k), lddb,
                                 d_one, dA(0,0), ldda);

                    magma_dsymm(MagmaRight, MagmaUpper,
                                k, kb,
                                c_half, dA(k,k), ldda,
                                dB(0,k), lddb,
                                c_one, dA(0, k), ldda);

                    magma_dtrmm(MagmaRight, MagmaUpper, MagmaTrans, MagmaNonUnit,
                                k, kb,
                                c_one, dB(k,k), lddb,
                                dA(0,k), ldda);

                }

                magma_queue_sync( stream[0] );

                lapackf77_dhegs2( &itype, uplo_, &kb, A(k, k), &lda, B(k, k), &ldb, info);

                magma_dsetmatrix_async( kb, kb,
                                        A(k, k),  lda,
                                        dA(k, k), ldda, stream[1] );

            }

            magma_queue_sync( stream[1] );
*/
        } else {

            /* Compute L'*A*L */

            printf("dsygst_m: type2 lower not implemented\n");
            exit(-1);

/*                        if (n > nb){

            magma_int_t nloc[N_MAX_GPU];
            for(igpu = 0; igpu < nrgpu; ++igpu)
                nloc[igpu] = 0;

            kb = min(nb, n);
            for (j = 0; j < nbl; ++j){
                igpu = j%nrgpu;
                magma_setdevice(igpu);
                jb = min(nb, n-j*nb);
                nloc[igpu] += jb;
                magma_dsetmatrix_async( jb, kb,
                                        A(j, 0),              lda,
                                        dA(igpu, j/nrgpu, 0), ldda, stream[igpu][0] );
            }
            for (igpu = 0; igpu < nrgpu; ++igpu){
                magma_setdevice(igpu);
                magma_dsetmatrix_async( kb, kb,
                                        B(0, 0),          ldb,
                                        dB_r(igpu, 0, 0), lddbr, stream[igpu][0] );
            }
            for (k = 0; k < nbl-1; ++k){
                nloc[k%nrgpu] -= nb;
                if (k < nbl-2){
                    kb = min(nb, n-(k+1)*nb);
                    for (j = k; j < nbl; ++j){
                        igpu = j%nrgpu;
                        magma_setdevice(igpu);
                        jb = min(nb, n-j*nb);
                        magma_dsetmatrix_async( jb, kb,
                                                A(j, k+1),              lda,
                                                dA(igpu, j/nrgpu, k+1), ldda, stream[igpu][(k+1)%2] );
                    }
                    for (igpu = 0; igpu < nrgpu; ++igpu){
                        magma_setdevice(igpu);
                        magma_dsetmatrix_async( kb, (k+1)*nb + kb,
                                                B(k+1, 0),              ldb,
                                                dB_r(igpu, (k+1)%2, 0), lddbr, stream[igpu][(k+1)%2] );
                    }
                }

                kb = min(nb, n-k*nb);

                if (k > 0){
                    for (igpu = 0; igpu < nrgpu; ++igpu){
                        magma_setdevice(igpu);
                        magmablasSetKernelStream(stream[igpu][k%2]);
                        magma_dgemm(MagmaNoTrans, MagmaNoTrans, nloc[igpu], k*nb, kb, c_one, dA(igpu, n-nloc[igpu], k), ldda,
                                    dB_r(igpu, k%2, 0), lddbr, c_one, dA(igpu, n-nloc[igpu], 0), ldda );
                    }
                }

                for (igpu = 0; igpu < nrgpu; ++igpu){
                    magma_setdevice(igpu);
                    magmablasSetKernelStream(stream[igpu][k%2]);
                    magma_dtrmm(MagmaRight, uplo, MagmaNoTrans, MagmaNonUnit, nloc[igpu], kb, c_one, dB_r(igpu, k%2, k), lddbr,
                                dA(igpu, n-nloc[igpu], k), ldda );
                }

                for (igpu = 0; igpu < nrgpu; ++igpu){
                    magma_queue_sync( stream[igpu][k%2] );
                }

            }

                        }

            /////////
            // put for loop!
            // put copies!

            magmablasSetKernelStream(stream[igpu][0]);

            magma_dsymm(MagmaRight, MagmaLower,
                        kb, k*nb,
                        c_half, dA(igpu, k/nrgpu, 0), ldda,
                        dB_r(igpu, 0, 0), lddbr,
                        c_one, dA(igpu, k/nrgpu, 0), ldda);

            magma_queue_sync( stream[igpu][0] );

            magma_dgetmatrix( kb, k*nb,
                              dA(igpu, k/nrgpu, 0), ldda,
                              A(k, 0),              lda );

            // send the partially updated panel of dA to each gpu in the second dB block
            // to overlap symm computation

            for (igpu = 0; igpu < nrgpu; ++igpu){
                magma_setdevice(igpu);
                magma_dsetmatrix_async( kb, // ERROR: missing dimension,
                                        A(k, 0),          lda,
                                        dB_r(igpu, 1, 0), lddbr, stream[igpu][0] );
            }

            igpu = k%nrgpu;
            magma_setdevice(igpu);
            magmablasSetKernelStream(stream[igpu][1]);

            magma_dsymm(MagmaRight, MagmaLower,
                        n-(k+1)*nb, kb,
                        c_neg_half, dA(igpu, k, k/nrgpu), ldda,
                        dB_c(igpu, k+1, 0), lddbc,
                        c_one, dA(igpu, k+1, k/nrgpu), ldda);

            for (igpu = 0; igpu < nrgpu; ++igpu){
                magma_queue_sync( stream[igpu][0] );
            }


            //copy B from mgpus
            for (j = 0; j < nbl; ++j){
                igpu = j%nrgpu;
                magma_setdevice(igpu);
                jb = min(nb, n-j*nb);
                magma_dgetmatrix_async( jb, n,
                                        dA(igpu, j/nrgpu, 0), ldda,
                                        A(j, 0),              lda, stream[igpu][0] );
            }

*//*            for(k = 0; k<n; k+=nb){
                kb= min(n-k,nb);

                magma_dgetmatrix_async( kb, kb,
                                        dA(k, k), ldda,
                                        A(k, k),  lda, stream[0] );

                // Update the lower triangle of A(1:k+kb-1,1:k+kb-1)
                if(k>0){

                    magma_dtrmm(MagmaRight, MagmaLower, MagmaNoTrans, MagmaNonUnit,
                                kb, k,
                                c_one ,dB(0,0), lddb,
                                dA(k,0), ldda);

                    magma_dsymm(MagmaLeft, MagmaLower,
                                kb, k,
                                c_half, dA(k,k), ldda,
                                dB(k,0), lddb,
                                c_one, dA(k, 0), ldda);

                    magma_queue_sync( stream[1] );

                    magma_dsyr2k(MagmaLower, MagmaTrans,
                                 k, kb,
                                 c_one, dA(k,0), ldda,
                                 dB(k,0), lddb,
                                 d_one, dA(0,0), ldda);

                    magma_dsymm(MagmaLeft, MagmaLower,
                                kb, k,
                                c_half, dA(k,k), ldda,
                                dB(k,0), lddb,
                                c_one, dA(k, 0), ldda);

                    magma_dtrmm(MagmaLeft, MagmaLower, MagmaTrans, MagmaNonUnit,
                                kb, k,
                                c_one, dB(k,k), lddb,
                                dA(k,0), ldda);
                }

                magma_queue_sync( stream[0] );

                lapackf77_dhegs2( &itype, uplo_, &kb, A(k,k), &lda, B(k,k), &ldb, info);

                magma_dsetmatrix_async( kb, kb,
                                        A(k, k),  lda,
                                        dA(k, k), ldda, stream[1] );
            }

            magma_queue_sync( stream[1] );
 */
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
} /* magma_dsygst_gpu */
