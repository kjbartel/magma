/*
 *   -- MAGMA (version 1.3.0) --
 *      Univ. of Tennessee, Knoxville
 *      Univ. of California, Berkeley
 *      Univ. of Colorado, Denver
 *      November 2012
 *
 * @generated d Wed Nov 14 22:52:27 2012
 */

#ifndef MAGMA_DBULGEINC_H
#define MAGMA_DBULGEINC_H

#define PRECISION_d

extern "C" {
    magma_int_t magma_dsytrd_sy2sb(char uplo, magma_int_t n, magma_int_t NB, double *a, magma_int_t lda,
                                   double *tau, double *work, magma_int_t lwork, double *dT, magma_int_t threads, magma_int_t *info);

    magma_int_t magma_dsytrd_hb2st(magma_int_t threads, char uplo, magma_int_t n, magma_int_t nb, magma_int_t Vblksiz,
                                   double *A, magma_int_t lda, double *D, double *E,
                                   double *V, magma_int_t ldv, double *TAU, magma_int_t compT, double *T, magma_int_t ldt);

    magma_int_t magma_dbulge_back(magma_int_t threads, char uplo, magma_int_t n, magma_int_t nb, magma_int_t ne, magma_int_t Vblksiz,
                                  double *Z, magma_int_t ldz, double *dZ, magma_int_t lddz,
                                  double *V, magma_int_t ldv, double *TAU, double *T, magma_int_t ldt, magma_int_t* info);

    magma_int_t magma_dormqr_gpu_2stages(char side, char trans, magma_int_t m, magma_int_t n, magma_int_t k, double *dA, magma_int_t ldda,
                                         double *dC, magma_int_t lddc, double *dT, magma_int_t nb, magma_int_t *info);

    magma_int_t magma_dbulge_get_lq2(magma_int_t n);

    magma_int_t magma_dbulge_get_Vblksiz(magma_int_t n, magma_int_t nb);
}
#endif
