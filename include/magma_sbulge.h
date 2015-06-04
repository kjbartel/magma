/*
 *   -- MAGMA (version 1.3.0) --
 *      Univ. of Tennessee, Knoxville
 *      Univ. of California, Berkeley
 *      Univ. of Colorado, Denver
 *      November 2012
 *
 * @generated s Wed Nov 14 22:52:27 2012
 */

#ifndef MAGMA_SBULGEINC_H
#define MAGMA_SBULGEINC_H

#define PRECISION_s

extern "C" {
    magma_int_t magma_ssytrd_sy2sb(char uplo, magma_int_t n, magma_int_t NB, float *a, magma_int_t lda,
                                   float *tau, float *work, magma_int_t lwork, float *dT, magma_int_t threads, magma_int_t *info);

    magma_int_t magma_ssytrd_hb2st(magma_int_t threads, char uplo, magma_int_t n, magma_int_t nb, magma_int_t Vblksiz,
                                   float *A, magma_int_t lda, float *D, float *E,
                                   float *V, magma_int_t ldv, float *TAU, magma_int_t compT, float *T, magma_int_t ldt);

    magma_int_t magma_sbulge_back(magma_int_t threads, char uplo, magma_int_t n, magma_int_t nb, magma_int_t ne, magma_int_t Vblksiz,
                                  float *Z, magma_int_t ldz, float *dZ, magma_int_t lddz,
                                  float *V, magma_int_t ldv, float *TAU, float *T, magma_int_t ldt, magma_int_t* info);

    magma_int_t magma_sormqr_gpu_2stages(char side, char trans, magma_int_t m, magma_int_t n, magma_int_t k, float *dA, magma_int_t ldda,
                                         float *dC, magma_int_t lddc, float *dT, magma_int_t nb, magma_int_t *info);

    magma_int_t magma_sbulge_get_lq2(magma_int_t n);

    magma_int_t magma_sbulge_get_Vblksiz(magma_int_t n, magma_int_t nb);
}
#endif
