/*
 *   -- MAGMA (version 1.3.0) --
 *      Univ. of Tennessee, Knoxville
 *      Univ. of California, Berkeley
 *      Univ. of Colorado, Denver
 *      November 2012
 *
 * @generated c Wed Nov 14 22:52:27 2012
 */

#ifndef MAGMA_CBULGEINC_H
#define MAGMA_CBULGEINC_H

#define PRECISION_c

extern "C" {
    magma_int_t magma_chetrd_he2hb(char uplo, magma_int_t n, magma_int_t NB, cuFloatComplex *a, magma_int_t lda,
                                   cuFloatComplex *tau, cuFloatComplex *work, magma_int_t lwork, cuFloatComplex *dT, magma_int_t threads, magma_int_t *info);

    magma_int_t magma_chetrd_hb2st(magma_int_t threads, char uplo, magma_int_t n, magma_int_t nb, magma_int_t Vblksiz,
                                   cuFloatComplex *A, magma_int_t lda, float *D, float *E,
                                   cuFloatComplex *V, magma_int_t ldv, cuFloatComplex *TAU, magma_int_t compT, cuFloatComplex *T, magma_int_t ldt);

    magma_int_t magma_cbulge_back(magma_int_t threads, char uplo, magma_int_t n, magma_int_t nb, magma_int_t ne, magma_int_t Vblksiz,
                                  cuFloatComplex *Z, magma_int_t ldz, cuFloatComplex *dZ, magma_int_t lddz,
                                  cuFloatComplex *V, magma_int_t ldv, cuFloatComplex *TAU, cuFloatComplex *T, magma_int_t ldt, magma_int_t* info);

    magma_int_t magma_cunmqr_gpu_2stages(char side, char trans, magma_int_t m, magma_int_t n, magma_int_t k, cuFloatComplex *dA, magma_int_t ldda,
                                         cuFloatComplex *dC, magma_int_t lddc, cuFloatComplex *dT, magma_int_t nb, magma_int_t *info);

    magma_int_t magma_cbulge_get_lq2(magma_int_t n);

    magma_int_t magma_cbulge_get_Vblksiz(magma_int_t n, magma_int_t nb);
}
#endif
