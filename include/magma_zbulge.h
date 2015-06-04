/*
 *   -- MAGMA (version 1.3.0) --
 *      Univ. of Tennessee, Knoxville
 *      Univ. of California, Berkeley
 *      Univ. of Colorado, Denver
 *      November 2012
 *
 * @precisions normal z -> s d c
 */

#ifndef MAGMA_ZBULGEINC_H
#define MAGMA_ZBULGEINC_H

#define PRECISION_z

extern "C" {
    magma_int_t magma_zhetrd_he2hb(char uplo, magma_int_t n, magma_int_t NB, cuDoubleComplex *a, magma_int_t lda,
                                   cuDoubleComplex *tau, cuDoubleComplex *work, magma_int_t lwork, cuDoubleComplex *dT, magma_int_t threads, magma_int_t *info);

    magma_int_t magma_zhetrd_hb2st(magma_int_t threads, char uplo, magma_int_t n, magma_int_t nb, magma_int_t Vblksiz,
                                   cuDoubleComplex *A, magma_int_t lda, double *D, double *E,
                                   cuDoubleComplex *V, magma_int_t ldv, cuDoubleComplex *TAU, magma_int_t compT, cuDoubleComplex *T, magma_int_t ldt);

    magma_int_t magma_zbulge_back(magma_int_t threads, char uplo, magma_int_t n, magma_int_t nb, magma_int_t ne, magma_int_t Vblksiz,
                                  cuDoubleComplex *Z, magma_int_t ldz, cuDoubleComplex *dZ, magma_int_t lddz,
                                  cuDoubleComplex *V, magma_int_t ldv, cuDoubleComplex *TAU, cuDoubleComplex *T, magma_int_t ldt, magma_int_t* info);

    magma_int_t magma_zunmqr_gpu_2stages(char side, char trans, magma_int_t m, magma_int_t n, magma_int_t k, cuDoubleComplex *dA, magma_int_t ldda,
                                         cuDoubleComplex *dC, magma_int_t lddc, cuDoubleComplex *dT, magma_int_t nb, magma_int_t *info);

    magma_int_t magma_zbulge_get_lq2(magma_int_t n);

    magma_int_t magma_zbulge_get_Vblksiz(magma_int_t n, magma_int_t nb);
}
#endif
