/**
 *
 * @file compute_z.h
 *
 *  MAGMA auxiliary routines
 *  MAGMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.3.1
 * @author Jakub Kurzak
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 **/

#define PRECISION_z

/***************************************************************************//**
 *  Macro for matrix conversion / Lapack interface
 **/
#define magma_zdesc_alloc( descA, mb, nb, lm, ln, i, j, m, n, free)     \
    descA = magma_desc_init(                                            \
        PlasmaComplexDouble, (mb), (nb), ((mb)*(nb)),                   \
        (m), (n), (i), (j), (m), (n));                                  \
    if ( magma_desc_mat_alloc( &(descA) ) ) {                           \
        magma_error( __func__, "magma_shared_alloc() failed");          \
        {free;};                                                        \
        return MAGMA_ERR_OUT_OF_RESOURCES;                              \
    }

#define magma_zdesc_alloc2( descA, mb, nb, lm, ln, i, j, m, n)          \
    descA = magma_desc_init(                                            \
        PlasmaComplexDouble, (mb), (nb), ((mb)*(nb)),                   \
        (m), (n), (i), (j), (m), (n));                                  \
    magma_desc_mat_alloc( &(descA) );

#define magma_zooplap2tile( descA, A, mb, nb, lm, ln, i, j, m, n, free) \
    descA = magma_desc_init(                                            \
        PlasmaComplexDouble, (mb), (nb), ((mb)*(nb)),                   \
        (lm), (ln), (i), (j), (m), (n));                                \
    if ( magma_desc_mat_alloc( &(descA) ) ) {                           \
        magma_error( __func__, "magma_shared_alloc() failed");          \
        {free;};                                                        \
        return MAGMA_ERR_OUT_OF_RESOURCES;                              \
    }                                                                   \
    magma_pzlapack_to_tile((A), (lm), &(descA), sequence, &request);
    
#define magma_ziplap2tile( descA, A, mb, nb, lm, ln, i, j, m, n)      \
    descA = magma_desc_init(                                          \
        PlasmaComplexDouble, (mb), (nb), ((mb)*(nb)),                 \
        (lm), (ln), (i), (j), (m), (n));                              \
    descA.desc.mat = A;                                               \
    MAGMA_zgecfi_Async((lm), (ln), (A), MagmaCM, (mb), (nb),          \
                       MagmaCCRB, (mb), (nb), sequence, &request);

#define magma_zooptile2lap( descA, A, mb, nb, lm, ln)  \
    magma_pztile_to_lapack( &(descA), (A), (lm),       \
                            sequence, &request);

#define magma_ziptile2lap( descA, A, mb, nb, lm, ln)                  \
    MAGMA_zgecfi_Async((lm), (ln), (A), PlasmaCCRB, (mb), (nb),       \
                       PlasmaCM, (mb), (nb), sequence, &request);

/***************************************************************************//**
 *  Declarations of parallel functions (dynamic scheduling) - alphabetical order
 **/
void magma_pzlapack_to_tile(PLASMA_Complex64_t *Af77, int lda, magma_desc_t *A,
                            magma_sequence_t *sequence, magma_request_t *request);
void magma_pztile_to_lapack(magma_desc_t *A, PLASMA_Complex64_t *Af77, int lda,
                            magma_sequence_t *sequence, magma_request_t *request);
void magma_pztile_zero(magma_desc_t *A, magma_sequence_t *sequence, magma_request_t *request);


void magma_pzaxpy(PLASMA_Complex64_t alpha, magma_desc_t *A, magma_desc_t *B, magma_sequence_t *sequence, magma_request_t *request);
void magma_pzgelqf(magma_desc_t *A, magma_desc_t *T, magma_sequence_t *sequence, magma_request_t *request);
void magma_pzgelqfrh(magma_desc_t *A, magma_desc_t *T, int BS, magma_sequence_t *sequence, magma_request_t *request);
void magma_pzgemm(MAGMA_enum transA, MAGMA_enum transB, PLASMA_Complex64_t alpha, magma_desc_t *A, magma_desc_t *B, PLASMA_Complex64_t beta, magma_desc_t *C, magma_sequence_t *sequence, magma_request_t *request);
void magma_pzgeqrf(magma_desc_t *A, magma_desc_t *T, magma_sequence_t *sequence, magma_request_t *request);
void magma_pzgeqrfrh(magma_desc_t *A, magma_desc_t *T, int BS, magma_sequence_t *sequence, magma_request_t *request);
void magma_pzgetmi2(MAGMA_enum idep, MAGMA_enum odep, MAGMA_enum storev, int m, int n, int mb, int nb, PLASMA_Complex64_t *A, magma_sequence_t *sequence, magma_request_t *request);
void magma_pzgetrf(magma_desc_t *A, magma_desc_t *L, int *IPIV, magma_sequence_t *sequence, magma_request_t *request);
#if defined(PRECISION_z) || defined(PRECISION_c)
void magma_pzhemm(MAGMA_enum side, MAGMA_enum uplo, PLASMA_Complex64_t alpha, magma_desc_t *A, magma_desc_t *B, PLASMA_Complex64_t beta, magma_desc_t *C, magma_sequence_t *sequence, magma_request_t *request);
void magma_pzherk(MAGMA_enum uplo, MAGMA_enum trans, double alpha, magma_desc_t *A, double beta, magma_desc_t *C, magma_sequence_t *sequence, magma_request_t *request);
void magma_pzher2k(MAGMA_enum uplo, MAGMA_enum trans, PLASMA_Complex64_t alpha, magma_desc_t *A, magma_desc_t *B, double beta, magma_desc_t *C, magma_sequence_t *sequence, magma_request_t *request);
#endif
void magma_pzlacpy(MAGMA_enum uplo, magma_desc_t *A, magma_desc_t *B, magma_sequence_t *sequence, magma_request_t *request);
void magma_pzlag2c(magma_desc_t *A, magma_desc_t *SB, magma_sequence_t *sequence, magma_request_t *request);
void magma_pzlange(MAGMA_enum norm, magma_desc_t *A, double *work, double *result, magma_sequence_t *sequence, magma_request_t *request);
#if defined(PRECISION_z) || defined(PRECISION_c)
void magma_pzlanhe(MAGMA_enum norm, MAGMA_enum uplo, magma_desc_t *A, double *work, double *result, magma_sequence_t *sequence, magma_request_t *request);
#endif
void magma_pzlansy(MAGMA_enum norm, MAGMA_enum uplo, magma_desc_t *A, double *work, double *result, magma_sequence_t *sequence, magma_request_t *request);
void magma_pzlauum(MAGMA_enum uplo, magma_desc_t *A, magma_sequence_t *sequence, magma_request_t *request);
#if defined(PRECISION_z) || defined(PRECISION_c)
void magma_pzplghe(double bump, magma_desc_t *A, unsigned long long int seed, magma_sequence_t *sequence, magma_request_t *request);
#endif
void magma_pzplgsy(PLASMA_Complex64_t bump, magma_desc_t *A, unsigned long long int seed, magma_sequence_t *sequence, magma_request_t *request );
void magma_pzplrnt(magma_desc_t *A, unsigned long long int seed, magma_sequence_t *sequence, magma_request_t *request );
void magma_pzpotrf(MAGMA_enum uplo, magma_desc_t *A, magma_sequence_t *sequence, magma_request_t *request);
void magma_pzshift(int, int, int, PLASMA_Complex64_t *, int *, int, int, magma_sequence_t*, magma_request_t*);
void magma_pzsymm(MAGMA_enum side, MAGMA_enum uplo, PLASMA_Complex64_t alpha, magma_desc_t *A, magma_desc_t *B, PLASMA_Complex64_t beta, magma_desc_t *C, magma_sequence_t *sequence, magma_request_t *request);
void magma_pzsyrk(MAGMA_enum uplo, MAGMA_enum trans, PLASMA_Complex64_t alpha, magma_desc_t *A, PLASMA_Complex64_t beta,  magma_desc_t *C, magma_sequence_t *sequence, magma_request_t *request);
void magma_pzsyr2k(MAGMA_enum uplo, MAGMA_enum trans, PLASMA_Complex64_t alpha, magma_desc_t *A, magma_desc_t *B, PLASMA_Complex64_t beta, magma_desc_t *C, magma_sequence_t *sequence, magma_request_t *request);
void magma_pztrmm(MAGMA_enum side, MAGMA_enum uplo, MAGMA_enum transA, MAGMA_enum diag, PLASMA_Complex64_t alpha, magma_desc_t *A, magma_desc_t *B, magma_sequence_t *sequence, magma_request_t *request);
void magma_pztrsm(MAGMA_enum side, MAGMA_enum uplo, MAGMA_enum transA, MAGMA_enum diag, PLASMA_Complex64_t alpha, magma_desc_t *A, magma_desc_t *B, magma_sequence_t *sequence, magma_request_t *request);
void magma_pztrsmpl(magma_desc_t *A, magma_desc_t *B, magma_desc_t *L, int *IPIV, magma_sequence_t *sequence, magma_request_t *request);
void magma_pztrtri(MAGMA_enum uplo, MAGMA_enum diag, magma_desc_t *A, magma_sequence_t *sequence, magma_request_t *request);
void magma_pzungqr(magma_desc_t *A, magma_desc_t *Q, magma_desc_t *T, magma_sequence_t *sequence, magma_request_t *request);
void magma_pzungqrrh(magma_desc_t *A, magma_desc_t *Q, magma_desc_t *T, int BS, magma_sequence_t *sequence, magma_request_t *request);
void magma_pzunglq(magma_desc_t *A, magma_desc_t *Q, magma_desc_t *T, magma_sequence_t *sequence, magma_request_t *request);
void magma_pzunglqrh(magma_desc_t *A, magma_desc_t *Q, magma_desc_t *T, int BS, magma_sequence_t *sequence, magma_request_t *request);
void magma_pzunmqr(MAGMA_enum side, MAGMA_enum trans, magma_desc_t *A, magma_desc_t *B, magma_desc_t *T, magma_sequence_t *sequence, magma_request_t *request);
void magma_pzunmqrrh(MAGMA_enum side, MAGMA_enum trans, magma_desc_t *A, magma_desc_t *B, magma_desc_t *T, int BS, magma_sequence_t *sequence, magma_request_t *request);
void magma_pzunmlq(MAGMA_enum side, MAGMA_enum trans, magma_desc_t *A, magma_desc_t *B, magma_desc_t *T, magma_sequence_t *sequence, magma_request_t *request);
void magma_pzunmlqrh(MAGMA_enum side, MAGMA_enum trans, magma_desc_t *A, magma_desc_t *B, magma_desc_t *T, int BS, magma_sequence_t *sequence, magma_request_t *request);

#undef PRECISION_z
