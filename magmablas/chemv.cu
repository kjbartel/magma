/*
    -- MAGMA (version 1.5.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date September 2014
       
       csymv.cu is nearly identical to chemv.cu, just change names and drop cuConjf.
       
       @generated from zhemv.cu normal z -> c, Wed Sep 17 15:08:23 2014
       
       @author Mark Gates
*/
#include "common_magma.h"

#define PRECISION_c

#define NB_X         64
#define NB_Y          4
#define bank_shift   33
#define quarter_NB_X 16
#define half_NB_X    32


/*******************************************************************************
    Lower case, compute block multiply, work = A*x, for any size n:
    
    [ A11*x1   A12*x2             A13*x3                    ]   [ A11 A12 A13 ]   [ x1 ]
    [  ---    (A21*x1 + A22*x2)   A23*x3                    ] = [ A21 A22 A23 ] * [ x2 ]
    [  ---      ---              (A31*x1 + A32*x2 + A33*x3) ]   [ A31 A32 A33 ]   [ x3 ]
    
    Uses a 64x4 thread block.
    For     diagonal tiles, covers a 64x64 tile using three 32x32 tiles (plus one gets transposed).
    For off-diagonal tiles, covers a 64x64 tile using four  64x16 tiles.
    In both cases, each thread multiplies 4 elements.
    
    For rows past the bottom of the matrix, the A pointer is adjusted to be the
    last valid row of A, which multiple threads will read.
    Extra rows are ignored when saving results to work.
    Columns past the right edge are explicitly ignored when loading.
    x values past the bottom are set to zero, thus, extra columns are zeroed
    when multiplying.
    ********************************************************************/
__global__ void
chemv_kernel_L(
    int n,
    const magmaFloatComplex * __restrict__ A, int lda,
    const magmaFloatComplex * __restrict__ x, int incx,
    magmaFloatComplex * __restrict__ work)
{
#if defined(PRECISION_s) || defined(PRECISION_d) || defined(PRECISION_c) || (__CUDA_ARCH__ >= 200)

    // treats sA as 16x64 block
    #define sA16(i_, j_) (sA[(i_)][(j_)])  // i.e., sA[ (i_)*(NB_X+3) + (j_) ]
    
    // treats sA as 32x32 block
    #define sA32(i_, j_) (sA[0][(i_) + bank_shift*(j_)])
    
    // 64x4 thread block
    const int tx  = threadIdx.x;
    const int ty  = threadIdx.y;
    const int blk = blockIdx.x;
    const int blk_ind = NB_X * blk;
    const int td  = NB_X * ty + tx;

    // 32x8 thread block
    const int tx2 = td % half_NB_X;
    const int ty2 = td / half_NB_X;

    // If this blk has fewer than NB_X rows, partial is the number of valid rows,
    // so tx = 0, ..., partial-1 are valid rows, and tx >= partial are invalid.
    // Else, partial == 0.
    const int partial = (blk == gridDim.x - 1 ? (n % NB_X) : 0);
    
    magmaFloatComplex psum, psum2;
    magmaFloatComplex total = MAGMA_C_ZERO;

    // sA is used as a 32x32 block, sA32(i,j),
    // and as a 16x64 block, sA16(i,j), in different parts of the code.
    // sA must be at least half_NB_X*bank_shift = 32x33 = 1056;
    // quarter_NB_X*(NB_X + 2) = 16*(64 + 2) = 1056
    __shared__ magmaFloatComplex sA [quarter_NB_X][NB_X + 3]; /* Why +3? seems it only needs +2. Does +3 reduce bank conflicts? */
    __shared__ magmaFloatComplex sx [NB_X];  // for x[ blk ]
    __shared__ magmaFloatComplex sx2[NB_X];  // for x[ blk2 ], which cycles over all blocks left of diag

    magmaFloatComplex rA[4];
    magmaFloatComplex psums[4];

    // --------------------
    // load 64x1 block x(blk_ind + 0:63) into sx
    x += (blk_ind + tx)*incx;  // x is x(blk_ind + tx)
    if ( ty == 0 ) {
        if ( partial && tx >= partial ) {
            sx[tx] = MAGMA_C_ZERO;
        }
        else {
            sx[tx] = x[0];
        }
    }

    // --------------------
    // move to 32x32 diag block
    A += blk_ind * (lda + 1);  // A is A(blk_ind, blk_ind)
    A += ty2*lda + tx2;        // A is A(blk_ind + tx2, blk_ind + ty2)

    // load 32x32 diag block A(blk_ind + 0:31, blk_ind + 0:31) into sA,
    // as four 32x8 sections one after another:
    // columns 0:7, then 8:15, then 16:23, then 24:31
    if ( partial ) {
        if ( tx2 >= partial ) {
            A = A - tx2 + (partial - 1);
        }
        #pragma unroll
        for(int j=0; j < half_NB_X; j += 8) {
            if ( ty2+j < partial ) {
                sA32(tx2, ty2 + j) = A[j*lda];
            }
        }
        if ( tx2 >= partial ) {
            A = A + tx2 - (partial - 1);
        }
    }
    else {
        #pragma unroll
        for(int j=0; j < half_NB_X; j += 8) {
            sA32(tx2, ty2 + j) = A[j*lda];
        }
    }
    __syncthreads();

    // symmetrize 32x32 diag block, copying lower to upper triangle,
    // as four 32x8 sections in parallel:
    // columns 0,4,8,12,16,20,24,28; then 1,5,...,29; then 2,6,...,30, then 3,7,...,31
    #pragma unroll
    for(int j=ty2*4; j < ty2*4 + 4; j++) {
        if ( j < tx2 )
            sA32(j, tx2) = cuConjf( sA32(tx2, j) );
    }
    __syncthreads();

    // multiply 32x32 diag block * x
    // each thread does partial row sA(tx2, ty2*4 : ty2*4 + 3)
    psum = MAGMA_C_ZERO;
    #pragma unroll
    for(int j=0; j < 4; j++) {
        psum += sA32(tx2, ty2*4 + j) * sx[ty2*4 + j];
    }
    __syncthreads();

    // store partial row sums
    sA32(ty2, tx2) = psum;
    __syncthreads();

    // sum up partial row sums, so thread (tx2,0) has total for row (blk_ind + tx2)
    if ( ty2 == 0 ) {
        total = sA32(0, tx2) + sA32(1, tx2)
              + sA32(2, tx2) + sA32(3, tx2)
              + sA32(4, tx2) + sA32(5, tx2)
              + sA32(6, tx2) + sA32(7, tx2);
    }
    __syncthreads();

    // --------------------
    // move to next 32x32 diag block, then repeat steps from first diag block
    A += half_NB_X + half_NB_X*lda;  // A is A(blk_ind + NB/2 + tx2, blk_ind + NB/2 + ty2)

    // load 32x32 diag block A[block + 0:31, block + 0:31] into sA
    if ( partial ) {
        if ( tx2 + half_NB_X >= partial ) {
            A = A - (tx2 + half_NB_X) + (partial - 1);
        }
        #pragma unroll
        for(int j=0; j < half_NB_X; j += 8) {
            if ( ty2+j + half_NB_X < partial ) {
                sA32(tx2, ty2 + j) = A[j*lda];
            }
        }
        if ( tx2 + half_NB_X >= partial ) {
            A = A + (tx2 + half_NB_X) - (partial - 1);
        }
    }
    else {
        #pragma unroll
        for(int j=0; j < half_NB_X; j += 8) {
            sA32(tx2, ty2 + j) = A[j*lda];
        }
    }
    __syncthreads();

    // symmetrize 32x32 diag block, copying lower to upper triangle
    #pragma unroll
    for(int j=ty2*4; j < ty2*4 + 4; j++) {
        if ( j < tx2 )
            sA32(j, tx2) = cuConjf( sA32(tx2, j) );
    }
    __syncthreads();

    // multiply 32x32 diag block * x
    psum = MAGMA_C_ZERO;
    #pragma unroll
    for(int j=0; j < 4; j++) {
        psum += sA32(tx2, ty2*4 + j) * sx[half_NB_X + ty2*4 + j];
    }
    __syncthreads();
    
    // store partial row sums
    sA32(ty2, tx2) = psum;
    __syncthreads();

    // sum up partial row sums, so thread (tx2,1) has total for row (blk_ind + NB/2 + tx2)
    if ( ty2 == 1 ) {
        total = sA32(0, tx2) + sA32(1, tx2)
              + sA32(2, tx2) + sA32(3, tx2)
              + sA32(4, tx2) + sA32(5, tx2)
              + sA32(6, tx2) + sA32(7, tx2);
    }
    __syncthreads();

    // --------------------
    // move to off-diag 32x32 block
    A -= half_NB_X*lda;  // A is A(blk_ind + NB/2 + tx2, blk_ind + ty2)

    // load 32x32 block of A into sA,
    // as four 32x8 sections one after another:
    // columns 0:7, then 8:15, then 16:23, then 24:31
    if ( partial ) {
        if ( tx2 + half_NB_X >= partial ) {
            A = A - (tx2 + half_NB_X) + (partial - 1);
        }
        #pragma unroll
        for(int j=0; j < half_NB_X; j += 8) {
            if ( ty2+j < partial ) {
                sA32(tx2, ty2 + j) = A[j*lda];
            }
        }
        if ( tx2 + half_NB_X >= partial ) {
            A = A + (tx2 + half_NB_X) - (partial - 1);
        }
    }
    else {
        #pragma unroll
        for(int j=0; j < half_NB_X; j += 8) {
            sA32(tx2, ty2 + j) = A[j*lda];
        }
    }
    __syncthreads();

    // multiply 32x32 block (below diag)
    psum = MAGMA_C_ZERO;
    #pragma unroll
    for(int j=0; j < 4; j++) {
        psum += sA32(tx2, ty2 + j*8) * sx[j*8 + ty2];
    }
    //__syncthreads();  // no sync needed here

    // multiply transposed 32x32 block (above diag)
    psum2 = MAGMA_C_ZERO;
    #pragma unroll
    for(int j=0; j < 4; j++) {
        psum2 += cuConjf( sA32(ty2*4 + j, tx2) ) * sx[half_NB_X + ty2*4 + j];
    }
    __syncthreads();

    // store partial sums for non-transposed 32x32 block
    sA32(ty2, tx2) = psum;
    __syncthreads();
    
    // sum up partial row sums, so thread (tx2,1) has total for row (blk_ind + NB/2 + tx2)
    if ( ty2 == 1 ) {
        total = total
              + sA32(0, tx2) + sA32(1, tx2)
              + sA32(2, tx2) + sA32(3, tx2)
              + sA32(4, tx2) + sA32(5, tx2)
              + sA32(6, tx2) + sA32(7, tx2);
    }
    __syncthreads();

    // store partial sums for transposed 32x32 block
    sA32(ty2, tx2) = psum2;
    __syncthreads();
    
    // sum up partial row sums, so thread (tx2,0) has total for row (blk_ind + tx2)
    if ( ty2 == 0 ) {
        total = total
              + sA32(0, tx2) + sA32(1, tx2)
              + sA32(2, tx2) + sA32(3, tx2)
              + sA32(4, tx2) + sA32(5, tx2)
              + sA32(6, tx2) + sA32(7, tx2);
    }
    __syncthreads();
    
    // --------------------
    // move to left most 64x64 block in block row, and
    // switch thread offset from (tx2,ty2) 32x8 block to (tx,ty) 64x4 block
    A -= half_NB_X;       // A is A(blk_ind + tx2, blk_ind + ty2)
    A -= ty2*lda + tx2;   // A is A(blk_ind, blk_ind)
    A -= blk_ind*lda;     // A is A(blk_ind, 0)
    A += 4*ty*lda + tx;   // A is A(blk_ind + tx, 4*ty)
    
    if ( partial && tx >= partial ) {
        A = A - tx + (partial - 1);
    }
    
    x -= blk_ind * incx;  // x is x(tx)

    // 16x16 thread block
    const int tx4 = td % quarter_NB_X;
    const int ty4 = td / quarter_NB_X;

    work += blk*lda + tx4;  // work is work(tx4, blk)
    
    for(int blk2=0; blk2 < blk; ++blk2) {
        // load 64x1 block x(blk2_ind + 0:63) into sx2
        // since this block is left of diagonal, x cannot be partial rows
        if ( ty == 0 ) {
            sx2[tx] = x[blk2*NB_X*incx];
        }
        __syncthreads();

        for( int k=0; k < 4; k++ ) {
            // load 64x16 block of A into rA, 4 elements per thread,
            // as four 64x4 sections in parallel:
            // columns 0,4,8,12; then 1,5,9,13; then 2,6,10,14; then 3,7,11,15
            // since this block is left of diagonal, it cannot be partial columns
            #pragma unroll
            for(int j=0; j < 4; j++) {
                rA[j] = A[j*lda];
            }

            // 1) multiply 64x16 block A * x2
            //    each thread does partial row rA(tx + 16*k, ty*4 + 16*k : ty*4 + 3 + 16*k)
            // 2) multiply transposed 16x64 block A**H * x,
            //    storing each product Aji*xi to sA(j,i)
            #pragma unroll
            for(int j=0; j < 4; j++) {
                total += rA[j] * sx2[quarter_NB_X*k + ty*4 + j];
                sA16(ty*4 + j, tx) = cuConjf( rA[j] ) * sx[tx];
            }
            __syncthreads();

            // do partial row sums for transposed 16x64 result
            // use 16x16 thread grid (tx4, ty4) instead of 64x4 (tx, ty)
            // sum sixteen 16x4 sections in parallel:
            // columns 0,4,8,...,60; then 1,5,...,61; then 2,6,...,62; then 3,7,...,63
            psum2 = MAGMA_C_ZERO;
            #pragma unroll
            for(int j=0; j < 4; j++) {
                psum2 += sA16(tx4, ty4*4 + j);
            }
            __syncthreads();

            // store partial row sums (locally)
            psums[k] = psum2;

            // move to next 64x16 block
            A += lda * quarter_NB_X;  // A is A(blk_ind + tx#, blk2*NB_x + k*NB_X/4 + 4*ty), # or partial
        }

        // store partial row sums
        #pragma unroll
        for(int k=0; k < 4; k++) {
            sA16(tx4, ty4 + quarter_NB_X*k) = psums[k];
        }
        __syncthreads();
        
        // sum up partial row sums and store final total to workspace
        // thread (tx4,ty4) where ty4 < 4 sums row tx4 + ty4*16
        // since this is the transposed block above the diagonal, it cannot be partial rows
        if ( ty4 < 4 ) {
            int k = ty4*quarter_NB_X;
            psum2 = sA16(tx4,  0 + k) + sA16(tx4,  1 + k)
                  + sA16(tx4,  2 + k) + sA16(tx4,  3 + k)
                  + sA16(tx4,  4 + k) + sA16(tx4,  5 + k)
                  + sA16(tx4,  6 + k) + sA16(tx4,  7 + k)
                  + sA16(tx4,  8 + k) + sA16(tx4,  9 + k)
                  + sA16(tx4, 10 + k) + sA16(tx4, 11 + k)
                  + sA16(tx4, 12 + k) + sA16(tx4, 13 + k)
                  + sA16(tx4, 14 + k) + sA16(tx4, 15 + k);
            work[blk2*NB_X + k] = psum2;  // store at work( blk2*NB_X + tx4 + ty4*16, blk )
        }
        __syncthreads();
    }

    work -= tx4;  // work is work(blk_ind)
    work += tx;   // work is work(blk_ind + tx)

    // store row sums
    sA16(ty, tx) = total;
    __syncthreads();
    
    // sum up final total for row tx
    if ( ty == 0 && (partial == 0 || tx < partial) ) {
        total = sA16(0, tx) + sA16(1, tx) + sA16(2, tx) + sA16(3, tx);
        work[blk*NB_X] = total;  // store at work( blk*NB_X + tx, blk )
    }
#endif  /* PRECISION_[sdc] || (__CUDA_ARCH__ >= 200) */
}


/**************************************************************
    Lower case, sum up final results
    
    On input:
           [ A11*x1   A12*x2             A13*x3                    ]
    work = [  ---    (A21*x1 + A22*x2)   A23*x3                    ]
           [  ---      ---              (A31*x1 + A32*x2 + A33*x3) ]
    
    On output:
              [ A11*x1 + A12*x2 + A13*x3 ]
    y = alpha*[ A11*x1 + A22*x2 + A23*x3 ] + beta*y
              [ A21*x1 + A22*x2 + A33*x3 ]
    
    
    Previously:
           [ A11*x1    ---                                         ]
    work = [ A12*x2  (A21*x1 + A22*x2)    ---                      ]
           [ A13*x3   A23*x3            (A31*x1 + A32*x2 + A33*x3) ]
    which doesn't work as well because A13*x3 has 64 rows,
    while A31*x1 has only n % NB rows. This is why it used to need
    lwork = lda*(blocks + 1) instead of lda*blocks.
    ********************************************************************/
__global__ void
chemv_kernel_L_sum(
    int n, magmaFloatComplex alpha,
    int lda,
    magmaFloatComplex beta,
    magmaFloatComplex * __restrict__ y, int incy,
    magmaFloatComplex * __restrict__ work )
{
    int tx  = threadIdx.x;
    int blk = blockIdx.x;
    int blk_ind = blk * NB_X;
    int ind     = blk_ind + tx;
    
    if ( ind < n ) {
        work += ind + blk*lda;
        magmaFloatComplex Ax = MAGMA_C_ZERO;
        for(int i = blk_ind; i < n; i += NB_X) {
            Ax += work[0];
            work += lda;
        }
        y[ind * incy] = beta*y[ind * incy] + alpha*Ax;
    }
}


/**************************************************************
 *  Lower case, launch kernels
 */
extern "C"
void magmablas_chemv_L(
    magma_int_t n, magmaFloatComplex alpha,
    const magmaFloatComplex *A, magma_int_t lda,
    const magmaFloatComplex *x, magma_int_t incx,
    magmaFloatComplex beta,
    magmaFloatComplex *y, magma_int_t incy,
    magmaFloatComplex *dwork)
{
    magma_int_t blocks = (n - 1)/NB_X + 1;
    dim3 grid( blocks, 1, 1 );

    dim3 threads( NB_X, NB_Y, 1 );
    chemv_kernel_L<<< grid, threads, 0, magma_stream >>>
        (n, A, lda, x, incx, dwork);

    dim3 threads_sum( NB_X, 1, 1 );
    chemv_kernel_L_sum<<< grid, threads_sum, 0, magma_stream >>>
        (n, alpha, lda, beta, y, incy, dwork);
}


/**
    Purpose
    -------
    magmablas_chemv_work performs the matrix-vector operation:

        y := alpha*A*x + beta*y,

    where alpha and beta are scalars, x and y are n element vectors and
    A is an n by n Hermitian matrix.

    Arguments
    ----------
    @param[in]
    uplo    magma_uplo_t.
            On entry, UPLO specifies whether the upper or lower
            triangular part of the array A is to be referenced as
            follows:
      -     = MagmaUpper:  Only the upper triangular part of A is to be referenced.
      -     = MagmaLower:  Only the lower triangular part of A is to be referenced.

    @param[in]
    n       INTEGER.
            On entry, N specifies the order of the matrix A.
            N must be at least zero.

    @param[in]
    alpha   COMPLEX.
            On entry, ALPHA specifies the scalar alpha.

    @param[in]
    A       COMPLEX array of DIMENSION ( LDA, n ).
            Before entry with UPLO = MagmaUpper, the leading n by n
            upper triangular part of the array A must contain the upper
            triangular part of the Hermitian matrix and the strictly
            lower triangular part of A is not referenced.
            Before entry with UPLO = MagmaLower, the leading n by n
            lower triangular part of the array A must contain the lower
            triangular part of the Hermitian matrix and the strictly
            upper triangular part of A is not referenced.
            Note that the imaginary parts of the diagonal elements need
            not be set and are assumed to be zero.

    @param[in]
    lda     INTEGER.
            On entry, LDA specifies the first dimension of A as declared
            in the calling (sub) program. LDA must be at least
            max( 1, n ).
            It is recommended that lda is multiple of 16. Otherwise
            performance would be deteriorated as the memory accesses
            would not be fully coalescent.

    @param[in]
    x       COMPLEX array of dimension at least
            ( 1 + ( n - 1 )*abs( INCX ) ).
            Before entry, the incremented array X must contain the n
            element vector x.

    @param[in]
    incx    INTEGER.
            On entry, INCX specifies the increment for the elements of
            X. INCX must not be zero.

    @param[in]
    beta    COMPLEX.
            On entry, BETA specifies the scalar beta. When BETA is
            supplied as zero then Y need not be set on input.

    @param[in, out]
    y       COMPLEX array of dimension at least
            ( 1 + ( n - 1 )*abs( INCY ) ).
            Before entry, the incremented array Y must contain the n
            element vector y. On exit, Y is overwritten by the updated
            vector y.

    @param[in]
    incy    INTEGER.
            On entry, INCY specifies the increment for the elements of
            Y. INCY must not be zero.

    @param[in]
    dwork   (workspace) COMPLEX array on the GPU, dimension (MAX(1, LWORK)),

    @param[in]
    lwork   INTEGER.
            The dimension of the array DWORK. LWORK >= LDA * ceil( N / NB_X ),
            where NB_X = 64.

    MAGMA implements chemv through two steps:
    1)  perform the multiplication in each thread block and put the
        intermediate value in dwork.
    2)  sum the intermediate values and store the final result in y.
    
    magamblas_chemv_work requires users to provide a workspace, while
    magmablas_chemv is a wrapper routine allocating the workspace inside the
    routine and provides the same interface as cublas.
    
    If users need to call chemv frequently, we suggest using
    magmablas_chemv_work instead of magmablas_chemv. As the overhead to
    allocate and free in device memory in magmablas_chemv would hurt performance.
    Our tests show that this penalty is about 10 Gflop/s when the matrix
    size is around 10000.

    @ingroup magma_cblas2
    ********************************************************************/
extern "C"
magma_int_t
magmablas_chemv_work(
    magma_uplo_t uplo, magma_int_t n,
    magmaFloatComplex alpha,
    const magmaFloatComplex *A, magma_int_t lda,
    const magmaFloatComplex *x, magma_int_t incx,
    magmaFloatComplex beta,
    magmaFloatComplex *y, magma_int_t incy,
    magmaFloatComplex *dwork, magma_int_t lwork)
{
#if defined(PRECISION_z)
    // z precision requires CUDA ARCH 2.x; call CUBLAS version instead.
    magma_int_t arch = magma_getdevice_arch();
    if ( arch < 200 ) {
        magma_chemv( uplo, n, alpha, A, lda, x, incx, beta, y, incy );
        return MAGMA_SUCCESS;
    }
#endif

    // --------------------
    // [sdc] precisions, or z precision with CUDA ARCH 2.x
    int upper = (uplo == MagmaUpper);

    magma_int_t blocks = (n - 1)/NB_X + 1;
    magma_int_t lwmin  = lda*blocks;

    /*
     * Test the input parameters.
     */
    magma_int_t info = 0;
    if ((! upper) && (uplo != MagmaLower)) {
        info = -1;
    } else if ( n < 0 ) {
        info = -2;
    } else if ( lda < max(1, n) ) {
        info = -5;
    } else if ( incx == 0 ) {
        info = -7;
    } else if ( incy == 0 ) {
        info = -10;
    } else if ( lwork < lwmin ) {
        info = -12;
    }
    
    if (info != 0) {
        magma_xerbla( __func__, -(info) );
        return info;
    }

    /*
     * Quick return if possible.
     */
    if ( (n == 0) || ( MAGMA_C_EQUAL(alpha, MAGMA_C_ZERO) && MAGMA_C_EQUAL(beta, MAGMA_C_ONE) ) )
        return info;

    /* TODO: Upper case is not implemented in MAGMA */
    if ( upper ) {
        magma_chemv( uplo, n, alpha, A, lda, x, incx, beta, y, incy);
    }
    else {
        magmablas_chemv_L(n, alpha, A, lda, x, incx, beta, y, incy, dwork);
    }
    return info;
}


/**
    Purpose
    -------
    magmablas_chemv performs the matrix-vector operation:

        y := alpha*A*x + beta*y,

    where alpha and beta are scalars, x and y are n element vectors and
    A is an n by n Hermitian matrix.

    Arguments
    ----------
    @param[in]
    uplo    magma_uplo_t.
            On entry, UPLO specifies whether the upper or lower
            triangular part of the array A is to be referenced as
            follows:
      -     = MagmaUpper:  Only the upper triangular part of A is to be referenced.
      -     = MagmaLower:  Only the lower triangular part of A is to be referenced.

    @param[in]
    n       INTEGER.
            On entry, N specifies the order of the matrix A.
            N must be at least zero.

    @param[in]
    alpha   COMPLEX.
            On entry, ALPHA specifies the scalar alpha.

    @param[in]
    A       COMPLEX array of DIMENSION ( LDA, n ).
            Before entry with UPLO = MagmaUpper, the leading n by n
            upper triangular part of the array A must contain the upper
            triangular part of the Hermitian matrix and the strictly
            lower triangular part of A is not referenced.
            Before entry with UPLO = MagmaLower, the leading n by n
            lower triangular part of the array A must contain the lower
            triangular part of the Hermitian matrix and the strictly
            upper triangular part of A is not referenced.
            Note that the imaginary parts of the diagonal elements need
            not be set and are assumed to be zero.

    @param[in]
    lda     INTEGER.
            On entry, LDA specifies the first dimension of A as declared
            in the calling (sub) program. LDA must be at least
            max( 1, n ).
            It is recommended that lda is multiple of 16. Otherwise
            performance would be deteriorated as the memory accesses
            would not be fully coalescent.

    @param[in]
    x       COMPLEX array of dimension at least
            ( 1 + ( n - 1 )*abs( INCX ) ).
            Before entry, the incremented array X must contain the n
            element vector x.

    @param[in]
    incx    INTEGER.
            On entry, INCX specifies the increment for the elements of
            X. INCX must not be zero.

    @param[in]
    beta    COMPLEX.
            On entry, BETA specifies the scalar beta. When BETA is
            supplied as zero then Y need not be set on input.

    @param[in, out]
    y       COMPLEX array of dimension at least
            ( 1 + ( n - 1 )*abs( INCY ) ).
            Before entry, the incremented array Y must contain the n
            element vector y. On exit, Y is overwritten by the updated
            vector y.

    @param[in]
    incy    INTEGER.
            On entry, INCY specifies the increment for the elements of
            Y. INCY must not be zero.

    @ingroup magma_cblas2
    ********************************************************************/
extern "C"
magma_int_t
magmablas_chemv(
    magma_uplo_t uplo, magma_int_t n,
    magmaFloatComplex alpha,
    const magmaFloatComplex *A, magma_int_t lda,
    const magmaFloatComplex *x, magma_int_t incx,
    magmaFloatComplex beta,
    magmaFloatComplex *y, magma_int_t incy)
{
#if defined(PRECISION_z)
    // z precision requires CUDA ARCH 2.x; call CUBLAS version instead.
    magma_int_t arch = magma_getdevice_arch();
    if ( arch < 200 ) {
        magma_chemv( uplo, n, alpha, A, lda, x, incx, beta, y, incy );
        return MAGMA_SUCCESS;
    }
#endif

    // --------------------
    // [sdc] precisions, or z precision with CUDA ARCH 2.x
    int upper = (uplo == MagmaUpper);

    /*
     * Test the input parameters.
     */
    magma_int_t info = 0;
    if ((! upper) && (uplo != MagmaLower)) {
        info = -1;
    } else if ( n < 0 ) {
        info = -2;
    } else if ( lda < max(1, n) ) {
        info = -5;
    } else if ( incx == 0 ) {
        info = -7;
    } else if ( incy == 0 ) {
        info = -10;
    }
    
    if (info != 0) {
        magma_xerbla( __func__, -(info) );
        return info;
    }

    /*
     * Quick return if possible.
     */
    if ( (n == 0) || ( MAGMA_C_EQUAL(alpha, MAGMA_C_ZERO) && MAGMA_C_EQUAL(beta, MAGMA_C_ONE) ) )
        return info;

    /* TODO: Upper case is not implemented in MAGMA */
    if ( upper ) {
        magma_chemv( uplo, n, alpha, A, lda, x, incx, beta, y, incy);
    }
    else {
        magmaFloatComplex *dwork;
        magma_int_t blocks = (n - 1)/NB_X + 1;
        magma_int_t lwork  = lda*blocks;

        magma_cmalloc( &dwork, lwork );
        if ( dwork == NULL ) {
            info = MAGMA_ERR_DEVICE_ALLOC;
            magma_xerbla( __func__, -(info) );
        }
        else {
            magmablas_chemv_L(n, alpha, A, lda, x, incx, beta, y, incy, dwork);
        }
        magma_free( dwork );
    }
    return info;
}
