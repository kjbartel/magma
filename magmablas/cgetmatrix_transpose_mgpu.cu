/*
    -- MAGMA (version 1.6.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date November 2014

       @generated from zgetmatrix_transpose_mgpu.cu normal z -> c, Sat Nov 15 19:53:59 2014
       @author Ichitaro Yamazaki
*/
#include "common_magma.h"

#define PRECISION_c

//
//    m, n - dimensions in the output (hA) matrix.
//             This routine copies the dAT matrix from the GPU
//             to hA on the CPU. In addition, the output matrix
//             is transposed. The routine uses a buffer of size
//             2*lddw*nb pointed to by dwork (lddw > m) on the GPU. 
//             Note that lda >= m and lddat >= n.
//
extern "C" void 
magmablas_cgetmatrix_transpose_mgpu(
    magma_int_t ngpu, magma_queue_t queues[][2],
    magmaFloatComplex_const_ptr const dAT[],   magma_int_t ldda,
    magmaFloatComplex                *hA,      magma_int_t lda,
    magmaFloatComplex_ptr             dwork[], magma_int_t lddw,
    magma_int_t m, magma_int_t n, magma_int_t nb)
{
#define hA(j)       (hA         + (j)*lda)
#define dwork(d, j) (dwork[(d)] + (j)*nb*lddw)
#define dAT(d, j)   (dAT[(d)]   + (j)*nb)

    magma_int_t nstreams = 2, d, j, j_local, id, ib;

    /* Quick return */
    if ( (m == 0) || (n == 0) )
        return;

    if (lda < m || ngpu*ldda < n || lddw < m){
        printf( "Wrong arguments in magmablas_cgetmatrix_transpose_mgpu (%d<%d), (%d*%d<%d), or (%d<%d).\n",
                (int) lda, (int) m, (int) ngpu, (int) ldda, (int) n, (int) lddw, (int) m );
        return;
    }
    
    /* Move data from GPU to CPU using two buffers; first transpose the data on the GPU */
    for(j=0; j<n; j+=nb){
       d       = (j/nb)%ngpu;
       j_local = (j/nb)/ngpu;
       id      = j_local%nstreams;
       magma_setdevice(d);

       ib = min(n-j, nb);
       magmablas_ctranspose_q( ib, m, dAT(d,j_local), ldda, dwork(d,id), lddw, queues[d][id] );
       magma_cgetmatrix_async( m, ib,
                               dwork(d, id), lddw,
                               hA(j),        lda, 
                               queues[d][id] );
    }
}
