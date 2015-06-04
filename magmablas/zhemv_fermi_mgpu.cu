/*
    -- MAGMA (version 1.2.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       June 2012

       @precisions normal z -> s d c

*/


#include "common_magma.h"
#define PRECISION_z

/*The version for tesla can be found in zhemv_tesla.cu */
#if (GPUSHMEM >= 200)

#define magmablas_zhemv_200_mgpu_64 magmablas_zhemv_mgpu_64

#define zhemv_bs         64
#define thread_x         64
#define thread_y          4
#define bank_shift       33
#define quarter_thread_x 16
#define half_thread_x    32

/*******************************************************************************
 *     Functions for each specific cases - Lower case
 */

__global__ void
magmablas_zhemv_200_L_special_mgpu( magma_int_t n, cuDoubleComplex alpha,
                               cuDoubleComplex *A, magma_int_t lda,
                               cuDoubleComplex *x, magma_int_t incx,
                               cuDoubleComplex  beta,
                               cuDoubleComplex *y, magma_int_t incy,
                               cuDoubleComplex *WC, 
                         magma_int_t my_gpu_id,
                         magma_int_t num_gpus, 
                         magma_int_t nb)
{
    magma_int_t tx   = threadIdx.x ;
    magma_int_t ty   = threadIdx.y ;
    magma_int_t blkc = blockIdx.x ;

    if(blkc < my_gpu_id)
    {
    return;
    }

    cuDoubleComplex res  = MAGMA_Z_ZERO;
    cuDoubleComplex res_ = MAGMA_Z_ZERO;
    cuDoubleComplex res1 = MAGMA_Z_ZERO;
    cuDoubleComplex res2 = MAGMA_Z_ZERO;

    __shared__ cuDoubleComplex la   [quarter_thread_x][thread_x+2];
    __shared__ cuDoubleComplex buff [thread_x];
    __shared__ cuDoubleComplex buff2 [thread_x];

    cuDoubleComplex tr[4];    
    cuDoubleComplex b[4];

    magma_int_t break_d   =  thread_x * blkc;
    const magma_int_t td  = (thread_x * ty ) + tx;
    magma_int_t       tx_ = td % half_thread_x;
    magma_int_t       ty_ = td / half_thread_x;

    WC +=  break_d + tx;
    x  += (break_d + tx ) * incx;
    A  +=  break_d ;
    A  += ty_* lda + tx_ ;

    if( ty == 0 ){
        buff[tx] = x[0];
    } // obtain the vector x store in buff;

    tx = tx_ ; ty = ty_ ;

    magma_int_t flag = 0;
    
    if ( (blkc % num_gpus) == my_gpu_id) 
    {
    A += lda * (blkc/num_gpus) * thread_x; // change

    #pragma unroll
    for(magma_int_t j =0; j<half_thread_x; j +=8)
        la[0][ bank_shift * (ty_+j) + tx_] =  A[ j * lda];
    __syncthreads();

    #pragma unroll
    for(magma_int_t  i=ty_*4; i<(ty_ * 4 + 4)  ; i++){
        if ( i < tx_ )   
        {
            la[0][bank_shift * tx_ + i] = cuConj( la[0][ i * bank_shift + tx_] ) ;
        }
        /*
        else
            la[0][bank_shift * tx_ + i] = la[0][ bank_shift * tx_ + i]  ;
        */
    }
    __syncthreads();

    #pragma unroll
    for(magma_int_t j=0; j < 4 ; j++)
        res+= cuConj( la[0][bank_shift * tx_ + j + ty_ * 4] ) * buff[j + ty_ * 4];
    __syncthreads();

    la[0][bank_shift*tx_+ty_]= res ;
    __syncthreads();

    if( ty_== 0 )
      res1 = la[0][tx_*bank_shift+0]+la[0][tx_*bank_shift+1]
        +    la[0][tx_*bank_shift+2]+la[0][tx_*bank_shift+3]
        +    la[0][tx_*bank_shift+4]+la[0][tx_*bank_shift+5]
        +    la[0][tx_*bank_shift+6]+la[0][tx_*bank_shift+7];
    /*
    else
        {
            MAGMA_Z_SET2REAL(res1,0);
        }

     */
    __syncthreads();


    MAGMA_Z_SET2REAL(res, 0) ; 

    A+= half_thread_x + half_thread_x *lda ;// the right bottom 1/4 corner

    #pragma unroll
    for(magma_int_t j =0; j<half_thread_x; j+=8)
        la[0][bank_shift*(ty_+j)+tx_] = A[ j * lda];
    __syncthreads();

    #pragma unroll
    for(magma_int_t  i=ty_*4; i<(4+ty_*4) ; i++){
        if ( i < tx_ )   {
            la[0][bank_shift*tx_+i] = cuConj( la[0][bank_shift*i+tx_] ) ;
        }
        /*
        else
            la[0][bank_shift*tx_+i] = la[0][bank_shift*tx_+i]  ;
        */
    }
    __syncthreads();

    #pragma unroll
    for(magma_int_t j=0; j < 4 ; j++)
        res+= cuConj( la[0][bank_shift*tx_+j+ty_*4] ) * buff[half_thread_x + j + 4 * ty_];
    __syncthreads();
    la[0][bank_shift*tx_+ty_]= res ;
    __syncthreads();

    if( ty_== 1 )
        res2 = la[0][tx_*bank_shift+0]+la[0][tx_*bank_shift+1]
          +    la[0][tx_*bank_shift+2]+la[0][tx_*bank_shift+3]
          +    la[0][tx_*bank_shift+4]+la[0][tx_*bank_shift+5]
          +    la[0][tx_*bank_shift+6]+la[0][tx_*bank_shift+7];
    /*
    else
    {
        MAGMA_Z_SET2REAL(res2,0);
    }
    */
    __syncthreads();

    MAGMA_Z_SET2REAL(res,0);

    A-=half_thread_x *lda ;// the left bottom 1/4 corner

    MAGMA_Z_SET2REAL(res_,0);

    #pragma unroll
    for(magma_int_t j=0; j<half_thread_x; j+=8)
        tr[j/8] = A[ j * lda];

    #pragma unroll
    for(magma_int_t j=0; j < 4 ; j++){
        res += tr[j] * buff[ j*8 + ty_];
        la[0][bank_shift*(ty_+j*8)+tx_] = tr[j];
    }
    __syncthreads();

    #pragma unroll
    for(magma_int_t j=0; j < 4 ; j++)
        res_+= cuConj(la[0][bank_shift*tx_+j+ty_*4]) * buff[half_thread_x +j+ty_*4];
    __syncthreads();

    la[0][bank_shift*tx_+ty_]= res ;
    __syncthreads();
    if( ty_ == 1 )
        res2 = res2 
            +  la[0][tx_*bank_shift+0]+la[0][tx_*bank_shift+1]
            +  la[0][tx_*bank_shift+2]+la[0][tx_*bank_shift+3]
            +  la[0][tx_*bank_shift+4]+la[0][tx_*bank_shift+5]
            +  la[0][tx_*bank_shift+6]+la[0][tx_*bank_shift+7];
    else
        {
            MAGMA_Z_SET2REAL(res2,0);
        }
    
    __syncthreads();

    la[0][bank_shift*tx_+ty_]= res_ ;
    __syncthreads();
    if( ty_ == 0 ) {
        res1 = res1
            +  la[0][tx_*bank_shift+0]+la[0][tx_*bank_shift+1]
            +  la[0][tx_*bank_shift+2]+la[0][tx_*bank_shift+3]
            +  la[0][tx_*bank_shift+4]+la[0][tx_*bank_shift+5]
            +  la[0][tx_*bank_shift+6]+la[0][tx_*bank_shift+7];
    }
    else
        {
            MAGMA_Z_SET2REAL(res1,0);
        }
    
    A-=half_thread_x; // go back to the original pointer

      A -= lda * (blkc/num_gpus) * thread_x; 
    
      flag = 1;
    }

    __syncthreads();
    tx = threadIdx.x ;
    ty = threadIdx.y ;

    if( ty_ == 0  && ty == 0  )
        res = res1 ;
    else if( ty_ == 1  && ty == 0  )
        res = res2 ;
    else
        {
            MAGMA_Z_SET2REAL(res,0);
        }

    A-=ty_* lda  ;
    A-=tx_;

    x= x - blkc * thread_x  *incx  ;

    //x= x- tx*incx;

    A+=4 * ty* lda  ;
    A+=tx;

    magma_int_t wc_c = my_gpu_id ;
    magma_int_t count = 0 ;

    tx_ = td % quarter_thread_x ;
    ty_ = td / quarter_thread_x ;

    WC-=tx ;
    WC+=tx_;


    magma_int_t num_blocks_iters = (blkc +1) /num_gpus - flag;
    
    if( my_gpu_id < ( (blkc+1) % num_gpus) )
    {
    num_blocks_iters += 1;
    }

    x += (my_gpu_id ) * nb ;


    if( blkc > my_gpu_id)

    for(magma_int_t s=0; s<num_blocks_iters; s++)
        {
            MAGMA_Z_SET2REAL(res_,0);
            count++;
            if(ty == 0)
            {
                buff2[tx] = x[0];
            }    
             __syncthreads();

            #pragma unroll
            for( magma_int_t k=0;k<4;k++)
            {


                #pragma unroll
                for(magma_int_t j=0; j < 4 ; j++)
                    tr[j] = A[j*lda];

                #pragma unroll
                for(magma_int_t j=0; j < 4 ; j++)
                {
                    res += tr[j] * buff2[ quarter_thread_x * k + ty * 4 + j];
                    la[( j + ty * 4)][tx] = cuConj(tr[j]) * buff[tx];
                }
                __syncthreads();


                MAGMA_Z_SET2REAL(res_,0);

                #pragma unroll
                for(magma_int_t j=0; j < 4 ; j++)
                {
                    res_+=la[tx_][ty_*4+j] ;
                }
                b[k] = res_ ;
                __syncthreads();

                A += lda * quarter_thread_x ;
            }

            #pragma unroll
            for(magma_int_t k=0; k < 4 ; k++){
                la[tx_][ty_+quarter_thread_x*k]= b[k] ;
            }
            __syncthreads();
            if( ty_ < 4 ) {
                magma_int_t k = ty_*quarter_thread_x;
                res_ = la[tx_][0+k] + la[tx_][1+k]
                    +  la[tx_][2+k] + la[tx_][3+k]
                    +  la[tx_][4+k] + la[tx_][5+k]
                    +  la[tx_][6+k] + la[tx_][7+k]
                    +  la[tx_][8+k] + la[tx_][9+k]
                    +  la[tx_][10+k]+ la[tx_][11+k]
                    +  la[tx_][12+k]+ la[tx_][13+k]
                    +  la[tx_][14+k]+ la[tx_][15+k];
                WC[k + wc_c*lda ] =   res_;
            }

            wc_c += num_gpus;
        x += num_gpus * nb;
            __syncthreads();

        }


    WC+=tx ;
    WC-=tx_;

    la[ty][tx]= res ;
    __syncthreads();
    if( ty == 0 ) {
        res = la[0][tx]+ la[1][tx]
            + la[2][tx]+ la[3][tx];
        WC[0+lda*(blkc)  ] =  res;
    }
}

/**************************************************************
 *    Lower case for generic sizes
 */
__global__ void
magmablas_zhemv_200_L_generic_mgpu(magma_int_t n, cuDoubleComplex alpha,
                              cuDoubleComplex *A, magma_int_t lda,
                              cuDoubleComplex *x, magma_int_t incx,
                              cuDoubleComplex beta,
                              cuDoubleComplex *y, magma_int_t incy,
                              cuDoubleComplex *WC,
                              magma_int_t m_mod_thread_x,
                         magma_int_t my_gpu_id,
                         magma_int_t num_gpus, 
                         magma_int_t nb)
{
    magma_int_t tx   = threadIdx.x ;
    magma_int_t ty   = threadIdx.y ;
    magma_int_t blkc = blockIdx.x ;

    if(blkc < my_gpu_id)
    {
    return;
    }

    cuDoubleComplex res  = MAGMA_Z_ZERO;
    cuDoubleComplex res_ = MAGMA_Z_ZERO;
    cuDoubleComplex res1 = MAGMA_Z_ZERO;
    cuDoubleComplex res2 = MAGMA_Z_ZERO;

    __shared__ cuDoubleComplex la   [quarter_thread_x][thread_x+2];
    __shared__ cuDoubleComplex buff [thread_x];
    __shared__ cuDoubleComplex buff2 [thread_x];


    cuDoubleComplex tr[4];
    cuDoubleComplex b[4];

    magma_int_t break_d   =  thread_x * blkc;
    const magma_int_t td  = (thread_x * ty ) + tx;
    magma_int_t       tx_ = td % half_thread_x;
    magma_int_t       ty_ = td / half_thread_x;

    WC+=  break_d + tx;
    x += (break_d + tx ) * incx;
    A +=  break_d ;
    A += lda * ty_;

    magma_int_t trackA ;
    if( blkc == ( gridDim.x - 1 ) ) {
        if( ty == 0 ){
            if( tx > m_mod_thread_x )
            {
                MAGMA_Z_SET2REAL(buff[tx],0);
            }
            else
                buff[tx]  = x[0];
        }
        if ( tx_ > m_mod_thread_x )
            trackA=m_mod_thread_x;
        else
            trackA=tx_;
        A += trackA ;
    }
    else {
        if( ty == 0 ){
            buff[tx]  = x[0];
        }
        trackA = tx_;
        A += trackA ;
    }


    magma_int_t flag = 0;
    
    if ( (blkc % num_gpus) == my_gpu_id) 
    {
    A += lda * (blkc/num_gpus) * thread_x; // change
    // Somehow merging these two if - else creates problem
    // It could be a potential bug -- from synchronization or from cuda or compiler
    if( blkc == ( gridDim.x - 1 ) ) {
        #pragma unroll
        for(magma_int_t j =0; j<half_thread_x; j+=8){
            if( ( ty_ + j ) > m_mod_thread_x )
            {
                MAGMA_Z_SET2REAL(la[0][bank_shift*(ty_+j)+tx_], 9999);
            }
            else
                la[0][bank_shift*(ty_+j)+tx_] =  A[ j * lda];
        }
        A-=trackA;
    }
    else {
        #pragma unroll
        for(magma_int_t j =0; j<half_thread_x; j+=8){
            la[0][bank_shift*(ty_+j)+tx_] = A[ j * lda];
        }
    }
    tx = tx_ ;
    ty = ty_ ;
    __syncthreads();

    #pragma unroll
    for(magma_int_t  i=ty_*4; i<(ty_*4+4)  ; i++){
        if ( i < tx_ )   {
            la[0][bank_shift*tx_+i] = cuConj(la[0][i*bank_shift+tx_]) ;
        }
        else
            la[0][bank_shift*tx_+i] = la[0][bank_shift*tx_+i]  ;
    }
    __syncthreads();

    #pragma unroll
    for(magma_int_t j=0; j < 4 ; j++)
        res += cuConj(la[0][bank_shift*tx_+j+ty_*4])* buff[j+ty_*4];
    __syncthreads();

    la[0][bank_shift*tx_+ty_]= res ;
    __syncthreads();
    if( ty_== 0 )
        res1 = la[0][tx_*bank_shift+0] 
            +  la[0][tx_*bank_shift+1]
            +  la[0][tx_*bank_shift+2]
            +  la[0][tx_*bank_shift+3]
            +  la[0][tx_*bank_shift+4]
            +  la[0][tx_*bank_shift+5]
            +  la[0][tx_*bank_shift+6]
            +  la[0][tx_*bank_shift+7];
    else
    {
        MAGMA_Z_SET2REAL(res1,0);
    }
    __syncthreads();


    MAGMA_Z_SET2REAL(res,0);

    if( blkc == ( gridDim.x - 1 ) ) {
        if ( (tx_+half_thread_x) > m_mod_thread_x )
            trackA = m_mod_thread_x;
        else
            trackA = tx_ + half_thread_x;
        A+= trackA+half_thread_x*lda ;

        #pragma unroll
        for(magma_int_t j =0; j<half_thread_x; j+=8){
            if( ( ty_ + j+half_thread_x ) > m_mod_thread_x )
            {
                MAGMA_Z_SET2REAL(la[0][bank_shift*(ty_+j)+tx_], 99999);
            }
            else
                la[0][bank_shift*(ty_+j)+tx_] =  A[ j * lda];
        }

        A-= trackA+half_thread_x*lda ;
        A+=tx_ ;
        A+= half_thread_x + half_thread_x *lda ;
    }
    else {
        A+= half_thread_x + half_thread_x *lda ;

        #pragma unroll
        for(magma_int_t j =0; j<half_thread_x; j+=8){
            la[0][bank_shift*(ty_+j)+tx_] = A[ j * lda];
        }
    }

    __syncthreads();
    #pragma unroll
    for(magma_int_t  i=ty_*4; i<(4+ty_*4) ; i++){
        if ( i < tx_ )   {
            la[0][bank_shift*tx_+i] = cuConj(la[0][bank_shift*i+tx_]) ;
        }
        else
            la[0][bank_shift*tx_+i] = la[0][bank_shift*tx_+i]  ;
    }
    __syncthreads();

    #pragma unroll
    for(magma_int_t j=0; j < 4 ; j++)
        res+= cuConj(la[0][bank_shift*tx_+j+ty_*4]) * buff[half_thread_x + j + 4 * ty_];
    __syncthreads();

    la[0][bank_shift*tx_+ty_]= res ;
    __syncthreads();

    if( ty_== 1 )
        res2 = la[0][tx_*bank_shift+0]
            +  la[0][tx_*bank_shift+1]
            +  la[0][tx_*bank_shift+2]
            +  la[0][tx_*bank_shift+3]
            +  la[0][tx_*bank_shift+4]
            +  la[0][tx_*bank_shift+5]
            +  la[0][tx_*bank_shift+6]
            +  la[0][tx_*bank_shift+7];
    else
    {
        MAGMA_Z_SET2REAL(res2,0);
    }
    __syncthreads();

    MAGMA_Z_SET2REAL(res,0);
    MAGMA_Z_SET2REAL(res_,0);

    A-=half_thread_x *lda ;
    if( blkc == ( gridDim.x - 1 ) ) {
        A-=tx_;
        if ( tx_ > m_mod_thread_x )
            trackA=m_mod_thread_x;
        else
            trackA=tx_;
        A+= trackA ;

        #pragma unroll
        for(magma_int_t j =0; j<half_thread_x; j+=8)
            if( ( ty_ + j ) > m_mod_thread_x )
            {
                MAGMA_Z_SET2REAL(tr[j/8], 99999);
            }
            else
                tr[j/8] = A[ j * lda];
        A-=trackA;
        A+=tx_;
    }
    else {
        #pragma unroll
        for(magma_int_t j =0; j<half_thread_x; j+=8)
            tr[j/8] = A[ j * lda];
    }
    __syncthreads();

    #pragma unroll
    for(magma_int_t j=0; j < 4 ; j++){
        res+= tr[j] * buff[ j*8 + ty_];
        la[0][bank_shift*(ty_+j*8)+tx_] = tr[j];
    }
    __syncthreads();

    #pragma unroll
    for(magma_int_t j=0; j < 4 ; j++)
        res_+= cuConj(la[0][bank_shift*tx_+j+ty_*4]) * buff[half_thread_x +j+ty_*4];
    __syncthreads();


    la[0][bank_shift*tx_+ty_]= res ;
    __syncthreads();
    if( ty_ == 1 )
        res2 = res2
            +  la[0][tx_*bank_shift+0]
            +  la[0][tx_*bank_shift+1]
            +  la[0][tx_*bank_shift+2]
            +  la[0][tx_*bank_shift+3]
            +  la[0][tx_*bank_shift+4]
            +  la[0][tx_*bank_shift+5]
            +  la[0][tx_*bank_shift+6]
            +  la[0][tx_*bank_shift+7];
    else
    {
        MAGMA_Z_SET2REAL(res2,0);
    }
    __syncthreads();

    la[0][bank_shift*tx_+ty_]= res_ ;
    __syncthreads();

    if( ty_ == 0 ) {
        res1 = res1
            +  la[0][tx_*bank_shift+0]
            +  la[0][tx_*bank_shift+1]
            +  la[0][tx_*bank_shift+2]
            +  la[0][tx_*bank_shift+3]
            +  la[0][tx_*bank_shift+4]
            +  la[0][tx_*bank_shift+5]
            +  la[0][tx_*bank_shift+6]
            +  la[0][tx_*bank_shift+7];
    }
    else
    {
        MAGMA_Z_SET2REAL(res1,0);
    }
    A-=half_thread_x;

      A -= lda * (blkc/num_gpus) * thread_x; 
    
      flag = 1;
    }

    __syncthreads();
    tx = threadIdx.x ;
    ty = threadIdx.y ;

    if( ty_ == 0  && ty == 0  )
        res = res1 ;
    else if( ty_ == 1  && ty == 0  )
        res = res2 ;
    else
    {
        MAGMA_Z_SET2REAL(res,0);
    }

    A-=ty_* lda  ;
    A-=tx_;

    x= x - break_d *incx  ;
    //x= x - tx * incx ;

    A+=4 * ty* lda  ;

    if( blkc  == ( gridDim.x - 1 ) ) {
        if(tx <= m_mod_thread_x )
            A+=tx;
        else
            A+=m_mod_thread_x;
    }
    else{
        A+=tx;
    }

    magma_int_t wc_c = my_gpu_id ;
    magma_int_t count = 0 ;

    tx_ = td % quarter_thread_x ;
    ty_ = td / quarter_thread_x ;

    WC-=tx ;
    WC+=tx_;


    magma_int_t num_blocks_iters = (blkc +1) /num_gpus - flag;
    
    if( my_gpu_id < ( (blkc+1) % num_gpus) )
    {
    num_blocks_iters += 1;
    }

    x += (my_gpu_id ) * nb ;


    if( blkc > my_gpu_id)

    for(magma_int_t s=0; s<num_blocks_iters; s++)
        {
            MAGMA_Z_SET2REAL(res_,0);
            count++;
            
            if( ty == 0)
            {
                  buff2[tx] = x[0];
            }
                
            __syncthreads();

            #pragma unroll
            for( magma_int_t k=0;k<4;k++){
                #pragma unroll
                for(magma_int_t j=0; j < 4 ; j++)
                    tr[j] = A[j*lda] ;

                #pragma unroll
                for(magma_int_t j=0; j < 4 ; j++){
                    res+=tr[j] * buff2[ quarter_thread_x*k + ty*4+(j)];
                    la[( (j)+ty*4)][tx] = cuConj(tr[j]) * buff[tx];
                }
                __syncthreads();

                MAGMA_Z_SET2REAL(res_, 0) ;

                #pragma unroll
                for(magma_int_t j=0; j < 4 ; j++)
                    res_+=la[tx_][ty_*4+j] ;
                b[k] = res_ ;
                __syncthreads();
                A+=lda* quarter_thread_x ;
            }

            #pragma unroll
            for(magma_int_t k=0; k < 4 ; k++){
                la[tx_][ty_+quarter_thread_x*k]= b[k] ;
            }
            __syncthreads();

            if( ty_ < 4 ) {
                magma_int_t k = ty_*quarter_thread_x;
                res_ = la[tx_][0+k] + la[tx_][1+k] 
                    +  la[tx_][2+k] + la[tx_][3+k]
                    +  la[tx_][4+k] + la[tx_][5+k]
                    +  la[tx_][6+k] + la[tx_][7+k]
                    +  la[tx_][8+k] + la[tx_][9+k]
                    +  la[tx_][10+k]+ la[tx_][11+k]
                    +  la[tx_][12+k]+ la[tx_][13+k]
                    +  la[tx_][14+k]+ la[tx_][15+k];
                WC[k + wc_c*lda ] =   res_;
            }
            
            wc_c += num_gpus;
        x += num_gpus * nb;

            __syncthreads();
        }


    WC+=tx ;
    WC-=tx_;
    la[ty][tx]= res ;
    __syncthreads();

    if( ty == 0 ) {
        res=la[0][tx]+ la[1][tx]+ la[2][tx]+ la[3][tx] ;
        WC[0+lda*(blkc)] = res;
    }
}

__global__ void
magmablas_zhemv_200_L_update_mgpu(magma_int_t n, cuDoubleComplex alpha,
                         cuDoubleComplex* A, magma_int_t lda,
                         cuDoubleComplex *x, magma_int_t incx,
                         cuDoubleComplex beta,
                         cuDoubleComplex *y, magma_int_t incy,
                         cuDoubleComplex *WC,
                         magma_int_t my_gpu_id,
                         magma_int_t num_gpus, 
                         magma_int_t nb )
{
    magma_int_t i;
    magma_int_t tx  = threadIdx.x ;
    magma_int_t ind = blockIdx.x * thread_x + tx ;
    cuDoubleComplex Ca;

    MAGMA_Z_SET2REAL(Ca, 0) ;
    WC+= ind + lda * blockIdx.x;

    for(i = blockIdx.x*thread_x; i<n; i+=thread_x){
        Ca += WC[0] ;
        WC += thread_x;
    }
    if( ind < n )
        y[ind * incy] = beta * y[ind * incy]  + alpha * Ca ;
}


extern "C"
void magmablas_zhemv_200_L_mgpu(magma_int_t m, cuDoubleComplex alpha,
                           cuDoubleComplex *A, magma_int_t lda,
                           cuDoubleComplex *X, magma_int_t incx,
                           cuDoubleComplex beta,
                           cuDoubleComplex *Y, magma_int_t incy,
                           cuDoubleComplex *dC_work,
                         magma_int_t my_gpu_id,
                         magma_int_t num_gpus, 
                         magma_int_t nb)
{
    magma_int_t blocks;

    if (m % zhemv_bs==0)
        blocks = m / zhemv_bs;
    else
        blocks = m / zhemv_bs + 1;

    dim3 grid(blocks, 1, 1);
    dim3 threads(thread_x, thread_y, 1);
    dim3 threads_u(zhemv_bs, 1, 1);

    /*
     * If matrix size is multiple of zhemv_bs, we use a specific code.
     * otherwise, we call the generic case.
     */
    if(m % zhemv_bs == 0 ) {
        magmablas_zhemv_200_L_special_mgpu <<< grid, threads, 0, magma_stream >>>(
            m, alpha, A, lda, X, incx, beta, Y, incy, dC_work, my_gpu_id, num_gpus, nb);
    }
    else{
        magma_int_t m_mod_thread_x = m%zhemv_bs - 1;
        magmablas_zhemv_200_L_generic_mgpu <<< grid, threads, 0, magma_stream >>> (
            m, alpha, A, lda, X, incx ,beta, Y, incy, dC_work, m_mod_thread_x, my_gpu_id, num_gpus, nb);
    }

    magmablas_zhemv_200_L_update_mgpu<<< grid, threads_u, 0, magma_stream >>>(
        m, alpha, A, lda, X, incx, beta, Y, incy, dC_work, my_gpu_id, num_gpus, nb);
}

/*************************************************************************

    Purpose
    =======

    magmablas_zhemv  performs the matrix-vector operation on fermi:

       y := alpha*A*x + beta*y,

    where alpha and beta are scalars, x and y are n element vectors and
    A is an n by n hermitian matrix.

    Arguments
    ==========

    UPLO   - CHARACTER*1.
             On entry, UPLO specifies whether the upper or lower
             triangular part of the array A is to be referenced as
             follows:

                UPLO = 'U' or 'u'   Only the upper triangular part of A
                                    is to be referenced.

                UPLO = 'L' or 'l'   Only the lower triangular part of A
                                    is to be referenced.

             Unchanged on exit.

    N      - INTEGER.
             On entry, N specifies the order of the matrix A.
             N must be at least zero.
             Unchanged on exit.

    ALPHA  - COMPLEX*16      .
             On entry, ALPHA specifies the scalar alpha.
             Unchanged on exit.

    A      - COMPLEX*16       array of DIMENSION ( LDA, n ).
             Before entry with  UPLO = 'U' or 'u', the leading n by n
             upper triangular part of the array A must contain the upper
             triangular part of the hermitian matrix and the strictly
             lower triangular part of A is not referenced.
             Before entry with UPLO = 'L' or 'l', the leading n by n
             lower triangular part of the array A must contain the lower
             triangular part of the hermitian matrix and the strictly
             upper triangular part of A is not referenced.
             Note that the imaginary parts of the diagonal elements need
             not be set and are assumed to be zero.
             Unchanged on exit.

    LDA    - INTEGER.
             On entry, LDA specifies the first dimension of A as declared
             in the calling (sub) program. LDA must be at least
             max( 1, n ).
             Unchanged on exit.
             It is recommended that lda is multiple of 16. Otherwise
             performance would be deteriorated as the memory accesses
             would not be fully coalescent.

    X      - COMPLEX*16       array of dimension at least
             ( 1 + ( n - 1 )*abs( INCX ) ).
             Before entry, the incremented array X must contain the n
             element vector x.
             Unchanged on exit.

    INCX   - INTEGER.
             On entry, INCX specifies the increment for the elements of
             X. INCX must not be zero.
             Unchanged on exit.

    BETA   - COMPLEX*16      .
             On entry, BETA specifies the scalar beta. When BETA is
             supplied as zero then Y need not be set on input.
             Unchanged on exit.

    Y      - COMPLEX*16       array of dimension at least
             ( 1 + ( n - 1 )*abs( INCY ) ).
             Before entry, the incremented array Y must contain the n
             element vector y. On exit, Y is overwritten by the updated
             vector y.

    INCY   - INTEGER.
             On entry, INCY specifies the increment for the elements of
             Y. INCY must not be zero.
             Unchanged on exit.

*/


extern "C"
magma_int_t
magmablas_zhemv_200_mgpu_64( char uplo, magma_int_t n,
                      cuDoubleComplex alpha,
                      cuDoubleComplex **A, magma_int_t lda,
                      cuDoubleComplex **X, magma_int_t incx,
                      cuDoubleComplex beta,
                      cuDoubleComplex **Y, magma_int_t incy,
                      cuDoubleComplex **work, magma_int_t lwork,
              magma_int_t num_gpus, 
              magma_int_t nb)

{
    char uplo_[2] = {uplo, 0};
    int  upper    = lapackf77_lsame(uplo_, "U");

    /*
     * Test the input parameters.
     */
    if ((! upper) && (! lapackf77_lsame(uplo_, "L"))) {
        return -1;
    } else if ( n < 0 ) {
        return -2;
    } else if ( lda < max(1,n) ) {
        return -5;
    } else if ( incx == 0 ) {
        return -7;
    } else if ( incy == 0 ) {
        return -10;
    }

    /*
     * Quick return if possible.
     */
    if ( (n == 0) || ( MAGMA_Z_EQUAL(alpha, MAGMA_Z_ZERO) && MAGMA_Z_EQUAL(beta, MAGMA_Z_ONE) ) )
        return MAGMA_SUCCESS;

    /* TODO: Upper case is not implemented in MAGMA */
    if ( upper)
    {
       printf("Upper case is not implemented on multi GPUs\n");
       exit(0);
    }
    else
    {
        magma_int_t blocks    = n / thread_x + (n % thread_x != 0);
        magma_int_t workspace = lda * (blocks + 1);

        if (lwork < workspace){
           printf("Not enough work space in magmablas_zhemv: passed %d, required %d\n",
                  lwork, workspace);
           exit(1);
        }
        if(nb != 64)
        {
        printf("Error in magmablas_zsymv_200_mgpu: nb != 64, program will exit! please reallocate your matrix among GPUs\n");
        exit(0);
        }
        magma_int_t i = 0;
        for(i=0; i<num_gpus; i++)
        {
             magma_setdevice(i);
             magmablas_zhemv_200_L_mgpu(n, alpha, A[i], lda, X[i], incx, beta, Y[i], incy, work[i], i, num_gpus, nb);     
        }
    }
    return MAGMA_SUCCESS;
}


__global__ void 
kernel_fillZero(cuDoubleComplex *A, magma_int_t size)
{
    magma_int_t id = blockIdx.x * blockDim.x + threadIdx.x;
    if(id < size)
    {
        MAGMA_Z_SET2REAL(A[id], 0.0); 
    }
}


void fillZero(cuDoubleComplex *A, magma_int_t size)
{

    magma_int_t blocks = (size-1)/512 + 1;
    
    dim3 grid(blocks, 1, 1);
    dim3 threads(512, 1, 1);
    
    kernel_fillZero<<<grid, threads>>>(A, size);

}


#endif /* (GPUSHMEM >= 200) */
