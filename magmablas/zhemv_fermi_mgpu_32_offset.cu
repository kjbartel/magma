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


#if (GPUSHMEM >= 200)

#define magmablas_zhemv_200_mgpu_32_offset magmablas_zhemv_mgpu_32_offset
#define magmablas_zhemv2_200_mgpu_32_offset magmablas_zhemv2_mgpu_32_offset

#define zhemv_bs         32
#define thread_x         64
#define thread_y          4
#define bank_shift       33
#define quarter_thread_x 16
#define half_thread_x    32

/*******************************************************************************
 *     Functions for each specific cases - Lower case
 */

#define SWITCH  1400



__global__ void
magmablas_zhemv_200_L_special_mgpu_32_offset_s( magma_int_t n, cuDoubleComplex alpha,
                               cuDoubleComplex *A, magma_int_t lda,
                               cuDoubleComplex *x, magma_int_t incx,
                               cuDoubleComplex  beta,
                               cuDoubleComplex *y, magma_int_t incy,
                               cuDoubleComplex *WC, 
                         magma_int_t my_gpu_id,
                         magma_int_t num_gpus, 
                         magma_int_t nb,
                         magma_int_t kstan)
{


    if(blockIdx.y > blockIdx.x) return;
    if(blockIdx.y % num_gpus != my_gpu_id) return;
    
    if(blockIdx.x < my_gpu_id)
    {
      return;
    }

    magma_int_t tx   = threadIdx.x ;
    magma_int_t ty   = threadIdx.y ;


    cuDoubleComplex res  = MAGMA_Z_ZERO;// used in scan the row
    cuDoubleComplex res_ = MAGMA_Z_ZERO;// used in scan the column

    __shared__ cuDoubleComplex la   [1056];
    __shared__ cuDoubleComplex buff [32];
    __shared__ cuDoubleComplex buff2 [32];


    magma_int_t break_d   =  32 * blockIdx.x;
                  
    A +=  break_d ;
    A +=  lda * ty + tx;
    A +=  lda * (blockIdx.y / num_gpus) * 32; // 

    x +=  tx;

/* ===========================================
    --------------------> threadIdx.y, blockIdx.y
    |
    |
    |
    |
    threadIdx.x blockIdx.x

=========================================== */

    if ( blockIdx.x == blockIdx.y && (blockIdx.x % num_gpus) == my_gpu_id) // diagonal 
    {    
    
        x  += (blockIdx.y * 32) * incx;    
        if( ty == 0 )
        {
            buff[tx] = x[0];
            if(blockIdx.y == 0 && my_gpu_id == 0 && tx < kstan)
            {
                 MAGMA_Z_SET2REAL(buff[tx], 0.0);
            }
        } // obtain the vector x store in buff;
       
        #pragma unroll
        for(magma_int_t j =0; j<half_thread_x; j +=8)
        la[ bank_shift * (ty+j) + tx] =  A[ j * lda];
        __syncthreads();

        #pragma unroll
        for(magma_int_t  i=ty*4; i<(ty * 4 + 4)  ; i++)
        {
            if ( i < tx )   
            {
                la[bank_shift * tx + i] = cuConj(la[ i * bank_shift + tx])  ;
            }
        }
        __syncthreads();

        #pragma unroll
        for(magma_int_t j=0; j < 4 ; j++)
            res += cuConj(la[bank_shift * tx + j + ty * 4])  * buff[j + ty * 4];
    
        __syncthreads();
    
    }
    else // non diagonal
    {
        x  += (blockIdx.x * 32) * incx;    
        if( ty == 0 )
        {
            buff[tx] = x[0];

        } // obtain the vector x and  store in buff; buff store its corresponding upper elements instead of buff2; 
            
        x  -= (blockIdx.x * 32 ) * incx;
        
        x  += (blockIdx.y * 32 ) * incx;    
            
        if( ty == 0 )
        {
            buff2[tx] = x[0]; 
            if(blockIdx.y == 0 && my_gpu_id == 0 && tx < kstan)
            {
                 MAGMA_Z_SET2REAL(buff2[tx], 0.0);
            }
        } // obtain the vector x store in buff2; 
    
        #pragma unroll
        for(magma_int_t j =0; j<half_thread_x; j +=8)
        {
            la[ bank_shift * (ty+j) + tx] =  A[ j * lda];
        }
            
        __syncthreads();

            #pragma unroll
            for(magma_int_t j=0; j < 4 ; j++)
            {
                    res += (la[bank_shift * (ty + j * 8) + tx] )* buff2[ ty + j * 8];
                    res_ += cuConj(la[bank_shift * tx + j + ty * 4]) * buff[j + ty * 4]; // 
            }
            __syncthreads();

            
            la[bank_shift*tx+ty]= res_ ;
            __syncthreads();

            if( ty== 0 )
            {
                  res_ = la[tx*bank_shift+0]+la[tx*bank_shift+1]
                +    la[tx*bank_shift+2]+la[tx*bank_shift+3]
                +    la[tx*bank_shift+4]+la[tx*bank_shift+5]
                +    la[tx*bank_shift+6]+la[tx*bank_shift+7];
            
                WC[ tx + blockIdx.y * 32 + lda * blockIdx.x ] =   res_; // write to its corresponding upper side position
            }
            __syncthreads();
            
      } // end if else 

        la[bank_shift*tx+ty]= res ;
        __syncthreads();

            if( ty== 0 )
            {
                  res = la[tx*bank_shift+0]+la[tx*bank_shift+1]
                +    la[tx*bank_shift+2]+la[tx*bank_shift+3]
                +    la[tx*bank_shift+4]+la[tx*bank_shift+5]
                +    la[tx*bank_shift+6]+la[tx*bank_shift+7];
            
                 WC[ tx + blockIdx.x * 32 + lda * blockIdx.y] =  res;
            }
        
}



__global__ void
magmablas_zhemv_200_L_generic_mgpu_32_offset_s(magma_int_t n, cuDoubleComplex alpha,
                              cuDoubleComplex *A, magma_int_t lda,
                              cuDoubleComplex *x, magma_int_t incx,
                              cuDoubleComplex beta,
                              cuDoubleComplex *y, magma_int_t incy,
                              cuDoubleComplex *WC,
                              magma_int_t m_mod_thread_x,
                         magma_int_t my_gpu_id,
                         magma_int_t num_gpus, 
                         magma_int_t nb,
                         magma_int_t kstan)
{


    if(blockIdx.y > blockIdx.x) return;
    if(blockIdx.y % num_gpus != my_gpu_id) return;
    
    if(blockIdx.x < my_gpu_id)
    {
      return;
    }

    magma_int_t tx   = threadIdx.x ;
    magma_int_t ty   = threadIdx.y ;


    cuDoubleComplex res  = MAGMA_Z_ZERO;// used in scan the row
    cuDoubleComplex res_ = MAGMA_Z_ZERO;// used in scan the column

    __shared__ cuDoubleComplex la   [1056];
    __shared__ cuDoubleComplex buff [32];
    __shared__ cuDoubleComplex buff2 [32];


    magma_int_t break_d   =  32 * blockIdx.x;
                  
    A +=  break_d ;
    A +=  lda * ty;
    A +=  lda * (blockIdx.y / num_gpus) * 32; // 
    x +=  tx;
    x  += (blockIdx.x * 32) * incx;    

    magma_int_t trackA ;
    if( blockIdx.x == ( gridDim.x - 1 ) ) {
        if( ty == 0 ){
            if( tx > m_mod_thread_x )
            {
                MAGMA_Z_SET2REAL(buff[tx],0);
            }
            else
                buff[tx]  = x[0];
        }
        if ( tx > m_mod_thread_x )
            trackA=m_mod_thread_x;
        else
            trackA=tx;
        A += trackA ;
    }
    else {
        if( ty == 0 ){
            buff[tx]  = x[0];
        }
        trackA = tx;
        A += trackA ;
    }

    //__syncthreads();

    if ( blockIdx.x == blockIdx.y && (blockIdx.x % num_gpus) == my_gpu_id) // diagonal 
    {    
    
        if( ty == 0 )
        {
            if(blockIdx.y == 0 && my_gpu_id == 0 && tx < kstan)
            {
                 MAGMA_Z_SET2REAL(buff[tx], 0.0);
            }
        } // obtain the vector x store in buff;
       
               if( blockIdx.x == ( gridDim.x - 1 ) ) {
        #pragma unroll
        for(magma_int_t j =0; j<half_thread_x; j+=8){
            if( ( ty + j ) > m_mod_thread_x )
            {
                MAGMA_Z_SET2REAL(la[bank_shift*(ty+j)+tx], 9999);
            }
            else
                la[bank_shift*(ty+j)+tx] =  A[ j * lda];
        }
        }
        else {
        #pragma unroll
        for(magma_int_t j =0; j<half_thread_x; j+=8){
            la[bank_shift*(ty+j)+tx] = A[ j * lda];
        }
        }
        __syncthreads();


        #pragma unroll
        for(magma_int_t  i=ty*4; i<(ty * 4 + 4)  ; i++)
        {
            if ( i < tx )   
            {
                la[bank_shift * tx + i] = cuConj(la[ i * bank_shift + tx])  ;
            }
        }
        __syncthreads();

        #pragma unroll
        for(magma_int_t j=0; j < 4 ; j++)
            res += cuConj(la[bank_shift * tx + j + ty * 4])  * buff[j + ty * 4];
    
        __syncthreads();
    
    }
    else // non diagonal
    {
            
            
        x  -= (blockIdx.x * 32 ) * incx;
        
        x  += (blockIdx.y * 32 ) * incx;    
            
        if( ty == 0 )
        {
            buff2[tx] = x[0]; 
            if(blockIdx.y == 0 && my_gpu_id == 0 && tx < kstan)
            {
                 MAGMA_Z_SET2REAL(buff2[tx], 0.0);
            }
        } // obtain the vector x store in buff2; 
    
        #pragma unroll
        for(magma_int_t j =0; j<half_thread_x; j +=8)
        {
            la[ bank_shift * (ty+j) + tx] =  A[ j * lda];
        }
            
        __syncthreads();

            #pragma unroll
            for(magma_int_t j=0; j < 4 ; j++)
            {
                    res += (la[bank_shift * (ty + j * 8) + tx] )* buff2[ ty + j * 8];
                    res_ += cuConj(la[bank_shift * tx + j + ty * 4]) * buff[j + ty * 4]; // 
            }
            __syncthreads();

            
            la[bank_shift*tx+ty]= res_ ;
            __syncthreads();

            if( ty== 0 )
            {
                  res_ = la[tx*bank_shift+0]+la[tx*bank_shift+1]
                +    la[tx*bank_shift+2]+la[tx*bank_shift+3]
                +    la[tx*bank_shift+4]+la[tx*bank_shift+5]
                +    la[tx*bank_shift+6]+la[tx*bank_shift+7];
            
                WC[ tx + blockIdx.y * 32 + lda * blockIdx.x ] =   res_; // write to its corresponding upper side position
            }
            __syncthreads();
            
      } // end if else 

        la[bank_shift*tx+ty]= res ;
        __syncthreads();

            if( ty== 0 )
            {
                  res = la[tx*bank_shift+0]+la[tx*bank_shift+1]
                +    la[tx*bank_shift+2]+la[tx*bank_shift+3]
                +    la[tx*bank_shift+4]+la[tx*bank_shift+5]
                +    la[tx*bank_shift+6]+la[tx*bank_shift+7];
            
                 WC[ tx + blockIdx.x * 32 + lda * blockIdx.y] =  res;
            }
        
}



__global__ void
magmablas_zhemv_200_L_special_mgpu_32_offset( magma_int_t n, cuDoubleComplex alpha,
                               cuDoubleComplex *A, magma_int_t lda,
                               cuDoubleComplex *x, magma_int_t incx,
                               cuDoubleComplex  beta,
                               cuDoubleComplex *y, magma_int_t incy,
                               cuDoubleComplex *WC, 
                         magma_int_t my_gpu_id,
                         magma_int_t num_gpus, 
                         magma_int_t nb,
                         magma_int_t kstan)
{
    magma_int_t tx   = threadIdx.x ;
    magma_int_t ty   = threadIdx.y ;
    magma_int_t blkc = blockIdx.x ;

    if(blkc < my_gpu_id)
    {
    return;
    }

    cuDoubleComplex res  = MAGMA_Z_ZERO;// used in scan the row
    cuDoubleComplex res_ = MAGMA_Z_ZERO;// used in scan the column
    cuDoubleComplex res1 = MAGMA_Z_ZERO;// tem for res
    cuDoubleComplex res2 = MAGMA_Z_ZERO;// tem for res_

    __shared__ cuDoubleComplex la   [16][64+2];
    __shared__ cuDoubleComplex sdata   [32][9];
    __shared__ cuDoubleComplex buff [32];
    __shared__ cuDoubleComplex buff2 [32];


    magma_int_t break_d   =  32 * blkc;

    x  += (break_d + tx ) * incx;
    A  +=  break_d ;
    A  +=  ty * lda + tx ;

    if( ty == 0 )
    {
        buff[tx] = x[0];
    if(blkc == 0 && my_gpu_id == 0 && tx < kstan)
    {
             MAGMA_Z_SET2REAL(buff[tx], 0.0);
        }
    } // obtain the vector x store in buff;
    

    magma_int_t flag = 0;
    
    if ( (blkc % num_gpus) == my_gpu_id) 
    {
        A += lda * (blkc/num_gpus) * 32; // change

        #pragma unroll
        for(magma_int_t j =0; j<half_thread_x; j +=8)
        la[0][ bank_shift * (ty+j) + tx] =  A[ j * lda];
        __syncthreads();

        #pragma unroll
        for(magma_int_t  i=ty*4; i<(ty * 4 + 4)  ; i++){
        if ( i < tx )   {
            la[0][bank_shift * tx + i] = cuConj( la[0][ i * bank_shift + tx] ) ;
        }
        }
        __syncthreads();

        #pragma unroll
        for(magma_int_t j=0; j < 4 ; j++)
            res += cuConj( la[0][bank_shift * tx + j + ty * 4] ) * buff[j + ty * 4];
            __syncthreads();
            

             A -= lda * (blkc/num_gpus) * 32; 
    
              flag = 1;
        }

        

        x -= blkc * 32  *incx  ;

        x= x- tx*incx;

        magma_int_t wc_c = my_gpu_id ;
        magma_int_t count = 0 ;

               WC +=  break_d + tx;
  
        magma_int_t num_blocks_iters = (blkc +1) /num_gpus - flag;
    
        if( my_gpu_id < ( (blkc+1) % num_gpus) )
        {
        num_blocks_iters += 1;
        }

        x += (my_gpu_id ) * 32 ;

        if( blkc > my_gpu_id)

        for(magma_int_t s=0; s<num_blocks_iters; s++)
        {
            MAGMA_Z_SET2REAL(res_,0);
            count++;

                     #pragma unroll
            for(magma_int_t j =0; j<half_thread_x; j +=8)
            la[0][ bank_shift * (ty+j) + tx] =  A[ j * lda];

            if( ty == 0 )
            {
                buff2[tx] = x[tx];
                if(my_gpu_id == 0 && tx < kstan && count==1)
                {
                     MAGMA_Z_SET2REAL(buff2[tx], 0.0);
                }
            } // obtain the vector x store in buff2;
            __syncthreads();

            #pragma unroll
            for(magma_int_t j=0; j < 4 ; j++)
            {
                        res += (la[0][bank_shift * (ty + j * 8) + tx] )* buff2[ ty + j * 8];
                        res_ += cuConj( la[0][bank_shift * tx + j + ty * 4] ) * buff[j + ty * 4]; //iterate colum
                }
                     //__syncthreads();

                    //la[0][bank_shift*tx+ty]= res_ ;
                    sdata[tx][ty]= res_ ;
            __syncthreads();
/*
            if( ty== 0 )
            {
              res2 = la[0][tx*bank_shift+0]+la[0][tx*bank_shift+1]
            +    la[0][tx*bank_shift+2]+la[0][tx*bank_shift+3]
            +    la[0][tx*bank_shift+4]+la[0][tx*bank_shift+5]
            +    la[0][tx*bank_shift+6]+la[0][tx*bank_shift+7];
                        
                WC[wc_c*lda ] =   res2;
            }
            __syncthreads();

*/ 
            if( ty== 1 )
            {
              res2 = sdata[tx][0]+sdata[tx][1]
              + sdata[tx][2]+sdata[tx][3]
              + sdata[tx][4]+sdata[tx][5]
              + sdata[tx][6]+sdata[tx][7];
                        
                WC[wc_c*lda ] =   res2;
            }


                    wc_c += num_gpus;
                x += num_gpus * 32;
                    A += lda * 32 ;

               }


        la[0][bank_shift*tx+ty]= res ;
        __syncthreads();

            if( ty== 0 )
            {
              res1 = la[0][tx*bank_shift+0]+la[0][tx*bank_shift+1]
            +    la[0][tx*bank_shift+2]+la[0][tx*bank_shift+3]
            +    la[0][tx*bank_shift+4]+la[0][tx*bank_shift+5]
            +    la[0][tx*bank_shift+6]+la[0][tx*bank_shift+7];
              
              WC[0+lda*(blkc)] =  res1;
            }
        
}

/*******************************************************************************
 *     Functions for each specific cases - Upper case
 */


__global__ void
magmablas_zhemv_200_U_special_mgpu_32_offset( magma_int_t n, cuDoubleComplex alpha,
                               cuDoubleComplex *A, magma_int_t lda,
                               cuDoubleComplex *x, magma_int_t incx,
                               cuDoubleComplex  beta,
                               cuDoubleComplex *y, magma_int_t incy,
                               cuDoubleComplex *WC, 
                         magma_int_t my_gpu_id,
                         magma_int_t num_gpus, 
                         magma_int_t nb,
                         magma_int_t kstan)
{
    magma_int_t tx   = threadIdx.x ;
    magma_int_t ty   = threadIdx.y ;
    magma_int_t blkc = blockIdx.x ;


    cuDoubleComplex res  = MAGMA_Z_ZERO;// used in scan the row
    cuDoubleComplex res_ = MAGMA_Z_ZERO;// used in scan the column
    cuDoubleComplex res1 = MAGMA_Z_ZERO;// tem for res
    cuDoubleComplex res2 = MAGMA_Z_ZERO;// tem for res_

    __shared__ cuDoubleComplex la   [16][66];
    __shared__ cuDoubleComplex buff [32];
    __shared__ cuDoubleComplex buff2 [32];


    magma_int_t break_d   =  32 * blkc;

    x  += (break_d + tx ) * incx;
    A  +=  break_d ;
    A  +=  ty * lda + tx ;

    if( ty == 0 )
    {
        buff[tx] = x[0];
    if(blkc == 0  && tx < kstan)
    {
             MAGMA_Z_SET2REAL(buff[tx], 0.0);
        }
    } // obtain the vector x store in buff;
 
    
    if ( (blkc % num_gpus) == my_gpu_id) 
    {
        A += lda * (blkc/num_gpus) * 32; // change

        #pragma unroll
        for(magma_int_t j =0; j<half_thread_x; j +=8)
        la[0][ bank_shift * (ty+j) + tx] =  A[ j * lda];
        __syncthreads();

        #pragma unroll
        for(magma_int_t  i=ty*4; i<(ty * 4 + 4)  ; i++){
            if ( i > tx )   
            {
                la[0][bank_shift * tx + i] = cuConj(la[0][ i * bank_shift + tx])  ;
            }
        }
        __syncthreads();

        #pragma unroll
        for(magma_int_t j=0; j < 4 ; j++)
            res += cuConj( la[0][bank_shift * tx + j + ty * 4] ) * buff[j + ty * 4];
    
            __syncthreads();

             A -= lda * (blkc/num_gpus) * 32; 
    
              
        }
        __syncthreads();

    
        x  -= (break_d + tx ) * incx;// return to the beginning

            x += (my_gpu_id ) * 32 ;// 

            magma_int_t wc_c = my_gpu_id ;

        magma_int_t total_blocks_gpu = gridDim.x /num_gpus;

        if( my_gpu_id < ( gridDim.x % num_gpus) )
        {
        total_blocks_gpu += 1;
        }

        magma_int_t shift = (blkc +1) /num_gpus ;
    
        if( my_gpu_id < ( (blkc+1) % num_gpus) )
        {
        shift += 1;
        }

            #pragma unroll
            for(magma_int_t s=0; s<shift; s++)
            {
                x += num_gpus * 32;
                    A += lda * 32 ;
            wc_c += num_gpus;
            }


               WC +=  break_d + tx;
           
           magma_int_t num_blocks_iters = total_blocks_gpu - shift;

       magma_int_t count = 0;


        for(magma_int_t s=0; s<num_blocks_iters; s++)
        {
            MAGMA_Z_SET2REAL(res_,0);
            count++;

                     #pragma unroll
            for(magma_int_t j =0; j<half_thread_x; j +=8)
            la[0][ bank_shift * (ty+j) + tx] =  A[ j * lda];

            if( ty == 0 )
            {
                buff2[tx] = x[tx];

            } // obtain the vector x store in buff;
            __syncthreads();

            #pragma unroll
            for(magma_int_t j=0; j < 4 ; j++)
            {
                        res += (la[0][bank_shift * (ty + j * 8) + tx] )* buff2[ ty + j * 8];
                        res_ += cuConj( la[0][bank_shift * tx + j + ty * 4] ) * buff[j + ty * 4]; //iterate colum
                }
                     __syncthreads();

                    la[0][bank_shift*tx+ty]= res_ ;
            __syncthreads();

            if( ty== 0 )
            {
              res2 = la[0][tx*bank_shift+0]+la[0][tx*bank_shift+1]
            +    la[0][tx*bank_shift+2]+la[0][tx*bank_shift+3]
            +    la[0][tx*bank_shift+4]+la[0][tx*bank_shift+5]
            +    la[0][tx*bank_shift+6]+la[0][tx*bank_shift+7];
                        
                WC[wc_c*lda ] =   res2;
            }

 
            __syncthreads();


                    wc_c += num_gpus;
                x += num_gpus * 32;
                    A += lda * 32 ;

               }


        la[0][bank_shift*tx+ty]= res ;
        __syncthreads();

        if( ty== 0 )
            {
              res1 = la[0][tx*bank_shift+0]+la[0][tx*bank_shift+1]
            +    la[0][tx*bank_shift+2]+la[0][tx*bank_shift+3]
            +    la[0][tx*bank_shift+4]+la[0][tx*bank_shift+5]
            +    la[0][tx*bank_shift+6]+la[0][tx*bank_shift+7];
              
              WC[0+lda*(blkc)] =  res1;
        }
}

/**************************************************************
 *    Lower case for generic sizes
 */
__global__ void
magmablas_zhemv_200_L_generic_mgpu_32_offset(magma_int_t n, cuDoubleComplex alpha,
                              cuDoubleComplex *A, magma_int_t lda,
                              cuDoubleComplex *x, magma_int_t incx,
                              cuDoubleComplex beta,
                              cuDoubleComplex *y, magma_int_t incy,
                              cuDoubleComplex *WC,
                              magma_int_t m_mod_thread_x,
                         magma_int_t my_gpu_id,
                         magma_int_t num_gpus, 
                         magma_int_t nb,
                         magma_int_t kstan)
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

    __shared__ cuDoubleComplex la   [16][64+2];
    __shared__ cuDoubleComplex sdata   [32][9];
    __shared__ cuDoubleComplex buff [32];
    __shared__ cuDoubleComplex buff2 [32];


    magma_int_t break_d   =  32 * blkc;

    x += (break_d + tx ) * incx;
    A +=  break_d ;
    A += lda * ty;

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
        if ( tx > m_mod_thread_x )
            trackA=m_mod_thread_x;
        else
            trackA=tx;
        A += trackA ;
    }
    else {
        if( ty == 0 ){
            buff[tx]  = x[0];
        }
        trackA = tx;
        A += trackA ;
    }

     if(ty == 0 )
     { 
           if(my_gpu_id == 0 && blkc ==0  && tx < kstan)//
       {
                 MAGMA_Z_SET2REAL(buff[tx], 0.0);
       }
     }

    magma_int_t flag = 0;
    
    if ( (blkc % num_gpus) == my_gpu_id) 
    {
        A += lda * (blkc/num_gpus) * 32; // change
        // Somehow merging these two if - else creates problem
        // It could be a potential bug -- from synchronization or from cuda or compiler
        if( blkc == ( gridDim.x - 1 ) ) {
        #pragma unroll
        for(magma_int_t j =0; j<half_thread_x; j+=8){
            if( ( ty + j ) > m_mod_thread_x )
            {
                MAGMA_Z_SET2REAL(la[0][bank_shift*(ty+j)+tx], 9999);
            }
            else
                la[0][bank_shift*(ty+j)+tx] =  A[ j * lda];
        }
        }
        else {
        #pragma unroll
        for(magma_int_t j =0; j<half_thread_x; j+=8){
            la[0][bank_shift*(ty+j)+tx] = A[ j * lda];
        }
        }
        __syncthreads();

        #pragma unroll
        for(magma_int_t  i=ty*4; i<(ty*4+4)  ; i++){
        if ( i < tx )   {
            la[0][bank_shift*tx+i] = cuConj(la[0][i*bank_shift+tx]) ;
        }
        }
        __syncthreads();

        #pragma unroll
        for(magma_int_t j=0; j < 4 ; j++)
        res += cuConj(la[0][bank_shift*tx+j+ty*4])* buff[j+ty*4];
        __syncthreads();

       
          A -= lda * (blkc/num_gpus) * 32; 
    
          flag = 1;
    }

    __syncthreads();


    x= x - break_d *incx  ;
    x= x - tx * incx ;


    magma_int_t wc_c = my_gpu_id ;
    magma_int_t count = 0 ;

    WC +=  break_d + tx;


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

                     #pragma unroll
            for(magma_int_t j =0; j<half_thread_x; j +=8)
            la[0][ bank_shift * (ty+j) + tx] =  A[ j * lda];

            if( ty == 0 )
            {
                buff2[tx] = x[tx];
                if(my_gpu_id == 0 && tx < kstan && count==1)//
                {
                     MAGMA_Z_SET2REAL(buff2[tx], 0.0);
                }
            } // obtain the vector x store in buff2;
            __syncthreads();

            #pragma unroll
            for(magma_int_t j=0; j < 4 ; j++)
            {
            
                        res += (la[0][bank_shift * (ty + j * 8) + tx] )* buff2[ ty + j * 8];
                        res_ += cuConj( la[0][bank_shift * tx + j + ty * 4] ) * buff[j + ty * 4]; //iterate colum
                }
                   //  __syncthreads();

                    //la[0][bank_shift*tx+ty]= res_ ;
                   sdata[tx][ty]= res_ ;
            __syncthreads();
/*
            if( ty== 0 )
            {
              res2 = la[0][tx*bank_shift+0]+la[0][tx*bank_shift+1]
            +    la[0][tx*bank_shift+2]+la[0][tx*bank_shift+3]
            +    la[0][tx*bank_shift+4]+la[0][tx*bank_shift+5]
            +    la[0][tx*bank_shift+6]+la[0][tx*bank_shift+7];
                        
                WC[wc_c*lda ] =   res2;
            }
 
            __syncthreads();
*/

            if( ty== 1 )
            {
              res2 = sdata[tx][0]+sdata[tx][1]
              + sdata[tx][2]+sdata[tx][3]
              + sdata[tx][4]+sdata[tx][5]
              + sdata[tx][6]+sdata[tx][7];
                        
                WC[wc_c*lda ] =   res2;
            }


                    wc_c += num_gpus;
                x += num_gpus * 32;
                    A += lda * 32 ;


               }


        la[0][bank_shift*tx+ty]= res ;
        __syncthreads();

            if( ty== 0 )
            {
              res1 = la[0][tx*bank_shift+0]+la[0][tx*bank_shift+1]
            +    la[0][tx*bank_shift+2]+la[0][tx*bank_shift+3]
            +    la[0][tx*bank_shift+4]+la[0][tx*bank_shift+5]
            +    la[0][tx*bank_shift+6]+la[0][tx*bank_shift+7];
              WC[0+lda*(blkc)] =  res1;
            }

}



__global__ void
magmablas_zhemv_200_U_generic_mgpu_32_offset(magma_int_t n, cuDoubleComplex alpha,
                              cuDoubleComplex *A, magma_int_t lda,
                              cuDoubleComplex *x, magma_int_t incx,
                              cuDoubleComplex beta,
                              cuDoubleComplex *y, magma_int_t incy,
                              cuDoubleComplex *WC,
                              magma_int_t m_mod_thread_x,
                         magma_int_t my_gpu_id,
                         magma_int_t num_gpus, 
                         magma_int_t nb,
                         magma_int_t kstan,
                                                 magma_int_t the_right_gpu)
{
    magma_int_t tx   = threadIdx.x ;
    magma_int_t ty   = threadIdx.y ;
    magma_int_t blkc = blockIdx.x ;


    cuDoubleComplex res  = MAGMA_Z_ZERO;
    cuDoubleComplex res_ = MAGMA_Z_ZERO;
    cuDoubleComplex res1 = MAGMA_Z_ZERO;
    cuDoubleComplex res2 = MAGMA_Z_ZERO;

    __shared__ cuDoubleComplex la   [16][64+2];
    __shared__ cuDoubleComplex buff [32];
    __shared__ cuDoubleComplex buff2 [32];


    magma_int_t break_d   =  32 * blkc;

    x += (break_d + tx ) * incx;
    A +=  break_d ;
    A += lda * ty;

    magma_int_t trackA ;
    if( blkc == ( gridDim.x - 1 )) 
    {
        if( ty == 0 ){
            if( tx > m_mod_thread_x )
            {
                MAGMA_Z_SET2REAL(buff[tx],0);
            }
            else
                buff[tx]  = x[0];
        }
        if ( tx > m_mod_thread_x )
            trackA=m_mod_thread_x;
        else
            trackA=tx;
        A += trackA ;
    }
    else
    {
        if( ty == 0 )
        {
            buff[tx]  = x[0];
        }
        
        A += tx ;
    }

     if(ty == 0 )
     { 
           if(blkc ==0  && tx < kstan)//
       {
                 MAGMA_Z_SET2REAL(buff[tx], 0.0);
       }
     }

    
    if ( (blkc % num_gpus) == my_gpu_id) 
    {
        A += lda * (blkc/num_gpus) * 32; // change

        if( blkc == ( gridDim.x - 1 ) ) {
        #pragma unroll
        for(magma_int_t j =0; j<half_thread_x; j+=8){
            if( ( ty + j ) > m_mod_thread_x )
            {
                MAGMA_Z_SET2REAL(la[0][bank_shift*(ty+j)+tx], 9999);
            }
            else
                la[0][bank_shift*(ty+j)+tx] =  A[ j * lda];
        }
        }
        else {
        #pragma unroll
        for(magma_int_t j =0; j<half_thread_x; j+=8){
            la[0][bank_shift*(ty+j)+tx] = A[ j * lda];
        }
        }
        __syncthreads();

        #pragma unroll
        for(magma_int_t  i=ty*4; i<(ty*4+4)  ; i++){
            if ( i > tx )   
            {
                la[0][bank_shift * tx + i] = cuConj(la[0][ i * bank_shift + tx])  ;
            }
        }
        __syncthreads();

        #pragma unroll
        for(magma_int_t j=0; j < 4 ; j++)
        res += cuConj(la[0][bank_shift*tx+j+ty*4])* buff[j+ty*4];
        __syncthreads();

       
          A -= lda * (blkc/num_gpus) * 32; 
    
      
    }


    
        x  -= (break_d + tx ) * incx;// return to the beginning

            x += (my_gpu_id ) * 32 ;// 

            magma_int_t wc_c = my_gpu_id ;

        magma_int_t total_blocks_gpu = gridDim.x /num_gpus;

        if( my_gpu_id < ( gridDim.x % num_gpus) )
        {
        total_blocks_gpu += 1;
        }

        magma_int_t shift = (blkc +1) /num_gpus ;
    
        if( my_gpu_id < ( (blkc+1) % num_gpus) )
        {
        shift += 1;
        }

            #pragma unroll
            for(magma_int_t s=0; s<shift; s++)
            {
                x += num_gpus * 32;
                    A += lda * 32 ;
            wc_c += num_gpus;
            }


               WC +=  break_d + tx;
           
           magma_int_t num_blocks_iters = total_blocks_gpu - shift;

       magma_int_t count = 0;


        for(magma_int_t s=0; s<num_blocks_iters; s++)
        {
            MAGMA_Z_SET2REAL(res_,0);
            count++;


           if(my_gpu_id == the_right_gpu && s==num_blocks_iters-1)
           {
            
                         if( ty == 0 )
                         {
                    if( tx > m_mod_thread_x )
                    {
                    MAGMA_Z_SET2REAL(buff2[tx],0);
                    }
                    else
                         
                    buff2[tx]  = x[tx];
                      }
              
                #pragma unroll
                  for(magma_int_t j =0; j<half_thread_x; j+=8)
                      {
                    if( ( ty + j ) > m_mod_thread_x )
                    {
                    MAGMA_Z_SET2REAL(la[0][bank_shift*(ty+j)+tx], 0);
                    }
                    else
                    la[0][bank_shift*(ty+j)+tx] =  A[ j * lda];
                }
                     __syncthreads();

                 }// end of the_right_gpu
                 else
                 {
                   #pragma unroll
            for(magma_int_t j =0; j<half_thread_x; j +=8)
            la[0][ bank_shift * (ty+j) + tx] =  A[ j * lda];
           
            if( ty == 0 )
            {
            buff2[tx] = x[tx];
            } // obtain the vector x store in buff;
            __syncthreads();
                 }


            #pragma unroll
            for(magma_int_t j=0; j < 4 ; j++)
            {
            
                        res += (la[0][bank_shift * (ty + j * 8) + tx] )* buff2[ ty + j * 8];
                        res_ += cuConj( la[0][bank_shift * tx + j + ty * 4] ) * buff[j + ty * 4]; //iterate colum
                }
                     __syncthreads();

                    la[0][bank_shift*tx+ty]= res_ ;
            __syncthreads();

            if( ty== 0 )
            {
              res2 = la[0][tx*bank_shift+0]+la[0][tx*bank_shift+1]
            +    la[0][tx*bank_shift+2]+la[0][tx*bank_shift+3]
            +    la[0][tx*bank_shift+4]+la[0][tx*bank_shift+5]
            +    la[0][tx*bank_shift+6]+la[0][tx*bank_shift+7];
                        
                WC[wc_c*lda ] =   res2;
            }
 
            __syncthreads();


                    wc_c += num_gpus;
                x += num_gpus * 32;
                    A += lda * 32 ;


               }


        la[0][bank_shift*tx+ty]= res ;
        __syncthreads();

            if( ty== 0 )
            {
              res1 = la[0][tx*bank_shift+0]+la[0][tx*bank_shift+1]
            +    la[0][tx*bank_shift+2]+la[0][tx*bank_shift+3]
            +    la[0][tx*bank_shift+4]+la[0][tx*bank_shift+5]
            +    la[0][tx*bank_shift+6]+la[0][tx*bank_shift+7];
              WC[0+lda*(blkc)] =  res1;
            }

}

/**************************************************************
 *    
 */

__global__ void
magmablas_zhemv_200_L_update_mgpu_32_offset_s(magma_int_t n, cuDoubleComplex alpha,
                         cuDoubleComplex* A, magma_int_t lda,
                         cuDoubleComplex *x, magma_int_t incx,
                         cuDoubleComplex beta,
                         cuDoubleComplex *y, magma_int_t incy,
                         cuDoubleComplex *WC,
                         magma_int_t my_gpu_id,
                         magma_int_t num_gpus, 
                         magma_int_t nb,
                                                 magma_int_t kstan )
{
    magma_int_t i;
    magma_int_t tx  = threadIdx.x ;
    magma_int_t ind = blockIdx.x * 32 + tx ;
    cuDoubleComplex Ca;

    MAGMA_Z_SET2REAL(Ca, 0) ;
    WC+= ind;

    for(i =0; i<n; i+=32){
        Ca += WC[i/32 * lda] ;
    }
    if( ind < n && ind >= kstan)
        y[ind * incy] = beta * y[ind * incy]  + alpha * Ca ;

}




__global__ void
magmablas_zhemv_200_L_update_mgpu_32_offset(magma_int_t n, cuDoubleComplex alpha,
                         cuDoubleComplex* A, magma_int_t lda,
                         cuDoubleComplex *x, magma_int_t incx,
                         cuDoubleComplex beta,
                         cuDoubleComplex *y, magma_int_t incy,
                         cuDoubleComplex *WC,
                         magma_int_t my_gpu_id,
                         magma_int_t num_gpus, 
                         magma_int_t nb,
                                                 magma_int_t kstan )
{
    magma_int_t i;
    magma_int_t tx  = threadIdx.x ;
    magma_int_t ind = blockIdx.x * 32 + tx ;
    cuDoubleComplex Ca;

    MAGMA_Z_SET2REAL(Ca, 0) ;
    WC+= ind + lda * blockIdx.x;

    for(i = blockIdx.x*32; i<n; i+=32){
        Ca += WC[0] ;
        WC += 32;
    }
    if( ind < n && ind >= kstan)
        y[ind * incy] = beta * y[ind * incy]  + alpha * Ca ;
}

__global__ void
magmablas_zhemv_200_U_update_mgpu_32_offset(magma_int_t n, cuDoubleComplex alpha,
                         cuDoubleComplex* A, magma_int_t lda,
                         cuDoubleComplex *x, magma_int_t incx,
                         cuDoubleComplex beta,
                         cuDoubleComplex *y, magma_int_t incy,
                         cuDoubleComplex *WC,
                         magma_int_t my_gpu_id,
                         magma_int_t num_gpus, 
                         magma_int_t nb,
                                                 magma_int_t kstan )
{
    magma_int_t i;
    magma_int_t tx  = threadIdx.x ;
    magma_int_t ind = blockIdx.x * 32 + tx ;
    cuDoubleComplex Ca;

    MAGMA_Z_SET2REAL(Ca, 0) ;
    WC+=  blockIdx.x * lda + tx;

    for(i = 0; i<(blockIdx.x+1)*32; i+=32)
    {

        Ca += WC[0] ;
        WC += 32 ;
    }
    if( ind < n && ind >= kstan)
        y[ind * incy] = beta * y[ind * incy]  + alpha * Ca ;
}



extern "C"
void magmablas_zhemv_200_L_mgpu_32_offset(magma_int_t m, cuDoubleComplex alpha,
                           cuDoubleComplex *A, magma_int_t lda,
                           cuDoubleComplex *X, magma_int_t incx,
                           cuDoubleComplex beta,
                           cuDoubleComplex *Y, magma_int_t incy,
                           cuDoubleComplex *dC_work,
                         magma_int_t my_gpu_id,
                         magma_int_t num_gpus, 
                         magma_int_t nb,
                             magma_int_t offset,
                         magma_int_t num_blocks_skipped)
{

    magma_int_t the_chosen_block_id = offset / 32; 
   
    magma_int_t kstan = offset % 32;

    A += lda * num_blocks_skipped * 32 + the_chosen_block_id * 32;
    X += the_chosen_block_id * 32;
    Y += the_chosen_block_id * 32;

    magma_int_t blocks;

    if (m % zhemv_bs==0)
        blocks = m / zhemv_bs;
    else
        blocks = m / zhemv_bs + 1;

    blocks -= the_chosen_block_id;

    dim3 grid(blocks, 1, 1);
    dim3 grid_s(blocks, blocks, 1);
    dim3 threads(32, 8, 1);
    dim3 threads_u(32, 1, 1);

        
    /*
         * If matrix size is multiple of zhemv_bs, we use a specific code.
         * otherwise, we call the generic case.
         */
      if(m % zhemv_bs == 0 ) 
      {
           if( m < SWITCH)
        magmablas_zhemv_200_L_special_mgpu_32_offset_s <<< grid_s, threads, 0, magma_stream >>>(
            m, alpha, A, lda, X, incx, beta, Y, incy, dC_work, my_gpu_id, num_gpus, nb, kstan);
           else
        magmablas_zhemv_200_L_special_mgpu_32_offset <<< grid, threads, 0, magma_stream >>>(
            m, alpha, A, lda, X, incx, beta, Y, incy, dC_work, my_gpu_id, num_gpus, nb, kstan);
        }
    else
        {
         magma_int_t m_mod_thread_x = m%zhemv_bs - 1;

           if( m  < SWITCH)
        magmablas_zhemv_200_L_generic_mgpu_32_offset_s <<< grid_s, threads, 0, magma_stream >>> (
            m, alpha, A, lda, X, incx ,beta, Y, incy, dC_work, m_mod_thread_x, my_gpu_id, num_gpus, nb, kstan);
           else

        magmablas_zhemv_200_L_generic_mgpu_32_offset <<< grid, threads, 0, magma_stream >>> (
            m, alpha, A, lda, X, incx ,beta, Y, incy, dC_work, m_mod_thread_x, my_gpu_id, num_gpus, nb, kstan);
        }
       if( m < SWITCH)
        magmablas_zhemv_200_L_update_mgpu_32_offset_s<<< grid, threads_u, 0, magma_stream >>>(
        m, alpha, A, lda, X, incx, beta, Y, incy, dC_work, my_gpu_id, num_gpus, nb, kstan);
       else
        magmablas_zhemv_200_L_update_mgpu_32_offset<<< grid, threads_u, 0, magma_stream >>>(
        m, alpha, A, lda, X, incx, beta, Y, incy, dC_work, my_gpu_id, num_gpus, nb, kstan);
    
}

extern "C"
void magmablas_zhemv_200_U_mgpu_32_offset(magma_int_t m, cuDoubleComplex alpha,
                           cuDoubleComplex *A, magma_int_t lda,
                           cuDoubleComplex *X, magma_int_t incx,
                           cuDoubleComplex beta,
                           cuDoubleComplex *Y, magma_int_t incy,
                           cuDoubleComplex *dC_work,
                         magma_int_t my_gpu_id,
                         magma_int_t num_gpus, 
                         magma_int_t nb,
                             magma_int_t offset,
                         magma_int_t num_blocks_skipped, 
                                                 magma_int_t the_right_gpu)
{

    magma_int_t the_chosen_block_id = offset / 32; 
    magma_int_t kstan = offset % 32;

    A += lda * num_blocks_skipped * 32 + the_chosen_block_id * 32;
    X += the_chosen_block_id * 32;
    Y += the_chosen_block_id * 32;

    magma_int_t blocks;

    if (m % zhemv_bs==0)
        blocks = m / zhemv_bs;
    else
        blocks = m / zhemv_bs + 1;

    blocks -= the_chosen_block_id;

    dim3 grid(blocks, 1, 1);
    dim3 threads(32, 8, 1);
    dim3 threads_u(32, 1, 1);

        
    /*
         * If matrix size is multiple of zhemv_bs, we use a specific code.
         * otherwise, we call the generic case.
         */
        if(m % zhemv_bs == 0 ) {
        magmablas_zhemv_200_U_special_mgpu_32_offset <<< grid, threads, 0, magma_stream >>>(
            m, alpha, A, lda, X, incx, beta, Y, incy, dC_work, my_gpu_id, num_gpus, nb, kstan);
        }
        else{
        magma_int_t m_mod_thread_x = m%zhemv_bs - 1;


        magmablas_zhemv_200_U_generic_mgpu_32_offset <<< grid, threads, 0, magma_stream >>> (
            m, alpha, A, lda, X, incx ,beta, Y, incy, dC_work, m_mod_thread_x, my_gpu_id, num_gpus, nb, kstan, the_right_gpu);
        }

        magmablas_zhemv_200_U_update_mgpu_32_offset<<< grid, threads_u, 0, magma_stream >>>(
        m, alpha, A, lda, X, incx, beta, Y, incy, dC_work, my_gpu_id, num_gpus, nb, kstan);
    
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
magmablas_zhemv_200_mgpu_32_offset( char uplo, magma_int_t n,
                      cuDoubleComplex alpha,
                      cuDoubleComplex **A, magma_int_t lda,
                      cuDoubleComplex **X, magma_int_t incx,
                      cuDoubleComplex beta,
                      cuDoubleComplex **Y, magma_int_t incy,
                      cuDoubleComplex **work, magma_int_t lwork,
              magma_int_t num_gpus, 
              magma_int_t nb,
                      magma_int_t offset,
                      cudaStream_t stream[][10])

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

        magma_int_t blocks    = n / zhemv_bs + (n % zhemv_bs != 0);
        magma_int_t workspace = lda * (blocks + 1);

        if (lwork < workspace){
           printf("Not enough work space in magmablas_zhemv: passed %d, required %d\n",
                  lwork, workspace);
           exit(1);
        }
        if(nb != 32)
        {
        printf("Error in magmablas_zsymv_200_mgpu: nb != 32, program will exit! please reallocate your matrix among GPUs\n");
        exit(0);
        }
        magma_int_t i = 0;
        for(i=0; i<num_gpus; i++)
        {
             magma_setdevice(i);
             magmablasSetKernelStream(stream[i][0]);

             magma_int_t the_chosen_block_id = offset / 32; 
         magma_int_t the_chosen_gpu_id = the_chosen_block_id % num_gpus; 

         magma_int_t  num_blocks_skipped = the_chosen_block_id / num_gpus;

         if(i < the_chosen_gpu_id)     
             {
         num_blocks_skipped += 1;
             }
              
             int new_gpu_id = ( i + num_gpus - the_chosen_gpu_id ) % num_gpus;
             

             magma_int_t the_right_block_id = n / nb ;
             magma_int_t the_right_gpu = the_right_block_id % num_gpus;

             the_right_gpu = ( the_right_gpu + num_gpus - the_chosen_gpu_id ) % num_gpus;
             // the_right_gpu is used in Upper generic case.

         if ( upper)
             { 
                 magmablas_zhemv_200_U_mgpu_32_offset(n, alpha, A[i], lda, X[i], incx, beta, Y[i], incy, work[i], 
                                                         new_gpu_id, num_gpus, nb, offset, num_blocks_skipped, the_right_gpu);     
             }
              else
             {
                 magmablas_zhemv_200_L_mgpu_32_offset(n, alpha, A[i], lda, X[i], incx, beta, Y[i], incy, work[i], 
                                                        new_gpu_id, num_gpus, nb, offset, num_blocks_skipped);
             }
        
      }



    return MAGMA_SUCCESS;
}

extern "C"
magma_int_t
magmablas_zhemv2_200_mgpu_32_offset( char uplo, magma_int_t n,
                      cuDoubleComplex alpha,
                      cuDoubleComplex **A, magma_int_t lda,
                      cuDoubleComplex **X, magma_int_t incx,
                      cuDoubleComplex beta,
                      cuDoubleComplex **Y, magma_int_t incy,
                      cuDoubleComplex **work, magma_int_t lwork,
              magma_int_t num_gpus, 
              magma_int_t nb,
                      magma_int_t offset)

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

        magma_int_t blocks    = n / thread_x + (n % thread_x != 0);
        magma_int_t workspace = lda * (blocks + 1);

        if (lwork < workspace){
           printf("Not enough work space in magmablas_zhemv: passed %d, required %d\n",
                  lwork, workspace);
           exit(1);
        }
        if(nb != 32)
        {
        printf("Error in magmablas_zsymv_200_mgpu: nb != 32, program will exit! please reallocate your matrix among GPUs\n");
        exit(0);
        }
        magma_int_t i = 0;
        for(i=0; i<num_gpus; i++)
        {
             magma_setdevice(i);
            // magmablasSetKernelStream(stream[i][0]);

             magma_int_t the_chosen_block_id = offset / 32; 
         magma_int_t the_chosen_gpu_id = the_chosen_block_id % num_gpus; 

         magma_int_t  num_blocks_skipped = the_chosen_block_id / num_gpus;

         if(i < the_chosen_gpu_id)     
             {
         num_blocks_skipped += 1;
             }
              
             int new_gpu_id = ( i + num_gpus - the_chosen_gpu_id ) % num_gpus;
             

             magma_int_t the_right_block_id = n / nb ;
             magma_int_t the_right_gpu = the_right_block_id % num_gpus;

             the_right_gpu = ( the_right_gpu + num_gpus - the_chosen_gpu_id ) % num_gpus;
             // the_right_gpu is used in Upper generic case.

         if ( upper)
             { 
                 magmablas_zhemv_200_U_mgpu_32_offset(n, alpha, A[i], lda, X[i], incx, beta, Y[i], incy, work[i], 
                                                         new_gpu_id, num_gpus, nb, offset, num_blocks_skipped, the_right_gpu);     
             }
              else
             {
                 magmablas_zhemv_200_L_mgpu_32_offset(n, alpha, A[i], lda, X[i], incx, beta, Y[i], incy, work[i], 
                                                        new_gpu_id, num_gpus, nb, offset, num_blocks_skipped);
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
