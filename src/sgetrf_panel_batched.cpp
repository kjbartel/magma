/*
    -- MAGMA (version 1.6.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2013
       
       @author Azzam Haidar
       @author Tingxing Dong

       @generated from zgetrf_panel_batched.cpp normal z -> s, Sat Nov 15 19:54:10 2014
*/
#include "common_magma.h"



////////////////////////////////////////////////////////////////////////////////////////
extern "C" magma_int_t
magma_sgetrf_recpanel_batched_q(
    magma_int_t m, magma_int_t n, magma_int_t min_recpnb,    
    float** dA_array,    magma_int_t ldda,
    magma_int_t** dipiv_array, magma_int_t** dpivinfo_array,
    float** dX_array,    magma_int_t dX_length,
    float** dinvA_array, magma_int_t dinvA_length,
    float** dW1_displ, float** dW2_displ,  
    float** dW3_displ, float** dW4_displ,
    float** dW5_displ,
    magma_int_t *info_array, magma_int_t gbstep,  
    magma_int_t batchCount, magma_queue_t stream, cublasHandle_t myhandle)
{

    //magma_int_t DEBUG = 3;
    // Quick return if possible
    if (m ==0 || n == 0) {
        return 0;
    }


    float **dA_displ  = NULL;
    magma_malloc((void**)&dA_displ,   batchCount * sizeof(*dA_displ));
    magma_int_t **dipiv_displ = NULL;
    magma_malloc((void**)&dipiv_displ, batchCount * sizeof(*dipiv_displ));
    
    magma_int_t panel_nb = n;
    if(panel_nb <= min_recpnb){
        //if(DEBUG>0)printf("calling bottom panel recursive with m=%d nb=%d\n",m,n);
        //  panel factorization
        //magma_sdisplace_pointers(dA_displ, dA_array, ldda, 0, 0, batchCount);
        magma_sgetf2_batched(
                           m, panel_nb,
                           dA_array, ldda,
                           dW1_displ, dW2_displ, dW3_displ,
                           dipiv_array, info_array, gbstep, batchCount, myhandle);
    }
    else{
        // split A over two [A A2]
        // panel on A1, update on A2 then panel on A1    
        magma_int_t n1 = n/2;
        magma_int_t n2 = n-n1;
        magma_int_t m1 = m;
        magma_int_t m2 = m-n1;
        magma_int_t p1 = 0;
        magma_int_t p2 = n1;
        // panel on A1
        //if(DEBUG>0)printf("calling recursive panel on A1 with m=%d nb=%d min_recpnb %d\n",m1,n1,min_recpnb);
        magma_sdisplace_pointers(dA_displ, dA_array, ldda, p1, p1, batchCount); 
        magma_idisplace_pointers(dipiv_displ, dipiv_array, 1, p1, 0, batchCount);
        magma_sgetrf_recpanel_batched_q(
                           m1, n1, min_recpnb,
                           dA_displ, ldda,
                           dipiv_displ, dpivinfo_array,
                           dX_array, dX_length,
                           dinvA_array, dinvA_length,
                           dW1_displ, dW2_displ,
                           dW3_displ, dW4_displ, dW5_displ,
                           info_array, gbstep, batchCount, stream, myhandle);

        // update A2
        //if(DEBUG>0)printf("calling TRSM  with             m=%d n=%d \n",m1,n2);
        
        // setup pivinfo 
        setup_pivinfo_batched_q(dpivinfo_array, dipiv_displ, m1, n1, stream, batchCount);
        magma_sdisplace_pointers(dW5_displ, dA_array, ldda, p1, p2, batchCount); 
        magma_slaswp_rowparallel_batched( n2, dW5_displ, ldda,
                           dX_array, n1,
                           0, n1,
                           dpivinfo_array, batchCount);
        magmablas_strsm_outofplace_batched(MagmaLeft, MagmaLower, MagmaNoTrans, MagmaUnit, 1,
                              n1, n2,
                              MAGMA_S_ONE,
                              dA_displ,    ldda, // dA
                              dX_array,  n1, // dB
                              dW5_displ,   ldda, // dX
                              dinvA_array, dinvA_length,
                              dW1_displ,   dW2_displ, 
                              dW3_displ,   dW4_displ,
                              0, batchCount);

        magma_sdisplace_pointers(dW1_displ, dA_array, ldda, p2, 0, batchCount); 
        magma_sdisplace_pointers(dA_displ, dA_array, ldda, p2, p2, batchCount); 

        //if(DEBUG>0)printf("calling update A2(%d,%d) -= A(%d,%d)*A(%d,%d)  with             m=%d n=%d k=%d ldda %d\n",p2,p2,p2,0,p1,p2,m2,n2,n1,ldda);

#if 0
        float neg_one = MAGMA_S_NEG_ONE;
        float one  = MAGMA_S_ONE;
        cublasSgemmBatched(myhandle, CUBLAS_OP_N, CUBLAS_OP_N, m2, n2, n1,
                                     &neg_one, (const float**) dW1_displ, ldda,
                                               (const float**) dW5_displ, ldda,
                                     &one,  dA_displ, ldda, batchCount );


#else

        magmablas_sgemm_batched( MagmaNoTrans, MagmaNoTrans, m2, n2, n1, 
                              MAGMA_S_NEG_ONE, dW1_displ, ldda, 
                              dW5_displ, ldda, 
                              MAGMA_S_ONE,  dA_displ, ldda, 
                              batchCount);
#endif

        // panel on A2
        //if(DEBUG>0)printf("calling recursive panel on A2 with m=%d nb=%d min_recpnb %d\n",m2,n2,min_recpnb);
        magma_idisplace_pointers(dipiv_displ, dipiv_array, 1, p2, 0, batchCount);
        magma_sgetrf_recpanel_batched_q(
                           m2, n2, min_recpnb,
                           dA_displ, ldda,
                           dipiv_displ, dpivinfo_array,
                           dX_array, dX_length,
                           dinvA_array, dinvA_length,
                           dW1_displ, dW2_displ,
                           dW3_displ, dW4_displ, dW5_displ,
                           info_array, gbstep+p2, batchCount, stream, myhandle);

        // setup pivinfo
        setup_pivinfo_batched_q(dpivinfo_array, dipiv_displ, m2, n2, stream, batchCount);
        adjust_ipiv_batched_q(dipiv_displ, n2, n1, magma_stream, batchCount);
        
        magma_sdisplace_pointers(dW1_displ, dA_array, ldda, p2, 0, batchCount); // no need since it is above
        magma_slaswp_rowparallel_batched( n1, dW1_displ, ldda,
                           dW1_displ, ldda,
                           n1, n,
                           dpivinfo_array, batchCount);

        
    }

    magma_free(dA_displ);
    magma_free(dipiv_displ);
    return 0;
}


////////////////////////////////////////////////////////////////////////////////////////

