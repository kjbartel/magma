/*
    -- MAGMA (version 1.2.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       May 2012

       @author Mark Gates
       @generated c Tue May 15 18:17:08 2012
*/
#include "common_magma.h"

// -------------------------
// Put 0s in the upper triangular part of a panel and 1s on the diagonal.
// Stores previous values in work array, to be restored later with cq_to_panel.
extern "C"
void cpanel_to_q(char uplo, int ib, cuFloatComplex *A, int lda, cuFloatComplex *work)
{
    int i, j, k = 0;
    cuFloatComplex *col;
    cuFloatComplex c_zero = MAGMA_C_ZERO;
    cuFloatComplex c_one  = MAGMA_C_ONE;
    
    if (uplo == 'U' || uplo == 'u'){
        for(i = 0; i < ib; ++i){
            col = A + i*lda;
            for(j = 0; j < i; ++j){
                work[k] = col[j];
                col [j] = c_zero;
                ++k;
            }
            
            work[k] = col[i];
            col [j] = c_one;
            ++k;
        }
    }
    else {
        for(i=0; i<ib; ++i){
            col = A + i*lda;
            work[k] = col[i];
            col [i] = c_one;
            ++k;
            for(j=i+1; j<ib; ++j){
                work[k] = col[j];
                col [j] = c_zero;
                ++k;
            }
        }
    }
}


// -------------------------
// Restores a panel, after call to cpanel_to_q.
extern "C"
void cq_to_panel(char uplo, int ib, cuFloatComplex *A, int lda, cuFloatComplex *work)
{
    int i, j, k = 0;
    cuFloatComplex *col;
    
    if (uplo == 'U' || uplo == 'u'){
        for(i = 0; i < ib; ++i){
            col = A + i*lda;
            for(j = 0; j <= i; ++j){
                col[j] = work[k];
                ++k;
            }
        }
    }
    else {
        for(i = 0; i < ib; ++i){
            col = A + i*lda;
            for(j = i; j < ib; ++j){
                col[j] = work[k];
                ++k;
            }
        }
    }
}
