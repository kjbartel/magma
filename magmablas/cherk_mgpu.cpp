/*
    -- MAGMA (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015
       
       @author Azzam Haidar
       @author Ichi Yamazaki

       @generated from zherk_mgpu.cpp normal z -> c, Fri Jan 30 19:00:10 2015

*/
#include "common_magma.h"
#include "trace.h"

/**
    Purpose
    -------
    This cherk_mgpu is internal routine used by cpotrf_mgpu_right.
    it has specific assumption on the block diagonal.
    
    @ingroup magma_cherk_comp
    ********************************************************************/

extern "C" void
magma_cherk_mgpu(
    magma_int_t ngpu,
    magma_uplo_t uplo, magma_trans_t trans, magma_int_t nb, magma_int_t n, magma_int_t k,
    float alpha,
    magmaFloatComplex_ptr dB[], magma_int_t lddb, magma_int_t b_offset,
    float beta,
    magmaFloatComplex_ptr dC[], magma_int_t lddc, magma_int_t c_offset,
    magma_int_t nqueue, magma_queue_t queues[][10])
{
#define dB(id, i, j)  (dB[(id)]+(j)*lddb + (i)+b_offset)
#define dC(id, i, j)  (dC[(id)]+(j)*lddc + (i))
#define STREAM_ID(i) (nqueue > 1 ? 1+((i)/nb)%(nqueue-1) : 0)

    magma_int_t i, id, ib, ii, kk, n1;
    magmaFloatComplex z_alpha = MAGMA_C_MAKE(alpha,0.0);
    magmaFloatComplex z_beta  = MAGMA_C_MAKE(beta, 0.0);

    magma_device_t orig_dev;
    magma_getdevice( &orig_dev );
    magma_queue_t orig_stream;
    magmablasGetKernelStream( &orig_stream );
    
    /* diagonal update */
    for( i=0; i < n; i += nb ) {
        id = ((i+c_offset)/nb)%ngpu;
        kk = STREAM_ID( i+c_offset );

        ib = min(nb, n-i);
        ii = nb*((i+c_offset)/(nb*ngpu));

        /* cher2k on diagonal block */
        magma_setdevice(id);
        magmablasSetKernelStream( queues[id][kk] );
        trace_gpu_start( id, kk, "syr2k", "syr2k" );
        magma_cherk(uplo, trans, ib, k,
                    alpha,  dB(id, i,          0 ), lddb,
                     beta,  dC(id, i+c_offset, ii), lddc);
        trace_gpu_end( id, kk );
    }

    /* off-diagonal update */
    if (uplo == MagmaUpper) {
        for( i=nb; i < n; i += nb ) {
            id = ((i+c_offset)/nb)%ngpu;
            kk = STREAM_ID( i+c_offset );

            ib = min(nb, n-i);
            ii = nb*((i+c_offset)/(nb*ngpu));

            magma_setdevice(id);
            magmablasSetKernelStream( queues[id][kk] );
            magma_cgemm(MagmaNoTrans, MagmaConjTrans, i, ib, k,
                        z_alpha, dB(id, 0, 0 ), lddb,
                                 dB(id, i, 0 ), lddb,
                        z_beta,  dC(id, 0, ii), lddc);
        }
    }
    else {
        for( i=0; i < n-nb; i += nb ) {
            id = ((i+c_offset)/nb)%ngpu;
            kk = STREAM_ID( i+c_offset );

            ib = min(nb, n-i);
            ii = nb*((i+c_offset)/(nb*ngpu));
            n1 = n-i-ib;

            /* cgemm on off-diagonal blocks */
            magma_setdevice(id);
            magmablasSetKernelStream( queues[id][kk] );
            trace_gpu_start( id, kk, "gemm_up", "gemm_up" );
            magma_cgemm(MagmaNoTrans, MagmaConjTrans, n1, ib, k,
                        z_alpha, dB(id, i+ib,           0 ), lddb,
                                 dB(id,  i,             0 ), lddb,
                        z_beta,  dC(id,  i+c_offset+ib, ii), lddc);
            trace_gpu_end( id, kk );
        }
    }

    // TODO why not sync?
    //for( id=0; id < ngpu; id++ ) {
    //    magma_setdevice(id);
    //    //for( kk=0; kk < nqueue; kk++ )
    //    //    magma_queue_sync( queues[id][kk] );
    //}
    magma_setdevice( orig_dev );
    magmablasSetKernelStream( orig_stream );
}
#undef dB
#undef dC
#undef STREAM_ID

// ----------------------------------------------------------------------
extern "C" void
magma_cherk_mgpu2(
    magma_int_t ngpu,
    magma_uplo_t uplo, magma_trans_t trans, magma_int_t nb, magma_int_t n, magma_int_t k,
    float alpha,
    magmaFloatComplex_ptr dB[], magma_int_t lddb, magma_int_t b_offset,
    float beta,
    magmaFloatComplex_ptr dC[], magma_int_t lddc, magma_int_t c_offset,
    magma_int_t nqueue, magma_queue_t queues[][10])
{
#define dB(id, i, j)  (dB[(id)]+(j)*lddb + (i)+b_offset)
#define dC(id, i, j)  (dC[(id)]+(j)*lddc + (i))
#define STREAM_ID(i) (nqueue > 1 ? 1+((i)/nb)%(nqueue-1) : 0)

    magma_int_t i, id, ib, ii, kk, n1;
    magmaFloatComplex z_alpha = MAGMA_C_MAKE(alpha,0.0);
    magmaFloatComplex z_beta  = MAGMA_C_MAKE(beta, 0.0);

    magma_device_t orig_dev;
    magma_getdevice( &orig_dev );
    magma_queue_t orig_stream;
    magmablasGetKernelStream( &orig_stream );
    
    /* diagonal update */
    for( i=0; i < n; i += nb ) {
        id = ((i+c_offset)/nb)%ngpu;
        kk = STREAM_ID( i+c_offset );

        ib = min(nb, n-i);
        ii = nb*((i+c_offset)/(nb*ngpu));
    }

    if (uplo == MagmaUpper) {
        for( i=0; i < n; i += nb ) {
            id = ((i+c_offset)/nb)%ngpu;
            kk = STREAM_ID( i+c_offset );

            ib = min(nb, n-i);
            ii = nb*((i+c_offset)/(nb*ngpu));
            n1 = i+ib;

            magma_setdevice(id);
            magmablasSetKernelStream( queues[id][kk] );

            /* cgemm on diag and off-diagonal blocks */
            magma_cgemm(MagmaNoTrans, MagmaConjTrans, n1, ib, k,
                        z_alpha, dB(id, 0, 0 ), lddb,
                                 dB(id, i, 0 ), lddb,
                        z_beta,  dC(id, 0, ii), lddc);
        }
    }
    else {
        for( i=0; i < n; i += nb ) {
            id = ((i+c_offset)/nb)%ngpu;
            kk = STREAM_ID( i+c_offset );

            ib = min(nb, n-i);
            ii = nb*((i+c_offset)/(nb*ngpu));
            n1 = n-i;

            magma_setdevice(id);
            magmablasSetKernelStream( queues[id][kk] );
            trace_gpu_start( id, kk, "gemm_up", "gemm_up" );
            /* cgemm on diag and off-diagonal blocks */
            magma_cgemm(MagmaNoTrans, MagmaConjTrans, n1, ib, k,
                        z_alpha, dB(id, i,           0), lddb,
                                 dB(id, i,           0), lddb,
                        z_beta,  dC(id, i+c_offset, ii), lddc);
            trace_gpu_end( id, kk );
        }
    }

    // TODO: why not sync?
    //for( id=0; id < ngpu; id++ ) {
    //    magma_setdevice(id);
    //    //for( kk=0; kk < nqueue; kk++ )
    //    //    magma_queue_sync( queues[id][kk] );
    //}
    magma_setdevice( orig_dev );
    magmablasSetKernelStream( orig_stream );
}

#undef dB
#undef dC
#undef STREAM_ID

