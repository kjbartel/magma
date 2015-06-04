/*
    -- MAGMA (version 1.5.0-beta2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2014

       @generated from testing_zhemm_mgpu.cpp normal z -> s, Fri May 30 10:41:21 2014
       
       @author Mark Gates
       @author Azzam Haidar
*/

// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <cuda_runtime_api.h>
#include <cublas.h>
#include <assert.h>

// includes, project
#include "flops.h"
#include "magma.h"
#include "magma_lapack.h"
#include "testings.h"

//#include "trace.h"

/* ////////////////////////////////////////////////////////////////////////////
   -- Testing magma_ssymm_mgpu
*/
int main( int argc, char** argv)
{
    TESTING_INIT();

    float c_neg_one = MAGMA_S_NEG_ONE;
    float calpha    = MAGMA_S_MAKE( 3.456, 5.678 );
    float cbeta     = MAGMA_S_MAKE( 1.234, 2.456 );
    
    real_Double_t    gflops, gpu_perf=0., cpu_perf=0., gpu_time=0., cpu_time=0.;
    real_Double_t    gpu_perf2=0., gpu_time2=0.;
    float           error=0., errorbis=0., work[1];
    float *hA, *hX, *hB, *hR;
    float *dA[MagmaMaxGPUs], *dX[MagmaMaxGPUs], *dB[MagmaMaxGPUs], *dwork[MagmaMaxGPUs], *hwork[MagmaMaxGPUs+1];
    float *dA2;
    magma_int_t M, N, size, lda, ldda, msize;
    magma_int_t ione     = 1;
    magma_int_t iseed[4] = {0,0,0,1};
    magma_int_t status = 0;
    
    magma_opts opts;
    parse_opts( argc, argv, &opts );
    
    float tol = opts.tolerance * lapackf77_slamch("E");
    
    magma_int_t gnode[MagmaMaxGPUs][MagmaMaxGPUs+2];
    magma_int_t nbcmplx = 0;
    magma_buildconnection_mgpu(gnode, &nbcmplx, opts.ngpu);
    printf("Initializing communication pattern... GPU-ncmplx %d\n\n", (int) nbcmplx);

    for (int i=0; i < nbcmplx; ++i) {
        int myngpu = gnode[i][MagmaMaxGPUs];
        printf("cmplx %d has %d gpu ", i, myngpu);
        for(int j=0; j < myngpu; ++j)
            printf("  %d", (int) gnode[i][j]);
        printf("\n");
    }

    magma_int_t nbevents = 2;
    magma_queue_t streams[MagmaMaxGPUs][20];
    magma_event_t redevents[MagmaMaxGPUs][20];
    magma_event_t redevents2[MagmaMaxGPUs][MagmaMaxGPUs*MagmaMaxGPUs+10];
    for( int d = 0; d < opts.ngpu; ++d ) {
        for( magma_int_t i = 0; i < opts.nstream; ++i ) {
            magma_queue_create( &streams[d][i] );
        }
        for( magma_int_t i = 0; i < nbevents; ++i ) {
            cudaEventCreateWithFlags(&redevents[d][i],  cudaEventDisableTiming);
            cudaEventCreateWithFlags(&redevents2[d][i], cudaEventDisableTiming);
        }
    }

    printf( "nb %d, ngpu %d, nstream %d version %d\n", (int) opts.nb, (int) opts.ngpu, (int) opts.nstream, (int) opts.version );
    printf("    M     N    nb offset  CPU GFlop/s (sec)   GPU GFlop/s (sec)   CUBLAS hemm (sec)   ||R|| / ||A||*||X||\n");
    printf("=========================================================================================================\n");
    for( int itest = 0; itest < opts.ntest; ++itest ) {
      for( int offset = 0; offset < 1; offset += min(N,opts.nb) ) {
        for( int iter = 0; iter < opts.niter; ++iter ) {
            M     = opts.msize[itest];
            N     = opts.nsize[itest];
            msize = M - offset;
            lda   = M;
            ldda  = ((M + 31)/32)*32;
            size  = lda*M;
            gflops = FLOPS_SSYMM( MagmaLeft, (float)msize, (float)N ) / 1e9;
            
            TESTING_MALLOC_CPU( hA, float, lda*M );
            TESTING_MALLOC_CPU( hX, float, lda*N );
            TESTING_MALLOC_CPU( hB, float, lda*N );
            
            TESTING_MALLOC_PIN( hR, float, lda*N );

            for( int d = 0; d < opts.ngpu; ++d ) {
                magma_int_t mlocal = ((M / opts.nb) / opts.ngpu + 1) * opts.nb;
                magma_setdevice( d );
                TESTING_MALLOC_DEV( dA[d],    float, ldda*mlocal );
                TESTING_MALLOC_DEV( dX[d],    float, ldda*N      );
                TESTING_MALLOC_DEV( dB[d],    float, ldda*N      );
                TESTING_MALLOC_DEV( dwork[d], float, ldda*N*3    );
                
                TESTING_MALLOC_PIN( hwork[d], float, lda*N       );
            }
            TESTING_MALLOC_PIN( hwork[opts.ngpu], float, lda*N );
        
            if ( opts.check ) {
                magma_setdevice( 0 );
                TESTING_MALLOC_DEV( dA2, float, ldda*M );
            }

            lapackf77_slarnv( &ione, iseed, &size, hA );
            magma_smake_symmetric( M, hA, lda );
            
            size = lda*N;
            lapackf77_slarnv( &ione, iseed, &size, hX );
            lapackf77_slarnv( &ione, iseed, &size, hB );
            lapackf77_slacpy( "Full", &M, &N, hB, &lda, hR, &lda );
            
            /* ====================================================================
               Performs operation using MAGMA
               =================================================================== */
            magma_ssetmatrix_1D_col_bcyclic( M, M, hA, lda, dA, ldda, opts.ngpu, opts.nb );
            for( int d = 0; d < opts.ngpu; ++d ) {
                magma_setdevice( d );
                //magmablasSetKernelStream( streams[ d ][  0 ] );
                magma_ssetmatrix( M, N, hX, lda, dX[d], ldda );
                //if (d == 0) magma_ssetmatrix( M, N, hB, lda, dB[d], ldda ); // this is wrong coz when offset != 0 the gpu who do the beta*C may be not 0 so this should be related to stdev(starting device who own i=0 first col)
                magma_ssetmatrix( M, N, hB, lda, dB[d], ldda );
            }
        
            //memset(hR, 0, lda*N*sizeof(float));
    
            //trace_init( 1, opts.ngpu, opts.nstream, (magma_queue_t*) streams );
    
            //magma_int_t offset = 0; //opts.nb;
    
            gpu_time = magma_sync_wtime(0);
            
            // 1GPU version light
            /*
            magmablas_ssymm_1gpu(
                MagmaLeft, MagmaLower, msize, N,
                calpha,    dA, ldda, offset,
                           dX, ldda,
                cbeta,     dB, ldda, hR, lda,
                opts.ngpu, opts.nb, streams, opts.nstream );
            */
            // multi gpu version
            
            if (opts.version == 21) {
                // TODO: not available?
                //magmablas_ssymm_mgpu(
                //    MagmaLeft, MagmaLower, msize, N,
                //    calpha,    dA, ldda, offset,
                //               dX, ldda,
                //    cbeta,     dB, ldda, dwork, ldda, hR, lda, hwork, lda,
                //    opts.ngpu, opts.nb, streams, opts.nstream, redevents, nbevents );
            }
            else {
                magmablas_ssymm_mgpu_com(
                    MagmaLeft, MagmaLower, msize, N,
                    calpha,    dA, ldda, offset,
                               dX, ldda,
                    cbeta,     dB, ldda, dwork, ldda, hR, lda, hwork, lda,
                    opts.ngpu, opts.nb, streams, opts.nstream, redevents2, nbevents, gnode, nbcmplx);
            }
           
            gpu_time = magma_sync_wtime(0) - gpu_time;
            gpu_perf = gflops / gpu_time;
                
            #ifdef TRACING
            char buf[80];
            snprintf( buf, sizeof(buf), "ssymm-m%d-n%d-nb%d-stream%d-ngpu%d-run%d.svg",
                      (int) M, (int) N, (int) opts.nb, (int) opts.nstream, (int) opts.ngpu, (int) j );
            trace_finalize( buf, "trace.css" );
            #endif
            
            /* ====================================================================
               Performs operation using CUBLAS
               =================================================================== */
            if ( opts.check && iter == 0 ) {
                magma_setdevice( 0 );
                magmablasSetKernelStream(  0  );
                magma_ssetmatrix( M, M, hA, lda, dA2, ldda );
                magma_ssetmatrix( M, N, hX, lda, dX[0], ldda );
                magma_ssetmatrix( M, N, hB, lda, dwork[0], ldda );
                
                gpu_time2 = magma_sync_wtime(0);
                magma_ssymm(
                    MagmaLeft, MagmaLower, msize, N,
                    calpha,    dA2+offset*ldda+offset, ldda,
                               dX[0],    ldda,
                    cbeta,     dwork[0], ldda );
                gpu_time2 = magma_sync_wtime(0) - gpu_time2;
                gpu_perf2 = gflops / gpu_time2;
            }
            
            /* =====================================================================
               Performs operation using LAPACK
               =================================================================== */
            if ( opts.check ) {
                // store ||A||*||X||
                errorbis  = lapackf77_slange("fro", &msize, &msize, hA+offset*lda+offset, &lda, work );
                errorbis *= lapackf77_slange("fro", &msize, &N, hX, &lda, work );
                
                //printf( "A =" ); magma_sprint( M, M, hA, lda );
                //printf( "X =" ); magma_sprint( M, N, hX, lda );
                //printf( "B =" ); magma_sprint( M, N, hB, lda );
                
                cpu_time = magma_wtime();
                blasf77_ssymm( "Left", "Lower", &msize, &N,
                                &calpha, hA+offset*lda+offset, &lda,
                                         hX, &lda,
                                &cbeta,  hB, &lda );
                cpu_time = magma_wtime() - cpu_time;
                cpu_perf = gflops / cpu_time;
                /*
                trace_file = fopen("AJETE/C", "w");
                for (int j = 0; j < N; j++)
                    for (int i = 0; i < siz; i++)
                        fprintf(trace_file, "%10d%10d%40.30e\n", i+1, j+1, hB[j*lda+i]);
                fclose(trace_file);
                */
                magma_int_t firstprint=0;
                for(magma_int_t dev=0; dev < opts.ngpu; ++dev) {
                    magma_setdevice( dev );
                    magma_sgetmatrix( M, N, dB[dev], ldda, hR, lda );
    
                    // compute relative error ||R||/||A||*||X||, where R := B_magma - B_lapack = R - B
                    size = lda*N;
                    blasf77_saxpy( &size, &c_neg_one, hB, &ione, hR, &ione );
                    error = lapackf77_slange("fro", &msize, &N, hR, &lda, work) / errorbis;
                    
                    //printf( "R ="  ); magma_sprint( M, N, hR, lda );
                    if (firstprint == 0) {
                        printf( "%5d %5d %5d %5d   %7.1f (%7.4f)   %7.1f (%7.4f)   %7.1f (%7.4f)   %8.2e   %s\n",
                                (int) M, (int) N, (int) opts.nb, (int) offset,
                                cpu_perf, cpu_time,
                                gpu_perf, gpu_time,
                                gpu_perf2, gpu_time2,
                                error, (error < tol ? "ok" : "failed") );
                    }
                    else {
                        printf( "%89s  %8.2e   %s\n", " ",
                                error, (error < tol ? "ok" : "failed") );
                    }
                    status += ! (error < tol);
                    firstprint =1;
                }
            } else {
                printf( "%5d %5d %5d %5d     ---   (  ---  )   %7.1f (%7.4f)     ---   (  ---  )   ---\n",
                        (int) M, (int) N, (int) opts.nb, (int) offset,
                        gpu_perf, gpu_time );
            }
    
            TESTING_FREE_CPU( hA );
            TESTING_FREE_CPU( hX );
            TESTING_FREE_CPU( hB );
            
            TESTING_FREE_PIN( hR );
        
            for( int d = 0; d < opts.ngpu; ++d ) {
                magma_setdevice( d );
                TESTING_FREE_DEV( dA[d]    );
                TESTING_FREE_DEV( dX[d]    );
                TESTING_FREE_DEV( dB[d]    );
                TESTING_FREE_DEV( dwork[d] );
                
                TESTING_FREE_PIN( hwork[d] );
            }
            TESTING_FREE_PIN( hwork[opts.ngpu] );
        
            if ( opts.check ) {
                magma_setdevice( 0 );
                TESTING_FREE_DEV( dA2 );
            }
            fflush( stdout );
        }
        if ( opts.niter > 1 ) {
            printf( "\n" );
        }
      }  // offset
    }

    for( int d = 0; d < opts.ngpu; ++d ) {
        magma_setdevice( d );
        for( magma_int_t i = 0; i < opts.nstream; ++i ) {
            magma_queue_destroy( streams[d][i] );
        }
        for( magma_int_t i = 0; i < nbevents; ++i ) {
            magma_event_destroy( redevents[d][i]  );
            magma_event_destroy( redevents2[d][i] );
        }
    }
    
    TESTING_FINALIZE();
    return status;
}
