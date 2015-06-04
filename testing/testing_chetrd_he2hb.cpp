/*
    -- MAGMA (version 1.3.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2010

       @generated c Wed Nov 14 22:54:25 2012

*/

// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <cuda_runtime_api.h>
#include <cublas.h>

// includes, project
#include "flops.h"
#include "magma.h"
#include "magma_lapack.h"
#include "testings.h"

#if defined(USEMKL)
#include <mkl_service.h>
#endif
#if defined(USEACML)
#include <omp.h>
#endif

// Flops formula
#define PRECISION_c
#if defined(PRECISION_z) || defined(PRECISION_c)
#define FLOPS(n) ( 6. * FMULS_HETRD(n) + 2. * FADDS_HETRD(n))
#else
#define FLOPS(n) (      FMULS_HETRD(n) +      FADDS_HETRD(n))
#endif


extern "C" magma_int_t
magma_chetrd_he2hb( char uplo, magma_int_t n, magma_int_t NB,
                    cuFloatComplex *a, magma_int_t lda,
                    cuFloatComplex *tau,
                    cuFloatComplex *work, magma_int_t lwork,
                    cuFloatComplex *dT, magma_int_t threads,  
                    magma_int_t *info);

extern "C" magma_int_t
magma_chetrd_bhe2trc( int THREADS, int WANTZ, char uplo, int NE, int n, int NB, 
                   cuFloatComplex *A, int LDA, float *D, float *E, cuFloatComplex *dT1, int ldt1);

extern "C" magma_int_t
magma_chetrd_bhe2trc_v5(magma_int_t threads, magma_int_t wantz, char uplo, magma_int_t ne, magma_int_t n, magma_int_t nb,
                        cuFloatComplex *A, magma_int_t lda, float *D, float *E,
                        cuFloatComplex *dT1, magma_int_t ldt1);

#if defined(PRECISION_z) || defined(PRECISION_d)
extern "C" void cmp_vals(int n, float *wr1, float *wr2, float *nrmI, float *nrm1, float *nrm2);
extern "C" void ccheck_eig_(char *JOBZ, int  *MATYPE, int  *N, int  *NB,
                       cuFloatComplex* A, int  *LDA, float *AD, float *AE, float *D1, float *EIG,
                    cuFloatComplex *Z, int  *LDZ, cuFloatComplex *WORK, float *RWORK, float *RESU);
#endif

/* ////////////////////////////////////////////////////////////////////////////
   -- Testing chetrd_he2hb
*/
int main( int argc, char** argv)
{
    TESTING_CUDA_INIT();

    magma_timestr_t       start, end;
    float           eps, flops, gpu_perf, gpu_time;
    cuFloatComplex *h_A, *h_R, *h_work, *dT1;
    cuFloatComplex *tau;
    float *D, *E;

    /* Matrix size */
    magma_int_t N = 0, n2, lda, lwork,ldt;
    magma_int_t size[10] = {1024,2048,3072,4032,5184,6016,7040,8064,9088,10112};

    magma_int_t i, j, k, info, checkres, once = 0;
    magma_int_t ione     = 1;
    magma_int_t ISEED[4] = {0,0,0,1};
    const char *uplo = MagmaLowerStr;

    magma_int_t WANTZ=0;
    magma_int_t THREADS=1;
    magma_int_t NE = 0;
    magma_int_t NB = 0;
    checkres  = 0; //getenv("MAGMA_TESTINGS_CHECK") != NULL;

    if (argc != 1){
        for(i = 1; i<argc; i++){
            if (strcmp("-N", argv[i])==0) {
                N = atoi(argv[++i]);
                once = 1;
            }
            else if (strcmp("-NB", argv[i])==0) {
                NB = atoi(argv[++i]);
            }
            else if (strcmp("-threads", argv[i])==0) {
                THREADS = atoi(argv[++i]);
            }
            else if (strcmp("-wantz", argv[i])==0) {
                WANTZ = atoi(argv[++i]);
            }
            else if (strcmp("-NE", argv[i])==0) {
                NE = atoi(argv[++i]);
            }
            else if ( strcmp("-c", argv[i]) == 0 ) {
                checkres = 1;
            }
            else if (strcmp("-U", argv[i])==0)
                uplo = MagmaUpperStr;
            else if (strcmp("-L", argv[i])==0)
                uplo = MagmaLowerStr;
        }
        if ( N > 0 )
            printf("  testing_chetrd_he2hb -L|U -N %d -NB %d   -wantz %d   -threads %d    check %d \n\n", N, NB, WANTZ, THREADS, checkres);
        else
        {
            printf("\nUsage: \n");
            printf("  testing_chetrd_he2hb -L|U -N %d -NB  -wantz -threads \n\n", 1024);
            exit(1);
        }
    }
    else {
        printf("\nUsage: \n");
        printf("  testing_chetrd_he2hb -L|U -N %d\n\n", 1024);
        N = size[9];
    }
        
 printf ("HELLOOOOOOOO\n");

    eps = lapackf77_slamch( "E" );
    lda = N;
    ldt = N;
    n2  = lda * N; 
    if(NB<1)
        NB  = 64; //64; //magma_get_chetrd_he2hb_nb(N);

    if(NE<1)
        NE  = N; //64; //magma_get_chetrd_he2hb_nb(N);

    /* We suppose the magma NB is bigger than lapack NB */
    lwork = N*NB; 

    /* Allocate host memory for the matrix */
    TESTING_MALLOC(    h_A,    cuFloatComplex, lda*N );
    TESTING_HOSTALLOC( h_R,    cuFloatComplex, lda*N );
    TESTING_HOSTALLOC( h_work, cuFloatComplex, lwork );
    TESTING_MALLOC(    tau,    cuFloatComplex, N-1   );
    TESTING_HOSTALLOC( D,    float, N );
    TESTING_HOSTALLOC( E,    float, N );
    //TESTING_DEVALLOC( dT1,  cuFloatComplex, (2*min(N,N)+(N+31)/32*32)*NB );
    TESTING_DEVALLOC( dT1,  cuFloatComplex, (N*NB) );

    printf("  N    GPU GFlop/s   \n");
    printf("=====================\n");
    for(i=0; i<10; i++){
        if ( !once ) {
            N = size[i];
        }
        lda  = N;
        n2   = N*lda;
        flops = FLOPS( (float)N ) / 1e6;
        if(WANTZ) flops = 2.0*flops;

        /* ====================================================================
           Initialize the matrix
           =================================================================== */
        lapackf77_clarnv( &ione, ISEED, &n2, h_A );

        // Make the matrix hermitian 
        {
            magma_int_t i, j;
            for(i=0; i<N; i++) {
                MAGMA_C_SET2REAL( h_A[i*lda+i], ( MAGMA_C_REAL(h_A[i*lda+i]) ) );
                for(j=0; j<i; j++)
                    h_A[i*lda+j] = cuConjf(h_A[j*lda+i]);
            }
        }
/*
            for(i=0; i<N; i++){ 
                for(j=0; j<N; j++){
                MAGMA_C_SET2REAL( h_A[i*lda+j], ( MAGMA_C_REAL(h_A[i*lda+j]) ) );
                }
            }
*/

/*
    FILE *trace_file;
    trace_file = fopen("AJETE/Ainit", "w");
    for (j = 0; j < N ; j++) 
          for (i = 0; i < N ; i++) 
                         fprintf(trace_file,"%10d %10d %25.15e %25.15e\n",i+1,j+1,MAGMA_C_REAL(h_A[j*lda+i]) ,  MAGMA_C_IMAG(h_A[j*lda+i])  );
    fclose(trace_file);
*/



        lapackf77_clacpy( MagmaUpperLowerStr, &N, &N, h_A, &lda, h_R, &lda );

/*
lapackf77_clarnv( &ione, ISEED, &N, D );
lapackf77_clarnv( &ione, ISEED, &N, E );
i= min(12,THREADS);
mkl_set_num_threads( i );
start = get_current_time();
dstedc_withZ('V', N, D, E, h_R, lda);
end = get_current_time();
printf("  Finish EIGEN   timing= %f  threads %d ---> 0000000 \n" ,GetTimerValue(start,end) / 1000., i);
mkl_set_num_threads( 1 );
return 0;
*/
/*

    FILE *trace_file;
    trace_file = fopen("AJETE/Ainit", "w");
    for (j = 0; j < N ; j++) 
          for (i = 0; i < N ; i++) 
                         fprintf(trace_file,"%10d %10d %25.15e %25.15e\n",i+1,j+1,MAGMA_C_REAL(h_R[j*lda+i]) ,  MAGMA_C_IMAG(h_R[j*lda+i])  );
    fclose(trace_file);
*/
    /*
    int pm,pn,indi,indj,n=N;
    i=1;
                      indi = i+NB;
                  indj = i;
                  pm   = n - i - NB + 1;
                  pn   = min(i+NB-1, n-NB) -i + 1;
*/
                  /*
                  printf("voici pm pn %d %d \n",pm,pn);
              lapackf77_cgeqrf(&pm, &pn, &h_R[NB], &lda, 
                             tau, h_work, &lwork, &info);
              printf("TOTOTOTO INFO %d\n",info);
              memset(h_work, 0, lwork*sizeof(cuFloatComplex));
            lapackf77_clarft( "F", "C",
                              &pm, &pn, &h_R[NB], &lda,
                              tau, h_work, &pn);

*/


        /* ====================================================================
           Performs operation using MAGMA
           =================================================================== */
        start = get_current_time();
        //magma_chetrd_he2hb(uplo[0], N, h_R, lda, tau, h_work, lwork, &info);
        magma_chetrd_he2hb(uplo[0], N, NB, h_R, lda, tau, h_work, lwork, dT1, THREADS, &info);
        end = get_current_time();
        printf("  Finish BAND    timing= %f \n" ,GetTimerValue(start,end) / 1000.);
// return 0;
/*        
    int   Vblksiz=-1, blkcnt=-1, LDV=-1, LDT =-1, INgrsiz=1, LDE=-1, BAND=6;
    Vblksiz = NB; //min(NB,64);
    LDT     = Vblksiz;
    findVTsiz(N, NB, Vblksiz, &blkcnt, &LDV);
    cuFloatComplex *dVV2, *dTT2, dV3;
    int dVVs;
           dVVs = max(N*N,blkcnt*LDV*Vblksiz);
           printf("dvsize %f \n",(16.0*(real_Double_t)dVVs)*1e-9);
           if( CUBLAS_STATUS_SUCCESS != cublasAlloc(dVVs, sizeof(cuFloatComplex), (void**)&dVV2) ) { 
               printf ("!!!! -------> cublasAlloc failed for: dVV2\n" );       
               exit(-1);                                                           
           }
    
           if( CUBLAS_STATUS_SUCCESS != cublasAlloc( dVVs, sizeof(cuFloatComplex), (void**)&dTT2) ) { 
              printf ("!!!! ---------> cublasAlloc failed for: dTT2\n" );       
              exit(-1);                                                           
           }
    
           if( CUBLAS_STATUS_SUCCESS != cublasAlloc( dVVs, sizeof(cuFloatComplex), (void**)&dV3) ) { 
              printf ("!!!! ---------> cublasAlloc failed for: dV3\n" );       
              exit(-1);                                                           
           }

           printf("done from alloc exit\n");
           */

            /*        
    trace_file = fopen("AJETE/Aafter", "w");
    for (j = 0; j < N ; j++) 
          for (i = 0; i < N ; i++) 
                         fprintf(trace_file,"%10d%10d%40.30e\n",i+1,j+1,h_R[j*lda+i]);
    fclose(trace_file);
*/
 /*
        memset(h_work, 0, lwork*sizeof(cuFloatComplex));
        magma_cgetmatrix( pn, pn, dT1, ldt, h_work, pn );
   trace_file = fopen("AJETE/T", "w");
    for (j = 0; j < pn ; j++) 
          for (i = 0; i < pn ; i++) 
                         fprintf(trace_file,"%10d%10d%40.30e\n",i+1,j+1,h_work[j*pn+i]);
    fclose(trace_file);
*/        
        
        //        dsytrd_bsy2trc(THREADS, uplo[0], N, NB, h_R, lda, D, E);
        magma_chetrd_bhe2trc_v5(THREADS, WANTZ, uplo[0], NE, N, NB, h_R, lda, D, E, dT1, ldt);
        end = get_current_time();
        if ( info < 0 )
            printf("Argument %d of magma_chetrd_he2hb had an illegal value\n", (int) -info);

        gpu_perf = flops / GetTimerValue(start,end);
        gpu_time = GetTimerValue(start,end) / 1000.;

        /* =====================================================================
           Check the factorization
           =================================================================== */
        /*
        if ( checkres ) {
            FILE        *fp ;

            printf("Writing input matrix in matlab_i_mat.txt ...\n");
            fp = fopen ("matlab_i_mat.txt", "w") ;
            if( fp == NULL ){ printf("Couldn't open output file\n"); exit(1);}

            for(j=0; j<N; j++)
                for(k=0; k<N; k++)
                    {
                        #if defined(PRECISION_z) || defined(PRECISION_c)
                        fprintf(fp, "%5d %5d %11.8f %11.8f\n", k+1, j+1, 
                                h_A[k+j*lda].x, h_A[k+j*lda].y);
                        #else
                        fprintf(fp, "%5d %5d %11.8f\n", k+1, j+1, h_A[k+j*lda]);
                        #endif
                    }
            fclose( fp ) ;

          printf("Writing output matrix in matlab_o_mat.txt ...\n");
          fp = fopen ("matlab_o_mat.txt", "w") ;
          if( fp == NULL ){ printf("Couldn't open output file\n"); exit(1);}

          for(j=0; j<N; j++)
            for(k=0; k<N; k++)
              {
                #if defined(PRECISION_z) || defined(PRECISION_c)
                fprintf(fp, "%5d %5d %11.8f %11.8f\n", k+1, j+1,
                        h_R[k+j*lda].x, h_R[k+j*lda].y);
                #else
                fprintf(fp, "%5d %5d %11.8f\n", k+1, j+1, h_R[k+j*lda]);
                #endif
              } 
          fclose( fp ) ;

        }*/



        /* =====================================================================
           Print performance and error.
           =================================================================== */
#if defined(CHECKEIG)
#if defined(PRECISION_z)  || defined(PRECISION_d)
        if ( checkres ) {
            printf("  Total N %5d  flops %6.2f  timing %6.2f seconds\n", (int) N, gpu_perf, gpu_time );
            char JOBZ;
            if(WANTZ==0) 
                    JOBZ='N';
            else
                    JOBZ = 'V';
            float nrmI=0.0, nrm1=0.0, nrm2=0.0;
            int    lwork2 = 256*N;
            cuFloatComplex *work2     = (cuFloatComplex *) malloc (lwork2*sizeof(cuFloatComplex));
            float *rwork2     = (float *) malloc (N*sizeof(float));
            float *D2          = (float *) malloc (N*sizeof(float));
            cuFloatComplex *AINIT    = (cuFloatComplex *) malloc (N*lda*sizeof(cuFloatComplex));
            memcpy(AINIT, h_A, N*lda*sizeof(cuFloatComplex));
            /* compute the eigenvalues using lapack routine to be able to compare to it and used as ref */
            start = get_current_time();
            i= min(12,THREADS);

#if defined(USEMKL)
            mkl_set_num_threads( i );
#endif
#if defined(USEACML)
            omp_set_num_threads(i);
#endif

#if defined(PRECISION_z) || defined (PRECISION_c)
            lapackf77_cheev( "N", "L", &N, h_A, &lda, D2, work2, &lwork2, rwork2, &info );
#else
            lapackf77_ssyev( "N", "L", &N, h_A, &lda, D2, work2, &lwork2, &info );
#endif
            ///* call eigensolver for our resulting tridiag [D E] and for Q */
            //dstedc_withZ('V', N, D, E, h_R, lda);
            ////ssterf_( &N, D, E, &info); 
            ////
            end = get_current_time();
            printf("  Finish CHECK - EIGEN   timing= %f  threads %d \n" ,GetTimerValue(start,end) / 1000., i);

            /*
            for(i=0;i<10;i++)
                printf(" voici lpk D[%d] %e\n",i,D2[i]);
            */

            //cuFloatComplex mydz=0.0,mydo=1.0;
            //cuFloatComplex *Z = (cuFloatComplex *) malloc(N*lda*sizeof(cuFloatComplex));
           // dgemm_("N","N",&N,&N,&N,&mydo,h_R,&lda,h_A,&lda,&mydz,Z,&lda);


           /* compare result */
           cmp_vals(N, D2, D, &nrmI, &nrm1, &nrm2);


           cuFloatComplex *WORKAJETER;
           float *RWORKAJETER, *RESU;
           WORKAJETER  = (cuFloatComplex *) malloc( (2* N * N + N) * sizeof(cuFloatComplex) );
           RWORKAJETER = (float *) malloc( N * sizeof(float) );
           RESU        = (float *) malloc(10*sizeof(float));
           int MATYPE;
           memset(RESU,0,10*sizeof(float));

 
           MATYPE=3;
           float NOTHING=0.0;
           start = get_current_time();
           // check results
           ccheck_eig_(&JOBZ, &MATYPE, &N, &NB, AINIT, &lda, &NOTHING, &NOTHING, D2 , D, h_R, &lda, WORKAJETER, RWORKAJETER, RESU );
           end = get_current_time();
           printf("  Finish CHECK - results timing= %f \n" ,GetTimerValue(start,end) / 1000.);
#if defined(USEMKL)
           mkl_set_num_threads( 1 );
#endif
#if defined(USEACML)
           omp_set_num_threads(1);
#endif

           printf("\n");
           printf(" ================================================================================================================\n");
           printf("   ==> INFO voici  threads=%d    N=%d    NB=%d   WANTZ=%d\n", (int) THREADS, (int) N, (int) NB, (int) WANTZ);
           printf(" ================================================================================================================\n");
           printf("            DSBTRD                : %15s \n", "STATblgv9withQ    ");
           printf(" ================================================================================================================\n");
           if(WANTZ>0)
              printf(" | A - U S U' | / ( |A| n ulp )   : %15.3E   \n",RESU[0]); 
           if(WANTZ>0)
              printf(" | I - U U' | / ( n ulp )         : %15.3E   \n", RESU[1]);
           printf(" | D1 - EVEIGS | / (|D| ulp)      : %15.3E   \n",  RESU[2]);
           printf(" max | D1 - EVEIGS |              : %15.3E   \n",  RESU[6]);
           printf(" ================================================================================================================\n\n\n");
       
           printf(" ****************************************************************************************************************\n");
           printf(" * Hello here are the norm  Infinite (max)=%e  norm one (sum)=%e   norm2(sqrt)=%e *\n",nrmI, nrm1, nrm2);
           printf(" ****************************************************************************************************************\n\n");

        } 
#endif         
#endif  

      printf("  Total N %5d  flops %6.2f        timing %6.2f seconds\n", (int) N, 0.0, gpu_time );
      printf("============================================================================\n\n\n");

      if ( once )
          break;
    }

    /* Memory clean up */
    TESTING_FREE( h_A );
    TESTING_FREE( tau ); 
    TESTING_HOSTFREE( h_R ); 
    TESTING_HOSTFREE( h_work ); 

    /* Shutdown */
    TESTING_CUDA_FINALIZE();
    return EXIT_SUCCESS;
}
