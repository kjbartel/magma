/*
    -- MAGMA (version 1.3.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2010

       @generated d Wed Nov 14 22:54:25 2012

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
#define PRECISION_d
#if defined(PRECISION_z) || defined(PRECISION_c)
#define FLOPS(n) ( 6. * FMULS_SYTRD(n) + 2. * FADDS_SYTRD(n))
#else
#define FLOPS(n) (      FMULS_SYTRD(n) +      FADDS_SYTRD(n))
#endif


extern "C" magma_int_t
magma_dsytrd_sy2sb( char uplo, magma_int_t n, magma_int_t NB,
                    double *a, magma_int_t lda,
                    double *tau,
                    double *work, magma_int_t lwork,
                    double *dT, magma_int_t threads,  
                    magma_int_t *info);

extern "C" magma_int_t
magma_dsytrd_bsy2trc( int THREADS, int WANTZ, char uplo, int NE, int n, int NB, 
                   double *A, int LDA, double *D, double *E, double *dT1, int ldt1);

extern "C" magma_int_t
magma_dsytrd_bsy2trc_v5(magma_int_t threads, magma_int_t wantz, char uplo, magma_int_t ne, magma_int_t n, magma_int_t nb,
                        double *A, magma_int_t lda, double *D, double *E,
                        double *dT1, magma_int_t ldt1);

#if defined(PRECISION_z) || defined(PRECISION_d)
extern "C" void cmp_vals(int n, double *wr1, double *wr2, double *nrmI, double *nrm1, double *nrm2);
extern "C" void dcheck_eig_(char *JOBZ, int  *MATYPE, int  *N, int  *NB,
                       double* A, int  *LDA, double *AD, double *AE, double *D1, double *EIG,
                    double *Z, int  *LDZ, double *WORK, double *RWORK, double *RESU);
#endif

/* ////////////////////////////////////////////////////////////////////////////
   -- Testing dsytrd_sy2sb
*/
int main( int argc, char** argv)
{
    TESTING_CUDA_INIT();

    magma_timestr_t       start, end;
    double           eps, flops, gpu_perf, gpu_time;
    double *h_A, *h_R, *h_work, *dT1;
    double *tau;
    double *D, *E;

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
            printf("  testing_dsytrd_sy2sb -L|U -N %d -NB %d   -wantz %d   -threads %d    check %d \n\n", N, NB, WANTZ, THREADS, checkres);
        else
        {
            printf("\nUsage: \n");
            printf("  testing_dsytrd_sy2sb -L|U -N %d -NB  -wantz -threads \n\n", 1024);
            exit(1);
        }
    }
    else {
        printf("\nUsage: \n");
        printf("  testing_dsytrd_sy2sb -L|U -N %d\n\n", 1024);
        N = size[9];
    }
        
 printf ("HELLOOOOOOOO\n");

    eps = lapackf77_dlamch( "E" );
    lda = N;
    ldt = N;
    n2  = lda * N; 
    if(NB<1)
        NB  = 64; //64; //magma_get_dsytrd_sy2sb_nb(N);

    if(NE<1)
        NE  = N; //64; //magma_get_dsytrd_sy2sb_nb(N);

    /* We suppose the magma NB is bigger than lapack NB */
    lwork = N*NB; 

    /* Allocate host memory for the matrix */
    TESTING_MALLOC(    h_A,    double, lda*N );
    TESTING_HOSTALLOC( h_R,    double, lda*N );
    TESTING_HOSTALLOC( h_work, double, lwork );
    TESTING_MALLOC(    tau,    double, N-1   );
    TESTING_HOSTALLOC( D,    double, N );
    TESTING_HOSTALLOC( E,    double, N );
    //TESTING_DEVALLOC( dT1,  double, (2*min(N,N)+(N+31)/32*32)*NB );
    TESTING_DEVALLOC( dT1,  double, (N*NB) );

    printf("  N    GPU GFlop/s   \n");
    printf("=====================\n");
    for(i=0; i<10; i++){
        if ( !once ) {
            N = size[i];
        }
        lda  = N;
        n2   = N*lda;
        flops = FLOPS( (double)N ) / 1e6;
        if(WANTZ) flops = 2.0*flops;

        /* ====================================================================
           Initialize the matrix
           =================================================================== */
        lapackf77_dlarnv( &ione, ISEED, &n2, h_A );

        // Make the matrix hermitian 
        {
            magma_int_t i, j;
            for(i=0; i<N; i++) {
                MAGMA_D_SET2REAL( h_A[i*lda+i], ( MAGMA_D_REAL(h_A[i*lda+i]) ) );
                for(j=0; j<i; j++)
                    h_A[i*lda+j] = (h_A[j*lda+i]);
            }
        }
/*
            for(i=0; i<N; i++){ 
                for(j=0; j<N; j++){
                MAGMA_D_SET2REAL( h_A[i*lda+j], ( MAGMA_D_REAL(h_A[i*lda+j]) ) );
                }
            }
*/

/*
    FILE *trace_file;
    trace_file = fopen("AJETE/Ainit", "w");
    for (j = 0; j < N ; j++) 
          for (i = 0; i < N ; i++) 
                         fprintf(trace_file,"%10d %10d %25.15e %25.15e\n",i+1,j+1,MAGMA_D_REAL(h_A[j*lda+i]) ,  MAGMA_D_IMAG(h_A[j*lda+i])  );
    fclose(trace_file);
*/



        lapackf77_dlacpy( MagmaUpperLowerStr, &N, &N, h_A, &lda, h_R, &lda );

/*
lapackf77_dlarnv( &ione, ISEED, &N, D );
lapackf77_dlarnv( &ione, ISEED, &N, E );
i= min(12,THREADS);
mkl_set_num_threads( i );
start = get_current_time();
dstedc_withZ('V', N, D, E, h_R, lda);
end = get_current_time();
printf("  Finish EIGEN   timing= %lf  threads %d ---> 0000000 \n" ,GetTimerValue(start,end) / 1000., i);
mkl_set_num_threads( 1 );
return 0;
*/
/*

    FILE *trace_file;
    trace_file = fopen("AJETE/Ainit", "w");
    for (j = 0; j < N ; j++) 
          for (i = 0; i < N ; i++) 
                         fprintf(trace_file,"%10d %10d %25.15e %25.15e\n",i+1,j+1,MAGMA_D_REAL(h_R[j*lda+i]) ,  MAGMA_D_IMAG(h_R[j*lda+i])  );
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
              lapackf77_dgeqrf(&pm, &pn, &h_R[NB], &lda, 
                             tau, h_work, &lwork, &info);
              printf("TOTOTOTO INFO %d\n",info);
              memset(h_work, 0, lwork*sizeof(double));
            lapackf77_dlarft( "F", "C",
                              &pm, &pn, &h_R[NB], &lda,
                              tau, h_work, &pn);

*/


        /* ====================================================================
           Performs operation using MAGMA
           =================================================================== */
        start = get_current_time();
        //magma_dsytrd_sy2sb(uplo[0], N, h_R, lda, tau, h_work, lwork, &info);
        magma_dsytrd_sy2sb(uplo[0], N, NB, h_R, lda, tau, h_work, lwork, dT1, THREADS, &info);
        end = get_current_time();
        printf("  Finish BAND    timing= %lf \n" ,GetTimerValue(start,end) / 1000.);
// return 0;
/*        
    int   Vblksiz=-1, blkcnt=-1, LDV=-1, LDT =-1, INgrsiz=1, LDE=-1, BAND=6;
    Vblksiz = NB; //min(NB,64);
    LDT     = Vblksiz;
    findVTsiz(N, NB, Vblksiz, &blkcnt, &LDV);
    double *dVV2, *dTT2, dV3;
    int dVVs;
           dVVs = max(N*N,blkcnt*LDV*Vblksiz);
           printf("dvsize %lf \n",(16.0*(real_Double_t)dVVs)*1e-9);
           if( CUBLAS_STATUS_SUCCESS != cublasAlloc(dVVs, sizeof(double), (void**)&dVV2) ) { 
               printf ("!!!! -------> cublasAlloc failed for: dVV2\n" );       
               exit(-1);                                                           
           }
    
           if( CUBLAS_STATUS_SUCCESS != cublasAlloc( dVVs, sizeof(double), (void**)&dTT2) ) { 
              printf ("!!!! ---------> cublasAlloc failed for: dTT2\n" );       
              exit(-1);                                                           
           }
    
           if( CUBLAS_STATUS_SUCCESS != cublasAlloc( dVVs, sizeof(double), (void**)&dV3) ) { 
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
        memset(h_work, 0, lwork*sizeof(double));
        magma_dgetmatrix( pn, pn, dT1, ldt, h_work, pn );
   trace_file = fopen("AJETE/T", "w");
    for (j = 0; j < pn ; j++) 
          for (i = 0; i < pn ; i++) 
                         fprintf(trace_file,"%10d%10d%40.30e\n",i+1,j+1,h_work[j*pn+i]);
    fclose(trace_file);
*/        
        
        //        dsytrd_bsy2trc(THREADS, uplo[0], N, NB, h_R, lda, D, E);
        magma_dsytrd_bsy2trc_v5(THREADS, WANTZ, uplo[0], NE, N, NB, h_R, lda, D, E, dT1, ldt);
        end = get_current_time();
        if ( info < 0 )
            printf("Argument %d of magma_dsytrd_sy2sb had an illegal value\n", (int) -info);

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
            double nrmI=0.0, nrm1=0.0, nrm2=0.0;
            int    lwork2 = 256*N;
            double *work2     = (double *) malloc (lwork2*sizeof(double));
            double *rwork2     = (double *) malloc (N*sizeof(double));
            double *D2          = (double *) malloc (N*sizeof(double));
            double *AINIT    = (double *) malloc (N*lda*sizeof(double));
            memcpy(AINIT, h_A, N*lda*sizeof(double));
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
            lapackf77_dsyev( "N", "L", &N, h_A, &lda, D2, work2, &lwork2, rwork2, &info );
#else
            lapackf77_dsyev( "N", "L", &N, h_A, &lda, D2, work2, &lwork2, &info );
#endif
            ///* call eigensolver for our resulting tridiag [D E] and for Q */
            //dstedc_withZ('V', N, D, E, h_R, lda);
            ////dsterf_( &N, D, E, &info); 
            ////
            end = get_current_time();
            printf("  Finish CHECK - EIGEN   timing= %lf  threads %d \n" ,GetTimerValue(start,end) / 1000., i);

            /*
            for(i=0;i<10;i++)
                printf(" voici lpk D[%d] %e\n",i,D2[i]);
            */

            //double mydz=0.0,mydo=1.0;
            //double *Z = (double *) malloc(N*lda*sizeof(double));
           // dgemm_("N","N",&N,&N,&N,&mydo,h_R,&lda,h_A,&lda,&mydz,Z,&lda);


           /* compare result */
           cmp_vals(N, D2, D, &nrmI, &nrm1, &nrm2);


           double *WORKAJETER;
           double *RWORKAJETER, *RESU;
           WORKAJETER  = (double *) malloc( (2* N * N + N) * sizeof(double) );
           RWORKAJETER = (double *) malloc( N * sizeof(double) );
           RESU        = (double *) malloc(10*sizeof(double));
           int MATYPE;
           memset(RESU,0,10*sizeof(double));

 
           MATYPE=3;
           double NOTHING=0.0;
           start = get_current_time();
           // check results
           dcheck_eig_(&JOBZ, &MATYPE, &N, &NB, AINIT, &lda, &NOTHING, &NOTHING, D2 , D, h_R, &lda, WORKAJETER, RWORKAJETER, RESU );
           end = get_current_time();
           printf("  Finish CHECK - results timing= %lf \n" ,GetTimerValue(start,end) / 1000.);
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
