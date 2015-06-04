/*
 * Copyright (c) 2011      The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 *
 *
 *     @author Azzam Haidar
 *     @author Stan Tomov
 *
 *     @precisions normal z -> c
 *
 */

#include "common_magma.h"
//#include "magma_zbulgeinc.h"
// === Define what BLAS to use ============================================

// === End defining what BLAS to use ======================================
 

//////////////////////////////////////////////////////////////
//          ZSTEDC          Divide and Conquer for tridiag
//////////////////////////////////////////////////////////////
extern "C" void  magma_zstedc_withZ(char JOBZ, magma_int_t N, double *D, double * E, cuDoubleComplex *Z, magma_int_t LDZ) {
  cuDoubleComplex *WORK;
  double *RWORK;
  magma_int_t *IWORK;
  magma_int_t LWORK, LIWORK, LRWORK;
  magma_int_t INFO;
  magma_int_t NxN=N*N;
   
  if(JOBZ=='V'){
      LWORK = N*N;
      LRWORK  = 1 + 3*N + 3*N*((magma_int_t)log2(N)+1) + 4*N*N+ 256*N; 
      LIWORK =  6 + 6*N + 6*N*((magma_int_t)log2(N)+1) + 256*N;
  }else if(JOBZ=='I'){
      LWORK = N;
      LRWORK  = 2*N*N+4*N+1+256*N; 
      LIWORK = 256*N;
  }else if(JOBZ=='N'){
      LWORK = N;
      LRWORK  = 256*N+1; 
      LIWORK = 256*N;  
  }else{
      printf("ERROR JOBZ %c\n",JOBZ);
      exit(-1);
  }

  RWORK  = (double*) malloc( LRWORK*sizeof( double) );
  WORK   = (cuDoubleComplex*) malloc( LWORK*sizeof( cuDoubleComplex) );
  IWORK  = (magma_int_t*) malloc( LIWORK*sizeof( magma_int_t) );

  lapackf77_zstedc(&JOBZ, &N, D, E, Z, &LDZ, WORK, &LWORK, RWORK, &LRWORK, IWORK, &LIWORK, &INFO);

  if(INFO!=0){
        printf("=================================================\n");
        printf("ZSTEDC ERROR OCCURED. HERE IS INFO %d \n ", (int) INFO);
        printf("=================================================\n");
          //assert(INFO==0);
  }


  free( IWORK );
  free( WORK );
  free( RWORK );
}
//////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////
//          ZSTEDC          Divide and Conquer for tridiag
//////////////////////////////////////////////////////////////
extern "C" void  magma_zstedx_withZ(magma_int_t N, magma_int_t NE, double *D, double * E, cuDoubleComplex *Z, magma_int_t LDZ) {
  double *RWORK;
  double *dwork;
  magma_int_t *IWORK;
  magma_int_t LWORK, LIWORK, LRWORK;
  magma_int_t INFO;
  magma_int_t NxN=N*N;
   
      LWORK = N;
      LRWORK  = 2*N*N+4*N+1+256*N; 
      LIWORK = 256*N;

  RWORK  = (double*) malloc( LRWORK*sizeof( double) );
  IWORK  = (magma_int_t*) malloc( LIWORK*sizeof( magma_int_t) );

  if (MAGMA_SUCCESS != magma_dmalloc( &dwork, 3*N*(N/2 + 1) )) {
     printf("=================================================\n");
     printf("ZSTEDC ERROR OCCURED IN CUDAMALLOC\n");
     printf("=================================================\n");
     return;
  }
  printf("using magma_zstedx\n");

  magma_zstedx('I', N, 0.,0., 1, NE, D, E, Z, LDZ, RWORK, LRWORK, IWORK, LIWORK, dwork, &INFO);

  if(INFO!=0){
        printf("=================================================\n");
        printf("ZSTEDC ERROR OCCURED. HERE IS INFO %d \n ", (int) INFO);
        printf("=================================================\n");
          //assert(INFO==0);
  }

  magma_free( dwork );
  free( IWORK );
  free( RWORK );
}

