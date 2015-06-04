/*
 * Copyright (c) 2011      The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 *
 *
 *     @author Azzam Haidar
 *     @author Stan Tomov
 *
 *     @generated c Tue May 15 18:17:54 2012
 *
 */

#include "common_magma.h"
//#include "magma_cbulgeinc.h"
// === Define what BLAS to use ============================================

// === End defining what BLAS to use ======================================
 

//////////////////////////////////////////////////////////////
//          DSTEDC          Divide and Conquer for tridiag
//////////////////////////////////////////////////////////////
extern "C" void  magma_cstedc_withZ(char JOBZ, magma_int_t N, float *D, float * E, cuFloatComplex *Z, magma_int_t LDZ) {
  cuFloatComplex *WORK;
  float *RWORK;
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

  RWORK  = (float*) malloc( LRWORK*sizeof( float) );
  WORK   = (cuFloatComplex*) malloc( LWORK*sizeof( cuFloatComplex) );
  IWORK  = (magma_int_t*) malloc( LIWORK*sizeof( magma_int_t) );

  lapackf77_cstedc(&JOBZ, &N, D, E, Z, &LDZ, WORK, &LWORK, RWORK, &LRWORK, IWORK, &LIWORK, &INFO);

  if(INFO!=0){
        printf("=================================================\n");
        printf("DSTEDC ERROR OCCURED. HERE IS INFO %d \n ",INFO);
        printf("=================================================\n");
          //assert(INFO==0);
  }


  free( IWORK );
  free( WORK );
  free( RWORK );
}
//////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////
//          DSTEDC          Divide and Conquer for tridiag
//////////////////////////////////////////////////////////////
extern "C" void  magma_cstedx_withZ(magma_int_t N, magma_int_t NE, float *D, float * E, cuFloatComplex *Z, magma_int_t LDZ) {
  float *RWORK;
  float *dwork;
  magma_int_t *IWORK;
  magma_int_t LWORK, LIWORK, LRWORK;
  magma_int_t INFO;
  magma_int_t NxN=N*N;
   
      LWORK = N;
      LRWORK  = 2*N*N+4*N+1+256*N; 
      LIWORK = 256*N;

  RWORK  = (float*) malloc( LRWORK*sizeof( float) );
  IWORK  = (magma_int_t*) malloc( LIWORK*sizeof( magma_int_t) );

  if (MAGMA_SUCCESS != magma_smalloc( &dwork, 3*N*(N/2 + 1) )) {
     printf("=================================================\n");
     printf("ZSTEDC ERROR OCCURED IN CUDAMALLOC\n");
     printf("=================================================\n");
     return;
  }
  printf("using magma_cstedx\n");

  magma_cstedx('I', N, 0.,0., 1, NE, D, E, Z, LDZ, RWORK, LRWORK, IWORK, LIWORK, dwork, &INFO);

  if(INFO!=0){
        printf("=================================================\n");
        printf("ZSTEDC ERROR OCCURED. HERE IS INFO %d \n ",INFO);
        printf("=================================================\n");
          //assert(INFO==0);
  }

  magma_free( dwork );
  free( IWORK );
  free( RWORK );
}

