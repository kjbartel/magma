/*
 * Copyright (c) 2011      The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 *
 *
 *     @author Azzam Haidar
 *     @author Stan Tomov
 *
 *     @generated s Tue May 15 18:17:54 2012
 *
 */

#include "common_magma.h"
#include "magma_sbulgeinc.h"
// === Define what BLAS to use ============================================

// === End defining what BLAS to use ======================================
 

//////////////////////////////////////////////////////////////
//          DSTEDC          Divide and Conquer for tridiag
//////////////////////////////////////////////////////////////
extern "C" void  magma_sstedc_withZ(char JOBZ, magma_int_t N, float *D, float * E, float *Z, magma_int_t LDZ) {
  float *WORK;
  magma_int_t *IWORK;
  magma_int_t LWORK, LIWORK;
  magma_int_t INFO;
  magma_int_t NxN=N*N;
   
  if(JOBZ=='V'){
        LWORK  = 1 + 3*N + 3*N*((magma_int_t)log2(N)+1) + 4*N*N+ 256*N; 
        LIWORK =  6 + 6*N + 6*N*((magma_int_t)log2(N)+1) + 256*N;
  }else if(JOBZ=='I'){
        LWORK  = 2*N*N+256*N+1; 
          LIWORK = 256*N;
  }else if(JOBZ=='N'){
        LWORK  = 256*N+1; 
          LIWORK = 256*N;  
  }else{
          printf("ERROR JOBZ %c\n",JOBZ);
          exit(-1);
  }

  WORK = (float*) malloc( LWORK*sizeof( float) );
  IWORK = (magma_int_t*) malloc( LIWORK*sizeof( magma_int_t) );

  lapackf77_sstedc(&JOBZ, &N, D, E, Z, &LDZ, WORK,&LWORK,IWORK,&LIWORK,&INFO);

  if(INFO!=0){
        printf("=================================================\n");
        printf("DSTEDC ERROR OCCURED. HERE IS INFO %d \n ",INFO);
        printf("=================================================\n");
          //assert(INFO==0);
  }


  free( IWORK );
  free( WORK );
}
//////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////
//          DSTEDX          Divide and Conquer for tridiag
//////////////////////////////////////////////////////////////
extern "C" void  magma_sstedx_withZ(magma_int_t N, magma_int_t NE, float *D, float * E, float *Z, magma_int_t LDZ) {
  float *WORK;
  float *dwork;
  magma_int_t *IWORK;
  magma_int_t LWORK, LIWORK;
  magma_int_t INFO;
   
  LWORK  = N*N+4*N+1; 
  LIWORK = 3 + 5*N;

  WORK = (float*) malloc( LWORK*sizeof( float) );
  IWORK = (magma_int_t*) malloc( LIWORK*sizeof( magma_int_t) );

  if (MAGMA_SUCCESS != magma_smalloc( &dwork, 3*N*(N/2 + 1) )) {
     printf("=================================================\n");
     printf("DSTEDC ERROR OCCURED IN CUDAMALLOC\n");
     printf("=================================================\n");
     return;
  }
  printf("using magma_sstedx\n");

  magma_sstedx('I', N, 0., 0., 1, NE, D, E, Z, LDZ, WORK,LWORK,IWORK,LIWORK,dwork,&INFO);

  if(INFO!=0){
        printf("=================================================\n");
        printf("DSTEDC ERROR OCCURED. HERE IS INFO %d \n ",INFO);
        printf("=================================================\n");
          //assert(INFO==0);
  }

  magma_free( dwork );
  free( IWORK );
  free( WORK );
}
//////////////////////////////////////////////////////////////
