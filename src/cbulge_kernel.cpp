/*
 * Copyright (c) 2011      The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 *
 *
 *     @author Azzam Haidar
 *     @author Stan Tomov
 *
 *     @generated c Tue Dec 17 13:18:36 2013
 *
 */

#include "common_magma.h"
#include <cblas.h>

#define PRECISION_c
 
#ifdef __cplusplus
extern "C" {
#endif

void findVTpos(magma_int_t N, magma_int_t NB, magma_int_t Vblksiz, magma_int_t sweep, magma_int_t st, magma_int_t *Vpos, magma_int_t *TAUpos, magma_int_t *Tpos, magma_int_t *myblkid);

void findVTsiz(magma_int_t N, magma_int_t NB, magma_int_t Vblksiz, magma_int_t *blkcnt, magma_int_t *LDV);

magma_int_t plasma_ceildiv(magma_int_t a, magma_int_t b);

void magma_ctrdtype1cbHLsym_withQ(magma_int_t N, magma_int_t NB, 
                                magmaFloatComplex *A, magma_int_t LDA, magmaFloatComplex *V, magmaFloatComplex *TAU, 
                                magma_int_t st, magma_int_t ed, magma_int_t sweep, magma_int_t Vblksiz);

void magma_ctrdtype2cbHLsym_withQ(magma_int_t N, magma_int_t NB, magmaFloatComplex *A, magma_int_t LDA, magmaFloatComplex *V, magmaFloatComplex *TAU, magma_int_t st, magma_int_t ed, magma_int_t sweep, magma_int_t Vblksiz);
   
void magma_ctrdtype3cbHLsym_withQ(magma_int_t N, magma_int_t NB, magmaFloatComplex *A, magma_int_t LDA, magmaFloatComplex *V, magmaFloatComplex *TAU, magma_int_t st, magma_int_t ed, magma_int_t sweep, magma_int_t Vblksiz);

void magma_clarfxsym(magma_int_t N, magmaFloatComplex *A, magma_int_t LDA, magmaFloatComplex *V, magmaFloatComplex *TAU);

#ifdef __cplusplus
}
#endif


////////////////////////////////////////////////////////////////////////////////////////////////////

extern "C" void 
magma_clarfxsym(magma_int_t N, magmaFloatComplex *A, magma_int_t LDA, magmaFloatComplex *V, magmaFloatComplex *TAU) {
  magma_int_t IONE=1; 
  magmaFloatComplex dtmp;
  magmaFloatComplex Z_ZERO =  MAGMA_C_ZERO;
  //magmaFloatComplex Z_ONE  =  MAGMA_C_ONE;
  magmaFloatComplex Z_MONE =  MAGMA_C_NEG_ONE;
  magmaFloatComplex Z_HALF =  MAGMA_C_HALF;
  //magmaFloatComplex WORK[N];
  magmaFloatComplex *WORK  = (magmaFloatComplex *) malloc( N * sizeof(magmaFloatComplex) );

  /* apply left and right on A(st:ed,st:ed)*/
  //magma_clarfxsym(len,A(st,st),LDX,V(st),TAU(st));
  /* X = AVtau */
  blasf77_chemv("L",&N, TAU, A, &LDA, V, &IONE, &Z_ZERO, WORK, &IONE);
  /* je calcul dtmp= X'*V */
#if defined(PRECISION_z) || defined(PRECISION_c)
   dtmp = Z_ZERO;
   for (magma_int_t j = 0; j < N ; j++)
      dtmp = dtmp + MAGMA_C_CNJG(WORK[j]) * V[j];
   //cblas_cdotc_sub(N, WORK, IONE, V, IONE, &dtmp);
#else
  dtmp = cblas_cdotc(N, WORK, IONE, V, IONE);
#endif
  /* je calcul 1/2 X'*V*t = 1/2*dtmp*tau  */
  dtmp = -dtmp * Z_HALF * (*TAU);
  /* je calcul W=X-1/2VX'Vt = X - dtmp*V */
  /*
  for (j = 0; j < N ; j++)
      WORK[j] = WORK[j] + (dtmp*V[j]); */
  blasf77_caxpy(&N, &dtmp, V, &IONE, WORK, &IONE);
  /* performs the symmetric rank 2 operation A := alpha*x*y' + alpha*y*x' + A */
  blasf77_cher2("L",&N,&Z_MONE,WORK,&IONE,V,&IONE,A,&LDA);
  
  magma_free_cpu(WORK);
}
////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////
//                  TYPE 1-BAND Householder
///////////////////////////////////////////////////////////
//// add -1 because of C
#define A(m,n)   &(A[((m)-(n)) + LDA*((n)-1)])
#define V(m)     &(V[(m)])
#define TAU(m)   &(TAU[(m)])
extern "C" void magma_ctrdtype1cbHLsym_withQ(magma_int_t N, magma_int_t NB, magmaFloatComplex *A, magma_int_t LDA, magmaFloatComplex *V, magmaFloatComplex *TAU, magma_int_t st, magma_int_t ed, magma_int_t sweep, magma_int_t Vblksiz) {
  //magma_int_t    J1, J2, J3, i, j;
  magma_int_t    len, LDX;
  magma_int_t    IONE=1;
  magma_int_t    blkid, vpos, taupos, tpos; 
  //magmaFloatComplex conjtmp;
  magmaFloatComplex Z_ONE  =  MAGMA_C_ONE;
  magmaFloatComplex *WORK  = (magmaFloatComplex *) malloc( N * sizeof(magmaFloatComplex) );


  findVTpos(N,NB,Vblksiz,sweep-1,st-1, &vpos, &taupos, &tpos, &blkid);
  //printf("voici vpos %d taupos %d  tpos %d  blkid %d \n", vpos, taupos, tpos, blkid);
  LDX     = LDA-1;
  len     = ed-st+1;
  *V(vpos)  = Z_ONE;
  memcpy(V(vpos+1), A(st+1, st-1), (len-1)*sizeof(magmaFloatComplex));
  memset(A(st+1, st-1), 0, (len-1)*sizeof(magmaFloatComplex));
  /* Eliminate the col  at st-1 */
  lapackf77_clarfg( &len, A(st, st-1), V(vpos+1), &IONE, TAU(taupos) );
  /* apply left and right on A(st:ed,st:ed)*/
  magma_clarfxsym(len,A(st,st),LDX,V(vpos),TAU(taupos));
  //conjtmp = MAGMA_C_CNJG(*TAU(taupos));
  //lapackf77_clarfy("L", &len, V(vpos), &IONE, &conjtmp, A(st,st), &LDX, WORK); //&(MAGMA_C_CNJG(*TAU(taupos)))
  magma_free_cpu(WORK);
}
#undef A
#undef V
#undef TAU
///////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////
//                  TYPE 1-LPK Householder
///////////////////////////////////////////////////////////
//// add -1 because of C
#define A(m,n)   &(A[((m)-(n)) + LDA*((n)-1)])
#define V(m)     &(V[(m)])
#define TAU(m)   &(TAU[(m)])
extern "C" void magma_ctrdtype2cbHLsym_withQ(magma_int_t N, magma_int_t NB, magmaFloatComplex *A, magma_int_t LDA, magmaFloatComplex *V, magmaFloatComplex *TAU, magma_int_t st, magma_int_t ed, magma_int_t sweep, magma_int_t Vblksiz) {
  magma_int_t    J1, J2, len, lem, LDX;
  //magma_int_t    i, j;
  magma_int_t    IONE=1;
  magma_int_t    blkid, vpos, taupos, tpos; 
  magmaFloatComplex conjtmp;
  magmaFloatComplex Z_ONE  =  MAGMA_C_ONE;
  //magmaFloatComplex WORK[NB];
  magmaFloatComplex *WORK  = (magmaFloatComplex *) malloc( NB * sizeof(magmaFloatComplex) );


  findVTpos(N,NB,Vblksiz,sweep-1,st-1, &vpos, &taupos, &tpos, &blkid);
  LDX    = LDA-1;
  J1     = ed+1;
  J2     = min(ed+NB,N);
  len    = ed-st+1;
  lem    = J2-J1+1;
  if(lem>0){
     /* apply remaining right commming from the top block */
     lapackf77_clarfx("R", &lem, &len, V(vpos), TAU(taupos), A(J1, st), &LDX, WORK);
  }
  if(lem>1){
     findVTpos(N,NB,Vblksiz,sweep-1,J1-1, &vpos, &taupos, &tpos, &blkid);
     /* remove the first column of the created bulge */
     *V(vpos)  = Z_ONE;
     memcpy(V(vpos+1), A(J1+1, st), (lem-1)*sizeof(magmaFloatComplex));
     memset(A(J1+1, st),0,(lem-1)*sizeof(magmaFloatComplex));
     /* Eliminate the col at st */
     lapackf77_clarfg( &lem, A(J1, st), V(vpos+1), &IONE, TAU(taupos) );
     /* apply left on A(J1:J2,st+1:ed) */
     len = len-1; /* because we start at col st+1 instead of st. col st is the col that has been revomved;*/
     conjtmp = MAGMA_C_CNJG(*TAU(taupos));
     lapackf77_clarfx("L", &lem, &len, V(vpos),  &conjtmp, A(J1, st+1), &LDX, WORK);
  }
  free (WORK);
}
#undef A
#undef V
#undef TAU
///////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////
//                  TYPE 1-LPK Householder
///////////////////////////////////////////////////////////
//// add -1 because of C
#define A(m,n)   &(A[((m)-(n)) + LDA*((n)-1)])
#define V(m)     &(V[(m)])
#define TAU(m)   &(TAU[(m)])
extern "C" void magma_ctrdtype3cbHLsym_withQ(magma_int_t N, magma_int_t NB, magmaFloatComplex *A, magma_int_t LDA, magmaFloatComplex *V, magmaFloatComplex *TAU, magma_int_t st, magma_int_t ed, magma_int_t sweep, magma_int_t Vblksiz) {
  //magma_int_t    J1, J2, J3, i, j;
  magma_int_t    len, LDX;
  //magma_int_t    IONE=1;
  magma_int_t    blkid, vpos, taupos, tpos; 
  //magmaFloatComplex conjtmp;
  magmaFloatComplex *WORK  = (magmaFloatComplex *) malloc( N * sizeof(magmaFloatComplex) );


  findVTpos(N,NB,Vblksiz,sweep-1,st-1, &vpos, &taupos, &tpos, &blkid);
  LDX    = LDA-1;
  len    = ed-st+1;

  /* apply left and right on A(st:ed,st:ed)*/
  magma_clarfxsym(len,A(st,st),LDX,V(vpos),TAU(taupos));
  //conjtmp = MAGMA_C_CNJG(*TAU(taupos));
  //lapackf77_clarfy("L", &len, V(vpos), &IONE,  &(MAGMA_C_CNJG(*TAU(taupos))), A(st,st), &LDX, WORK);
  magma_free_cpu(WORK);
}
#undef A
#undef V
#undef TAU
///////////////////////////////////////////////////////////





