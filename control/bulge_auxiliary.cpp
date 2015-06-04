/*
 * Copyright (c) 2011      The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 *
 *
 *     @author Azzam Haidar
 *     @author Stan Tomov
 *
 *
 */

#include "common_magma.h"


#ifdef __cplusplus
extern "C" {
#endif

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
#include <sys/time.h>

 real_Double_t get_time_azz(void)
{
    struct timeval  time_val;
    struct timezone time_zone;

    gettimeofday(&time_val, &time_zone);

    return (real_Double_t)(time_val.tv_sec) + (real_Double_t)(time_val.tv_usec) / 1000000.0;
}
////////////////////////////////////////////////////////////////////////////////////////////////////
        
////////////////////////////////////////////////////////////////////////////////////////////////////
 magma_int_t plasma_ceildiv(magma_int_t a, magma_int_t b)
{
  real_Double_t r = (real_Double_t)a/(real_Double_t)b;
  r = (r-(magma_int_t)r)==0? (magma_int_t)r:(magma_int_t)r+1;
  return (magma_int_t) r;
}
////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////
 void cmp_vals(int n, double *wr1, double *wr2, double *nrmI, double *nrm1, double *nrm2) {
  int i;
  double curv, maxv, sumv;

  maxv = 0.0;
  sumv = 0.0;
  for (i = 0; i < n; ++i) {

    curv = fabs( wr1[i] - wr2[i]);
    sumv += curv;
    if (maxv < curv) maxv = curv;
  }

    *nrmI = maxv;
    *nrm1 = sumv;
    *nrm2 = sqrt( sumv );
}
//////////////////////////////////////////////////////////////////////////////////////////////////// 




////////////////////////////////////////////////////////////////////////////////////////////////////
 void findVTpos(magma_int_t N, magma_int_t NB, magma_int_t Vblksiz, magma_int_t sweep, magma_int_t st, magma_int_t *Vpos, magma_int_t *TAUpos, magma_int_t *Tpos, magma_int_t *myblkid)
{
  magma_int_t prevcolblknb, prevblkcnt, prevcolblkid;
  magma_int_t curcolblknb, nbprevcolblk, mastersweep;
  magma_int_t blkid, locj, LDV;
  prevcolblknb = 0;
  prevblkcnt   = 0;
  curcolblknb  = 0;

  nbprevcolblk = sweep/Vblksiz;
  for (prevcolblkid = 0; prevcolblkid < nbprevcolblk; prevcolblkid++)
  {
       mastersweep  = prevcolblkid * Vblksiz;
       prevcolblknb = plasma_ceildiv((N-(mastersweep+2)),NB);
       prevblkcnt   = prevblkcnt + prevcolblknb;
  }
  curcolblknb = plasma_ceildiv((st-sweep),NB);
  blkid       = prevblkcnt + curcolblknb -1;
  locj        = sweep%Vblksiz;
  LDV         = NB + Vblksiz -1;

  *myblkid= blkid;
  *Vpos   = blkid*Vblksiz*LDV  + locj*LDV + locj;
  *TAUpos = blkid*Vblksiz + locj;
  *Tpos   = blkid*Vblksiz*Vblksiz + locj*Vblksiz + locj;
  //printf("voici  blkid  %d  locj %d  vpos %d tpos %d \n",blkid,locj,*Vpos,*Tpos);
}
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
 void findVTsiz(magma_int_t N, magma_int_t NB, magma_int_t Vblksiz, magma_int_t *blkcnt, magma_int_t *LDV)
{
  magma_int_t colblk, nbcolblk;
  magma_int_t curcolblknb, mastersweep;

  *blkcnt   = 0;
  nbcolblk = plasma_ceildiv((N-1),Vblksiz);
  for (colblk = 0; colblk<nbcolblk; colblk++)
  {
        mastersweep = colblk * Vblksiz;
        curcolblknb = plasma_ceildiv((N-(mastersweep+2)),NB);
        *blkcnt      = *blkcnt + curcolblknb;
        //printf("voici  nbcolblk %d    master sweep %d     blkcnt %d \n",nbcolblk, mastersweep,*blkcnt); 
  }
  *blkcnt = *blkcnt +1;
  *LDV= NB+Vblksiz-1;
}
////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __cplusplus
}
#endif


