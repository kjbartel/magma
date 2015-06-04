/*
 * Copyright (c) 2011      The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 *
 *
 *     @author Azzam Haidar
 *     @author Stan Tomov
 *     @author Raffaele Solca
 *
 */

#include "common_magma.h"
#include "magma_bulge.h"
#include <sys/time.h>

extern "C" {

    magma_int_t magma_bulge_get_nb(magma_int_t n)
    {
        return 64;
    }

    void cmp_vals(int n, double *wr1, double *wr2, double *nrmI, double *nrm1, double *nrm2)
    {
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

    void magma_bulge_findpos(magma_int_t n, magma_int_t nb, magma_int_t Vblksiz, magma_int_t sweep, magma_int_t st, magma_int_t *myblkid)
    {
        magma_int_t prevcolblknb, prevblkcnt, prevcolblkid;
        magma_int_t curcolblknb, nbprevcolblk, mastersweep;

        prevcolblknb = 0;
        prevblkcnt   = 0;
        curcolblknb  = 0;

        nbprevcolblk = sweep/Vblksiz;
        for (prevcolblkid = 0; prevcolblkid < nbprevcolblk; prevcolblkid++)
        {
            mastersweep  = prevcolblkid * Vblksiz;
            prevcolblknb = magma_ceildiv((n-(mastersweep+2)),nb);
            prevblkcnt   = prevblkcnt + prevcolblknb;
        }
        curcolblknb = magma_ceildiv((st-sweep),nb);
        *myblkid    = prevblkcnt + curcolblknb -1;

    }

    void magma_bulge_findVTAUpos(magma_int_t n, magma_int_t nb, magma_int_t Vblksiz, magma_int_t sweep, magma_int_t st, magma_int_t ldv,
                                 magma_int_t *Vpos, magma_int_t *TAUpos)
    {
        magma_int_t myblkid;
        magma_int_t locj = sweep%Vblksiz;

        magma_bulge_findpos(n, nb, Vblksiz, sweep, st, &myblkid);

        *Vpos   = myblkid*Vblksiz*ldv + locj*ldv + locj;
        *TAUpos = myblkid*Vblksiz + locj;
    }

    void magma_bulge_findVTpos(magma_int_t n, magma_int_t nb, magma_int_t Vblksiz, magma_int_t sweep, magma_int_t st, magma_int_t ldv, magma_int_t ldt,
                               magma_int_t *Vpos, magma_int_t *Tpos)
    {
        magma_int_t myblkid;
        magma_int_t locj = sweep%Vblksiz;

        magma_bulge_findpos(n, nb, Vblksiz, sweep, st, &myblkid);

        *Vpos   = myblkid*Vblksiz*ldv + locj*ldv + locj;
        *Tpos   = myblkid*Vblksiz*ldt + locj*ldt + locj;
    }

    void magma_bulge_findVTAUTpos(magma_int_t n, magma_int_t nb, magma_int_t Vblksiz, magma_int_t sweep, magma_int_t st, magma_int_t ldv, magma_int_t ldt,
                               magma_int_t *Vpos, magma_int_t *TAUpos, magma_int_t *Tpos, magma_int_t *blkid)
    {
        magma_int_t myblkid;
        magma_int_t locj = sweep%Vblksiz;

        magma_bulge_findpos(n, nb, Vblksiz, sweep, st, &myblkid);

        *Vpos   = myblkid*Vblksiz*ldv + locj*ldv + locj;
        *TAUpos = myblkid*Vblksiz     + locj;
        *Tpos   = myblkid*Vblksiz*ldt + locj*ldt + locj;
        *blkid  = myblkid;
    }

    magma_int_t magma_bulge_get_blkcnt(magma_int_t n, magma_int_t nb, magma_int_t Vblksiz)
    {
        magma_int_t colblk, nbcolblk;
        magma_int_t curcolblknb, mastersweep;

        magma_int_t blkcnt = 0;
        nbcolblk = magma_ceildiv((n-1),Vblksiz);
        for (colblk = 0; colblk<nbcolblk; colblk++)
        {
            mastersweep = colblk * Vblksiz;
            curcolblknb = magma_ceildiv((n-(mastersweep+2)),nb);
            blkcnt      = blkcnt + curcolblknb;
            //printf("voici  nbcolblk %d    master sweep %d     blkcnt %d \n",nbcolblk, mastersweep,*blkcnt);
        }
        return blkcnt +1;
    }

    ///////////////////
    // Old functions //
    ///////////////////

    magma_int_t plasma_ceildiv(magma_int_t a, magma_int_t b)
    {
        return magma_ceildiv(a,b);
    }

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


}



