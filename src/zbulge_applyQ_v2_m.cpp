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
 *     @precisions normal z -> s d c
 *
 */

#include "common_magma.h"
#include "magma_bulge.h"
#include <cblas.h>
// === Define what BLAS to use ============================================
#define PRECISION_z

// === End defining what BLAS to use ======================================

////////////////////////////////////////////////////////////////////////////////////////////////////

#define dE(ind,i,j)  (dE[ind]+(i) + ldde*(j))
#define V(j)     (V+(j))
#define T(j)     (T+(j))

extern "C" magma_int_t magma_zbulge_applyQ_v2_m(magma_int_t nrgpu, char side, magma_int_t NE, magma_int_t N, magma_int_t NB, magma_int_t Vblksiz, cuDoubleComplex *E, magma_int_t lde, cuDoubleComplex *V, magma_int_t ldv, cuDoubleComplex *T, magma_int_t ldt, magma_int_t *info)
{
    //%===========================
    //%   local variables
    //%===========================
    magma_int_t rownbm, m, n;
    magma_int_t st,ed,fst,vlen,vnb,colj;
    magma_int_t vpos,tpos;
    magma_int_t cur_blksiz,avai_blksiz, ncolinvolvd;
    magma_int_t nbgr, colst, coled;

    *info=0;

    cuDoubleComplex *dE[MagmaMaxGPUs], *dwork[MagmaMaxGPUs], *dT[MagmaMaxGPUs], *dV[MagmaMaxGPUs];
    magma_int_t ldde = N;
    magma_int_t lddv = NB + Vblksiz + 1;
    magma_int_t lddt = Vblksiz;
    magma_int_t ldwork =NE;

    magma_int_t ne_loc = magma_ceildiv(NE, nrgpu);

    if (ne_loc<256)
       ne_loc=256;

    nrgpu = min(nrgpu, magma_ceildiv(NE,ne_loc)); // Don't use GPU that will not have data.

    if (N==0 || NE==0){
        return MAGMA_SUCCESS;
    }

    int gpu_b;
    magma_getdevice(&gpu_b);

    for (magma_int_t igpu=0; igpu < nrgpu; ++igpu){

        magma_setdevice(igpu);

        if(MAGMA_SUCCESS != magma_zmalloc( &dE[igpu], ldde * ne_loc)) {
            printf ("!!!!  magma_zbulge_applyQ magma_alloc failed for: dE\n" );
            exit(-1);
        }
        if(MAGMA_SUCCESS != magma_zmalloc( &dwork[igpu], 2*ldwork*Vblksiz + Vblksiz * (lddv+lddt))) {
            printf ("!!!!  magma_zbulge_applyQ magma_alloc failed for: dwork\n" );
            exit(-1);
        }
        dV[igpu] = dwork[igpu] + 2*ldwork*Vblksiz;
        dT[igpu] = dV[igpu] + Vblksiz*lddv;

        magma_int_t ie_loc = min(ne_loc, NE - ne_loc*igpu);

        magma_zsetmatrix_async( N, ie_loc, E+lde*ne_loc*igpu, lde, dE(igpu, 0, 0), ldde, NULL );
    }

    magma_setdevice(gpu_b);

    /* SIDE LEFT  meaning apply E = Q*E = (q_1*q_2*.....*q_n) * E ==> so traverse Vs in reverse order (forward) from q_n to q_1
     *            Also E is splitten by row meaning each apply consist in a block of row (horizontal block) */
    /* SIDE RIGHT meaning apply E = E*Q = E * (q_1*q_2*.....*q_n) ==> so tarverse Vs in normal  order (forward) from q_1 to q_n
     *            Also E is splitten by col meaning each apply consist in a block of col (vertical block) */


    printf("  APPLY Q_v2_m nrGPU %d, with  N %d, NE %d,  NB %d, Vblksiz %d, SIDE %c\n",
           nrgpu, N, NE, NB, Vblksiz, side);

    if(side=='L'){
        rownbm    = magma_ceildiv((N-1),NB);
        for (m = rownbm; m>0; m--)
        {
           ncolinvolvd = min(N-1, m*NB);
           avai_blksiz=min(Vblksiz,ncolinvolvd);
           nbgr = magma_ceildiv(ncolinvolvd,avai_blksiz);
           for (n = nbgr; n>0; n--)
           {
               vlen = 0;
               vnb  = 0;
               cur_blksiz = min(ncolinvolvd-(n-1)*avai_blksiz, avai_blksiz);
               colst = (n-1)*avai_blksiz;
               coled = colst + cur_blksiz -1;
               fst   = (rownbm -m)*NB+colst +1;
               for (colj=colst; colj<=coled; colj++)
               {
                   st       = (rownbm -m)*NB+colj +1;
                   ed       = min(st+NB-1,N-1);
                   if(st>ed)break;
                   if((st==ed)&&(colj!=N-2))break;
                   vlen=ed-fst+1;
                   vnb=vnb+1;
               }
               if((vlen>0)&&(vnb>0)){
                   magma_bulge_findVTpos(N, NB, Vblksiz, colst, fst, ldv, ldt, &vpos, &tpos);

                   for (magma_int_t igpu=0; igpu < nrgpu; ++igpu){

                       magma_setdevice(igpu);

                       magma_int_t ie_loc = min(ne_loc, NE - ne_loc*igpu);

                       magma_zsetmatrix_async(vlen, vnb, V(vpos), ldv, dV[igpu], lddv, NULL);
                       magma_zsetmatrix_async(vnb,  vnb, T(tpos), ldt, dT[igpu], lddt, NULL);

                       // performance loss if the reflector are applied to a big number of eigenvectors (~10000)
                       // => apply the reflectors to blocks of eigenvectors.
                       magma_int_t nr_bl = magma_ceildiv(ie_loc,10000);       //nr of blocks
                       magma_int_t sz_bl = magma_ceildiv(ie_loc,nr_bl*64)*64; //maximum size of blocks (to have blocks of around the same size and multiple of 64)
                       magma_int_t ib;                                        //size of current block

                       for(magma_int_t i=0; i<ie_loc; i+= sz_bl){
                           ib = min(sz_bl, ie_loc-i);
                           magma_zlarfb_gpu( MagmaLeft, MagmaNoTrans, MagmaForward, MagmaColumnwise, vlen, ib, vnb, dV[igpu], lddv, dT[igpu], lddt, dE(igpu,fst,i), ldde, dwork[igpu], ldwork);
                       }
                   }
               }
           }
        }
    }else if (side=='R'){
        printf("Side 'R' not implemented in zbulge_applyQ_v2_m\n");
        exit(-1);
/*         rownbm    = magma_ceildiv((N-1),NB);
         for (m = 1; m<=rownbm; m++)
         {
            ncolinvolvd = min(N-1, m*NB);
            avai_blksiz=min(Vblksiz,ncolinvolvd);
            nbgr = magma_ceildiv(ncolinvolvd,avai_blksiz);
            for (n = 1; n<=nbgr; n++)
            {
                vlen = 0;
                vnb  = 0;
                cur_blksiz = min(ncolinvolvd-(n-1)*avai_blksiz, avai_blksiz);
                colst = (n-1)*avai_blksiz;
                coled = colst + cur_blksiz -1;
                fst   = (rownbm -m)*NB+colst +1;
                for (colj=colst; colj<=coled; colj++)
                {
                    st       = (rownbm -m)*NB+colj +1;
                    ed       = min(st+NB-1,N-1);
                    if(st>ed)break;
                    if((st==ed)&&(colj!=N-2))break;
                    vlen=ed-fst+1;
                    vnb=vnb+1;
                }
                magma_bulge_findVTpos(N, NB, Vblksiz ,colst, fst, ldv, ldt, &vpos, &tpos);
                magma_zsetmatrix_async(vlen, vnb, V(vpos), ldv, dV, lddv, NULL);
                magma_zsetmatrix_async(vnb,  vnb, T(tpos), ldt, dT, lddt, NULL);
                if((vlen>0)&&(vnb>0)){
                   magma_zlarfb_gpu( MagmaRight, MagmaNoTrans, MagmaForward, MagmaColumnwise, NE, vlen, vnb, dV, lddv, dT, lddt, dE(0, fst), ldde, dwork, ldwork);

               }
             }
         }*/
    }else{
            printf("ERROR SIDE %d \n",side);
    }

    for (magma_int_t igpu=0; igpu < nrgpu; ++igpu){

        magma_setdevice(igpu);

        magma_int_t ie_loc = min(ne_loc, NE - ne_loc*igpu);

        magma_zgetmatrix_async( N, ie_loc, dE(igpu, 0, 0), ldde, E+lde*ne_loc*igpu, lde, NULL );
    }
    for (magma_int_t igpu=0; igpu < nrgpu; ++igpu){

        magma_setdevice(igpu);

        magma_device_sync();

        magma_free(dwork[igpu]);
        magma_free(dE[igpu]);
    }

    magma_setdevice(gpu_b);

    return MAGMA_SUCCESS;
}
#undef E
#undef V
#undef TAU
#undef T
#undef dE
#undef dV
#undef dT
////////////////////////////////////////////////////////////////////////////////////////////////////



