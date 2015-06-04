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

#define dE(i,j)  (dE+(i) + ldde*(j))
#define V(j)     (V+(j))
#define T(j)     (T+(j))

extern "C" void magma_zbulge_applyQ_v2(char side, magma_int_t NE, magma_int_t N, magma_int_t NB, magma_int_t Vblksiz, cuDoubleComplex *dE, magma_int_t ldde, cuDoubleComplex *V, magma_int_t ldv, cuDoubleComplex *T, magma_int_t ldt, magma_int_t *info)
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

    cuDoubleComplex *dwork, *dT, *dV;
    magma_int_t lddv = NB + Vblksiz + 1;
    magma_int_t lddt = Vblksiz;
    magma_int_t ldwork;
    ldwork  = NE;
    if(MAGMA_SUCCESS != magma_zmalloc( &dwork, 2*ldwork*Vblksiz + Vblksiz * (lddv+lddt))) {
       printf ("!!!!  magma_zbulge_applyQ magma_alloc failed for: dwork\n" );
       exit(-1);
    }
    dV = dwork + 2*ldwork*Vblksiz;
    dT = dV + Vblksiz*lddv;


    /* SIDE LEFT  meaning apply E = Q*E = (q_1*q_2*.....*q_n) * E ==> so traverse Vs in reverse order (forward) from q_n to q_1
     *            Also E is splitten by row meaning each apply consist in a block of row (horizontal block) */
    /* SIDE RIGHT meaning apply E = E*Q = E * (q_1*q_2*.....*q_n) ==> so tarverse Vs in normal  order (forward) from q_1 to q_n
     *            Also E is splitten by col meaning each apply consist in a block of col (vertical block) */


    printf("  APPLY Q_v2 GPU with  N %d, NE %d,  NB %d, Vblksiz %d, SIDE %c\n",
           N, NE, NB, Vblksiz, side);

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
                   magma_zsetmatrix_async(vlen, vnb, V(vpos), ldv, dV, lddv, NULL);
                   magma_zsetmatrix_async(vnb,  vnb, T(tpos), ldt, dT, lddt, NULL);
                   //printf("voici bg %d m %d  vlen %d  vnb %d fcolj %d vpos %d taupos %d \n",bg,m,vlen, vnb,colst+1,vpos+1,taupos+1);

                   // performance loss if the reflector are applied to a big number of eigenvectors (~10000)
                   // => apply the reflectors to blocks of eigenvectors.
                   magma_int_t nr_bl = magma_ceildiv(NE,10000);        //nr of blocks
                   magma_int_t sz_bl = magma_ceildiv(NE,nr_bl*64)*64; //maximum size of blocks (to have blocks of around the same size and multiple of 64)
                   magma_int_t ib;                                      //size of current block

                   for(magma_int_t i=0; i<NE; i+= sz_bl){
                       ib = min(sz_bl, NE-i);
                       magma_zlarfb_gpu( MagmaLeft, MagmaNoTrans, MagmaForward, MagmaColumnwise, vlen, ib, vnb, dV, lddv, dT, lddt, dE(fst,i), ldde, dwork, NE);
                   }
               }
           }
        }
    }else if (side=='R'){
         rownbm    = magma_ceildiv((N-1),NB);
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
                   magma_zlarfb_gpu( MagmaRight, MagmaNoTrans, MagmaForward, MagmaColumnwise, NE, vlen, vnb, dV, lddv, dT, lddt, dE(0, fst), ldde, dwork, NE);

               }
            }
         }
    }else{
            printf("ERROR SIDE %d \n",side);
    }

}
#undef E
#undef V
#undef TAU
#undef T
#undef dE
#undef dV
#undef dT
////////////////////////////////////////////////////////////////////////////////////////////////////



