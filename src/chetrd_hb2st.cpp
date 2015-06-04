/*
 -- MAGMA (version 1.3.0) --
 Univ. of Tennessee, Knoxville
 Univ. of California, Berkeley
 Univ. of Colorado, Denver
 November 2012

 @author Azzam Haidar
 @author Stan Tomov

 @generated c Wed Nov 14 22:53:26 2012

 */
#include "common_magma.h"
#include "magma_bulge.h"
#include <cblas.h>

#if defined(USEMKL)
#include <mkl_service.h>
#endif
#if defined(USEACML)
#include <omp.h>
#endif
// === Define what BLAS to use ============================================
#define PRECISION_c
#if (defined(PRECISION_s) || defined(PRECISION_d))
//#define magma_cgemm magmablas_cgemm
#endif
// === End defining what BLAS to use =======================================
extern "C" {
    
    void magma_ctrdtype1cbHLsym_withQ_v2(magma_int_t n, magma_int_t nb, cuFloatComplex *A, magma_int_t lda, cuFloatComplex *V, magma_int_t ldv, cuFloatComplex *TAU,
                                         magma_int_t st, magma_int_t ed, magma_int_t sweep, magma_int_t Vblksiz, cuFloatComplex *work);
    void magma_ctrdtype2cbHLsym_withQ_v2(magma_int_t n, magma_int_t nb, cuFloatComplex *A, magma_int_t lda, cuFloatComplex *V, magma_int_t ldv, cuFloatComplex *TAU,
                                         magma_int_t st, magma_int_t ed, magma_int_t sweep, magma_int_t Vblksiz, cuFloatComplex *work);
    void magma_ctrdtype3cbHLsym_withQ_v2(magma_int_t n, magma_int_t nb, cuFloatComplex *A, magma_int_t LDA, cuFloatComplex *V, magma_int_t ldv, cuFloatComplex *TAU,
                                         magma_int_t st, magma_int_t ed, magma_int_t sweep, magma_int_t Vblksiz, cuFloatComplex *work);
        
}

static void *magma_chetrd_hb2st_parallel_section(void *arg);

static void magma_ctile_bulge_parallel(magma_int_t my_core_id, magma_int_t cores_num, cuFloatComplex *A, magma_int_t lda,
                                       cuFloatComplex *V, magma_int_t ldv, cuFloatComplex *TAU, magma_int_t n, magma_int_t nb, magma_int_t nbtiles,
                                       magma_int_t grsiz, magma_int_t Vblksiz, volatile magma_int_t *prog);

static void magma_ctile_bulge_computeT_parallel(magma_int_t my_core_id, magma_int_t cores_num, cuFloatComplex *V, magma_int_t ldv, cuFloatComplex *TAU,
                                                cuFloatComplex *T, magma_int_t ldt, magma_int_t n, magma_int_t nb, magma_int_t Vblksiz);

//////////////////////////////////////////////////////////////////////////////////////////////////////////

class magma_cbulge_data {

public:

    magma_cbulge_data(magma_int_t threads_num_, magma_int_t n_, magma_int_t nb_, magma_int_t nbtiles_,
                      magma_int_t grsiz_, magma_int_t Vblksiz_, magma_int_t compT_,
                      cuFloatComplex *A_, magma_int_t lda_, cuFloatComplex *V_, magma_int_t ldv_, cuFloatComplex *TAU_,
                      cuFloatComplex *T_, magma_int_t ldt_, volatile magma_int_t* prog_)
    :
    threads_num(threads_num_),
    n(n_),
    nb(nb_),
    nbtiles(nbtiles_),
    grsiz(grsiz_),
    Vblksiz(Vblksiz_),
    compT(compT_),
    A(A_),
    lda(lda_),
    V(V_),
    ldv(ldv_),
    TAU(TAU_),
    T(T_),
    ldt(ldt_),
    prog(prog_)
    {
        pthread_barrier_init(&barrier, NULL, threads_num);
    }

    ~magma_cbulge_data()
    {
        pthread_barrier_destroy(&barrier);
    }

    const magma_int_t threads_num;
    const magma_int_t n;
    const magma_int_t nb;
    const magma_int_t nbtiles;
    const magma_int_t grsiz;
    const magma_int_t Vblksiz;
    const magma_int_t compT;
    cuFloatComplex* const A;
    const magma_int_t lda;
    cuFloatComplex* const V;
    const magma_int_t ldv;
    cuFloatComplex* const TAU;
    cuFloatComplex* const T;
    const magma_int_t ldt;
    volatile magma_int_t *prog;
    pthread_barrier_t barrier;

private:

    magma_cbulge_data(magma_cbulge_data& data); // disable copy

};

class magma_cbulge_id_data {
    
public:
    
    magma_cbulge_id_data()
    : id(-1), data(NULL)
    {}
    
    magma_cbulge_id_data(magma_int_t id_, magma_cbulge_data* data_)
    : id(id_), data(data_)
    {}
    
    magma_int_t id;
    magma_cbulge_data* data;
    
};


///////////////////////////////////////////////////////////////////////////////////////////////////////////

extern "C" magma_int_t magma_chetrd_hb2st(magma_int_t threads, char uplo, magma_int_t n, magma_int_t nb, magma_int_t Vblksiz,
                                          cuFloatComplex *A, magma_int_t lda, float *D, float *E,
                                          cuFloatComplex *V, magma_int_t ldv, cuFloatComplex *TAU, magma_int_t compT, cuFloatComplex *T, magma_int_t ldt)
{
    /*  -- MAGMA (version 1.3.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2012

    Purpose
    =======


    Arguments
    =========
    THREADS (input) INTEGER
            Specifies the number of pthreads used.
            THREADS > 0

    UPLO    (input) CHARACTER*1
            = 'U':  Upper triangles of A is stored;
            = 'L':  Lower triangles of A is stored.

    N       (input) INTEGER
            The order of the matrix A.  N >= 0.
     
    NB      (input) INTEGER
            The order of the band matrix A.  N >= NB >= 0.
     
    VBLKSIZ (input) INTEGER
            The size of the block of householder vectors applied at once.

    A       (input/workspace) COMPLEX*16 array, dimension (LDA, N)
            On entry the band matrix stored in the following way:

    LDA     (input) INTEGER
            The leading dimension of the array A.  LDA >= 2*NB.

    D       (output) DOUBLE array, dimension (N)   
            The diagonal elements of the tridiagonal matrix T:   
            D(i) = A(i,i).   

    E       (output) DOUBLE array, dimension (N-1)   
            The off-diagonal elements of the tridiagonal matrix T:   
            E(i) = A(i,i+1) if UPLO = 'U', E(i) = A(i+1,i) if UPLO = 'L'.

    V       (output) COMPLEX*16 array, dimension (BLKCNT, LDV, VBLKSIZ)
            On exit it contains the blocks of householder reflectors
            BLKCNT is the number of block and it is returned by the funtion MAGMA_BULGE_GET_BLKCNT.

    LDV     (input) INTEGER
            The leading dimension of V.
            LDV > NB + VBLKSIZ + 1

    TAU     (output) COMPLEX*16 dimension(BLKCNT, VBLKSIZ)
            ???
     
    COMPT   (input) INTEGER
            if COMPT = 0 T is not computed
            if COMPT = 1 T is computed

    T       (output) COMPLEX*16 dimension(LDT *)
            if COMPT = 1 on exit contains the matrices T needed for Q2
            if COMPT = 0 T is not referenced
     
    LDT     (input) INTEGER
            The leading dimension of T.
            LDT > Vblksiz
     
    INFO    (output) INTEGER ????????????????????????????????????????????????????????????????????????????????????
            = 0:  successful exit
            

    =====================================================================  */
    
    char uplo_[2] = {uplo, 0};
    float timeblg=0.0;

    magma_int_t mklth = threads;

    magma_int_t INgrsiz=1;

    magma_int_t blkcnt = magma_bulge_get_blkcnt(n, nb, Vblksiz);

    magma_int_t nbtiles = magma_ceildiv(n, nb);
    
    memset(T,   0, blkcnt*ldt*Vblksiz*sizeof(cuFloatComplex));
    memset(TAU, 0, blkcnt*Vblksiz*sizeof(cuFloatComplex));
    memset(V,   0, blkcnt*ldv*Vblksiz*sizeof(cuFloatComplex));
    
    magma_int_t* prog = new magma_int_t[2*nbtiles+threads+10];
    memset(prog, 0, (2*nbtiles+threads+10)*sizeof(magma_int_t));
    
    magma_cbulge_id_data* arg = new magma_cbulge_id_data[threads];
    pthread_t* thread_id = new pthread_t[threads];

    pthread_attr_t thread_attr;
    
#if defined(USEMKL)
    mkl_set_num_threads( 1 );
#endif
#if defined(USEACML)
    omp_set_num_threads(1);
#endif

    magma_cbulge_data data_bulge(threads, n, nb, nbtiles, INgrsiz, Vblksiz, compT,
                                 A, lda, V, ldv, TAU, T, ldt, prog);

    // Set one thread per core
    pthread_attr_init(&thread_attr);
    pthread_attr_setscope(&thread_attr, PTHREAD_SCOPE_SYSTEM);
    pthread_setconcurrency(threads);

    //timing
    timeblg = magma_wtime();

    // Launch threads
    for (magma_int_t thread = 1; thread < threads; thread++)
    {
        arg[thread] = magma_cbulge_id_data(thread, &data_bulge);
        pthread_create(&thread_id[thread], &thread_attr, magma_chetrd_hb2st_parallel_section, &arg[thread]);
    }
    arg[0] = magma_cbulge_id_data(0, &data_bulge);
    magma_chetrd_hb2st_parallel_section(&arg[0]);

    // Wait for completion
    for (magma_int_t thread = 1; thread < threads; thread++)
    {
        void *exitcodep;
        pthread_join(thread_id[thread], &exitcodep);
    }
    
    // timing
    timeblg = magma_wtime()-timeblg;

    
    delete[] thread_id;
    delete[] arg;
    delete[] prog;
    
    printf("time BULGE+T = %f \n" ,timeblg);

#if defined(USEMKL)
    mkl_set_num_threads( mklth );
#endif
#if defined(USEACML)
    omp_set_num_threads(mklth);
#endif
    
    /*================================================
     *  store resulting diag and lower diag D and E
     *  note that D and E are always real
     *================================================*/

    /* Make diagonal and superdiagonal elements real,
     * storing them in D and E
     */
    /* In complex case, the off diagonal element are
     * not necessary real. we have to make off-diagonal
     * elements real and copy them to E.
     * When using HouseHolder elimination,
     * the CLARFG give us a real as output so, all the
     * diagonal/off-diagonal element except the last one are already
     * real and thus we need only to take the abs of the last
     * one.
     *  */

#if defined(PRECISION_z) || defined(PRECISION_c)
    if(uplo==MagmaLower){
        for (magma_int_t i=0; i < n-1 ; i++)
        {
            D[i] = MAGMA_C_REAL(A[i*lda  ]);
            E[i] = MAGMA_C_REAL(A[i*lda+1]);
        }
        D[n-1] = MAGMA_C_REAL(A[(n-1)*lda]);
    } else { /* MagmaUpper not tested yet */
        for (magma_int_t i=0; i<n-1; i++)
        {
            D[i]  =  MAGMA_C_REAL(A[i*lda+nb]);
            E[i] = MAGMA_C_REAL(A[i*lda+nb-1]);
        }
        D[n-1] = MAGMA_C_REAL(A[(n-1)*lda+nb]);
    } /* end MagmaUpper */
#else
    if( uplo == MagmaLower ){
        for (magma_int_t i=0; i < n-1; i++) {
            D[i] = A[i*lda];   // diag
            E[i] = A[i*lda+1]; //lower diag
        }
        D[n-1] = A[(n-1)*lda];
    } else {
        for (magma_int_t i=0; i < n-1; i++) {
            D[i] = A[i*lda+nb];   // diag
            E[i] = A[i*lda+nb-1]; //lower diag
        }
        D[n-1] = A[(n-1)*lda+nb];
    }
#endif
    return MAGMA_SUCCESS;
    
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////

static void *magma_chetrd_hb2st_parallel_section(void *arg)
{
    magma_int_t my_core_id  = ((magma_cbulge_id_data*)arg) -> id;
    magma_cbulge_data* data = ((magma_cbulge_id_data*)arg) -> data;

    magma_int_t allcores_num   = data -> threads_num;
    magma_int_t n              = data -> n;
    magma_int_t nb             = data -> nb;
    magma_int_t nbtiles        = data -> nbtiles;
    magma_int_t grsiz          = data -> grsiz;
    magma_int_t Vblksiz        = data -> Vblksiz;
    magma_int_t compT          = data -> compT;
    cuFloatComplex *A         = data -> A;
    magma_int_t lda            = data -> lda;
    cuFloatComplex *V         = data -> V;
    magma_int_t ldv            = data -> ldv;
    cuFloatComplex *TAU       = data -> TAU;
    cuFloatComplex *T         = data -> T;
    magma_int_t ldt            = data -> ldt;
    volatile magma_int_t* prog = data -> prog;

    pthread_barrier_t* barrier = &(data -> barrier);

    magma_int_t sys_corenbr    = 1;

    float timeB=0.0, timeT=0.0;

#if defined(SETAFFINITY)
    // bind threads
    cpu_set_t set;
    // bind threads
    CPU_ZERO( &set );
    CPU_SET( my_core_id, &set );
    sched_setaffinity( 0, sizeof(set), &set) ;
#endif

    if(compT==1)
    {
        /* compute the Q1 overlapped with the bulge chasing+T.
         * if all_cores_num=1 it call Q1 on GPU and then bulgechasing.
         * otherwise the first thread run Q1 on GPU and
         * the other threads run the bulgechasing.
         * */

        if(allcores_num==1)
        {

            //=========================
            //    bulge chasing
            //=========================
            timeB = magma_wtime();

            magma_ctile_bulge_parallel(0, 1, A, lda, V, ldv, TAU, n, nb, nbtiles, grsiz, Vblksiz, prog);

            timeB = magma_wtime()-timeB;
            printf("  Finish BULGE   timing= %f \n" ,timeB);


            //=========================
            // compute the T's to be used when applying Q2
            //=========================
            timeT = magma_wtime();
            magma_ctile_bulge_computeT_parallel(0, 1, V, ldv, TAU, T, ldt, n, nb, Vblksiz);

            timeT = magma_wtime()-timeT;
            printf("  Finish T's     timing= %f \n" ,timeT);

        }else{ // allcore_num > 1

            magma_int_t id  = my_core_id;
            magma_int_t tot = allcores_num;


                //=========================
                //    bulge chasing
                //=========================
                if(id == 0)timeB = magma_wtime();

                magma_ctile_bulge_parallel(id, tot, A, lda, V, ldv, TAU, n, nb, nbtiles, grsiz, Vblksiz, prog);
                pthread_barrier_wait(barrier);

                if(id == 0){
                    timeB = magma_wtime()-timeB;
                    printf("  Finish BULGE   timing= %f \n" ,timeB);
                }

                //=========================
                // compute the T's to be used when applying Q2
                //=========================
                if(id == 0)timeT = magma_wtime();

                magma_ctile_bulge_computeT_parallel(id, tot, V, ldv, TAU, T, ldt, n, nb, Vblksiz);
                pthread_barrier_wait(barrier);

                if (id == 0){
                    timeT = magma_wtime()-timeT;
                    printf("  Finish T's     timing= %f \n" ,timeT);
                }

        } // allcore == 1

    }else{ // WANTZ = 0

        //=========================
        //    bulge chasing
        //=========================
        if(my_core_id == 0)
            timeB = magma_wtime();
        
        magma_ctile_bulge_parallel(my_core_id, allcores_num, A, lda, V, ldv, TAU, n, nb, nbtiles, grsiz, Vblksiz, prog);

        pthread_barrier_wait(barrier);
        
        if(my_core_id == 0){
            timeB = magma_wtime()-timeB;
            printf("  Finish BULGE   timing= %f \n" ,timeB);
        }
    } // WANTZ > 0

#if defined(SETAFFINITY)
    // unbind threads
    sys_corenbr = sysconf(_SC_NPROCESSORS_ONLN);
    CPU_ZERO( &set );
    for(magma_int_t i=0; i<sys_corenbr; i++)
        CPU_SET( i, &set );
    sched_setaffinity( 0, sizeof(set), &set) ;
#endif

    return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

static void magma_ctile_bulge_parallel(magma_int_t my_core_id, magma_int_t cores_num, cuFloatComplex *A, magma_int_t lda,
                                       cuFloatComplex *V, magma_int_t ldv, cuFloatComplex *TAU, magma_int_t n, magma_int_t nb, magma_int_t nbtiles,
                                       magma_int_t grsiz, magma_int_t Vblksiz, volatile magma_int_t *prog)
{
    magma_int_t sweepid, myid, shift, stt, st, ed, stind, edind;
    magma_int_t blklastind, colpt;
    magma_int_t stepercol;
    magma_int_t i,j,m,k;
    magma_int_t thgrsiz, thgrnb, thgrid, thed;
    magma_int_t coreid;
    magma_int_t colblktile,maxrequiredcores,colpercore,mycoresnb;
    magma_int_t fin;
    cuFloatComplex *work;
    
    if(n<=0)
        return ;

    //printf("=================> my core id %d of %d \n",my_core_id, cores_num);

    /* As I store V in the V vector there are overlap between
     * tasks so shift is now 4 where group need to be always
     * multiple of 2, because as example if grs=1 task 2 from
     * sweep 2 can run with task 6 sweep 1., but task 2 sweep 2
     * will overwrite the V of tasks 5 sweep 1 which are used by
     * task 6, so keep in mind that group need to be multiple of 2,
     * and thus tasks 2 sweep 2 will never run with task 6 sweep 1.
     * However, when storing V in A, shift could be back to 3.
     * */

    magma_cmalloc_cpu(&work, n);
    
    mycoresnb = cores_num;

    shift   = 5;
    if(grsiz==1)
        colblktile=1;
    else
        colblktile=grsiz/2;

    maxrequiredcores = nbtiles/colblktile;
    if(maxrequiredcores<1)maxrequiredcores=1;
    colpercore  = colblktile*nb;
    if(mycoresnb > maxrequiredcores)
    {
        if(my_core_id==0)printf("==================================================================================\n");
        if(my_core_id==0)printf("  WARNING only %3d threads are required to run this test optimizing cache reuse\n",maxrequiredcores);
        if(my_core_id==0)printf("==================================================================================\n");
        mycoresnb = maxrequiredcores;
    }
    thgrsiz = n;
    
    if(my_core_id==0) printf("  Static bulgechasing version v9_9col threads  %4d      N %5d      NB %5d    grs %4d thgrsiz %4d \n",cores_num, n, nb, grsiz,thgrsiz);

    stepercol = magma_ceildiv(shift, grsiz);

    thgrnb  = magma_ceildiv(n-1, thgrsiz);
    
    for (thgrid = 1; thgrid<=thgrnb; thgrid++){
        stt  = (thgrid-1)*thgrsiz+1;
        thed = min( (stt + thgrsiz -1), (n-1));
        for (i = stt; i <= n-1; i++){
            ed=min(i,thed);
            if(stt>ed)break;
            for (m = 1; m <=stepercol; m++){
                st=stt;
                for (sweepid = st; sweepid <=ed; sweepid++){

                    for (k = 1; k <=grsiz; k++){
                        myid = (i-sweepid)*(stepercol*grsiz) +(m-1)*grsiz + k;
                        if(myid%2 ==0){
                            colpt      = (myid/2)*nb+1+sweepid-1;
                            stind      = colpt-nb+1;
                            edind      = min(colpt,n);
                            blklastind = colpt;
                            if(stind>=edind){
                                printf("ERROR---------> st>=ed  %d  %d \n\n",stind, edind);
                                exit(-10);
                            }
                        }else{
                            colpt      = ((myid+1)/2)*nb + 1 +sweepid -1 ;
                            stind      = colpt-nb+1;
                            edind      = min(colpt,n);
                            if( (stind>=edind-1) && (edind==n) )
                                blklastind=n;
                            else
                                blklastind=0;
                            if(stind>edind){
                                printf("ERROR---------> st>=ed  %d  %d \n\n",stind, edind);
                                exit(-10);
                            }
                        }

                        coreid = (stind/colpercore)%mycoresnb;

                        if(my_core_id==coreid)
                        {

                            fin=0;
                            while(fin==0)
                            {
                                if(myid==1)
                                {
                                    if( (prog[myid+shift-1]== (sweepid-1)) )
                                    {
                                        magma_ctrdtype1cbHLsym_withQ_v2(n, nb, A, lda, V, ldv, TAU, stind, edind, sweepid, Vblksiz, work);                                        

                                        fin=1;
                                        prog[myid]= sweepid;
                                        if(blklastind >= (n-1))
                                        {
                                            for (j = 1; j <= shift; j++)
                                                prog[myid+j]=sweepid;
                                        }
                                    } // END progress condition

                                }else{
                                    if( (prog[myid-1]==sweepid) && (prog[myid+shift-1]== (sweepid-1)) )
                                    {
                                        if(myid%2 == 0)
                                            magma_ctrdtype2cbHLsym_withQ_v2(n, nb, A, lda, V, ldv, TAU, stind, edind, sweepid, Vblksiz, work);
                                        else
                                            magma_ctrdtype3cbHLsym_withQ_v2(n, nb, A, lda, V, ldv, TAU, stind, edind, sweepid, Vblksiz, work);                                      

                                        fin=1;
                                        prog[myid]= sweepid;
                                        if(blklastind >= (n-1))
                                        {
                                            for (j = 1; j <= shift+mycoresnb; j++)
                                                prog[myid+j]=sweepid;
                                        }
                                    } // END progress condition
                                } // END if myid==1
                            } // END while loop

                        } // END if my_core_id==coreid

                        if(blklastind >= (n-1))
                        {
                            stt=stt+1;
                            break;
                        }
                    }   // END for k=1:grsiz
                } // END for sweepid=st:ed
            } // END for m=1:stepercol
        } // END for i=1:n-1
    } // END for thgrid=1:thgrnb
    
    magma_free_cpu(work);

} // END FUNCTION
////////////////////////////////////////////////////////////////////////////////////////////////////

#define V(m)     &(V[(m)])
#define TAU(m)   &(TAU[(m)])
#define T(m)   &(T[(m)])
static void magma_ctile_bulge_computeT_parallel(magma_int_t my_core_id, magma_int_t cores_num, cuFloatComplex *V, magma_int_t ldv, cuFloatComplex *TAU,
                                                cuFloatComplex *T, magma_int_t ldt, magma_int_t n, magma_int_t nb, magma_int_t Vblksiz)
{
    //%===========================
    //%   local variables
    //%===========================
    magma_int_t firstcolj;
    magma_int_t rownbm;
    magma_int_t st,ed,fst,vlen,vnb,colj;
    magma_int_t blkid,vpos,taupos,tpos;
    magma_int_t blkpercore, myid;
    
    if(n<=0)
        return ;
    
    magma_int_t blkcnt = magma_bulge_get_blkcnt(n, nb, Vblksiz);
    
    blkpercore = blkcnt/cores_num;
    
    magma_int_t nbGblk  = magma_ceildiv(n-1, Vblksiz);
    
    if(my_core_id==0) printf("  COMPUTE T parallel threads %d with  N %d   NB %d   Vblksiz %d \n",cores_num,n,nb,Vblksiz);
    
    for (magma_int_t bg = nbGblk; bg>0; bg--)
    {
        firstcolj = (bg-1)*Vblksiz + 1;
        rownbm    = magma_ceildiv(n-(firstcolj+1), nb);
        if(bg==nbGblk) 
            rownbm    = magma_ceildiv(n-firstcolj ,nb);  // last blk has size=1 used for complex to handle A(N,N-1)

        for (magma_int_t m = rownbm; m>0; m--)
        {
            vlen = 0;
            vnb  = 0;
            colj      = (bg-1)*Vblksiz; // for k=0;I compute the fst and then can remove it from the loop
            fst       = (rownbm -m)*nb+colj +1;
            for (magma_int_t k=0; k<Vblksiz; k++)
            {
                colj     = (bg-1)*Vblksiz + k;
                st       = (rownbm -m)*nb+colj +1;
                ed       = min(st+nb-1,n-1);
                if(st>ed)
                    break;
                if((st==ed)&&(colj!=n-2))
                    break;
                
                vlen=ed-fst+1;
                vnb=k+1;
            }        
            colj     = (bg-1)*Vblksiz;
            magma_bulge_findVTAUTpos(n, nb, Vblksiz, colj, fst, ldv, ldt, &vpos, &taupos, &tpos, &blkid);
            myid = blkid/blkpercore;
            if(my_core_id==(myid%cores_num)){
                if((vlen>0)&&(vnb>0))
                    lapackf77_clarft( "F", "C", &vlen, &vnb, V(vpos), &ldv, TAU(taupos), T(tpos), &ldt);
            }
        }
    }
}
#undef V
#undef TAU
#undef T
////////////////////////////////////////////////////////////////////////////////////////////////////
