/*
 -- MAGMA (version 1.3.0) --
 Univ. of Tennessee, Knoxville
 Univ. of California, Berkeley
 Univ. of Colorado, Denver
 November 2012

 @author Azzam Haidar
 @author Stan Tomov

 @precisions normal z -> s d c

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
#define PRECISION_z
#if (defined(PRECISION_s) || defined(PRECISION_d))
//#define magma_zgemm magmablas_zgemm
#endif
// === End defining what BLAS to use =======================================
extern "C" {
    void magma_zbulge_applyQ_v2(char side, magma_int_t NE, magma_int_t N, magma_int_t NB, magma_int_t Vblksiz, cuDoubleComplex *dE,
                                magma_int_t ldde, cuDoubleComplex *V, magma_int_t ldv, cuDoubleComplex *T, magma_int_t ldt, magma_int_t *info);

}

static void *magma_zapplyQ_parallel_section(void *arg);

static void magma_ztile_bulge_applyQ(char side, magma_int_t n_loc, magma_int_t n, magma_int_t nb, magma_int_t Vblksiz,
                                     cuDoubleComplex *E, magma_int_t lde, cuDoubleComplex *V, magma_int_t ldv,
                                     cuDoubleComplex *TAU, cuDoubleComplex *T, magma_int_t ldt);

//////////////////////////////////////////////////////////////////////////////////////////////////////////

class magma_zapplyQ_data {
    
public:
    
    magma_zapplyQ_data(magma_int_t threads_num_, magma_int_t n_, magma_int_t ne_, magma_int_t n_gpu_,
                magma_int_t nb_, magma_int_t Vblksiz_, cuDoubleComplex *E_, magma_int_t lde_,
                cuDoubleComplex *V_, magma_int_t ldv_, cuDoubleComplex *TAU_,
                cuDoubleComplex *T_, magma_int_t ldt_, cuDoubleComplex *dE_, magma_int_t ldde_)
    :
    threads_num(threads_num_),
    n(n_),
    ne(ne_),
    n_gpu(n_gpu_),
    nb(nb_),
    Vblksiz(Vblksiz_),
    E(E_),
    lde(lde_),
    V(V_),
    ldv(ldv_),
    TAU(TAU_),
    T(T_),
    ldt(ldt_),
    dE(dE_),
    ldde(ldde_)
    {
        magma_int_t count = threads_num;
        
        if(threads_num > 1)
            --count;
        
        pthread_barrier_init(&barrier, NULL, count);
    }
    
    ~magma_zapplyQ_data()
    {
        pthread_barrier_destroy(&barrier);
    }
    
    const magma_int_t threads_num;
    const magma_int_t n;
    const magma_int_t ne;
    const magma_int_t n_gpu;
    const magma_int_t nb;
    const magma_int_t Vblksiz;
    cuDoubleComplex* const E;
    const magma_int_t lde;
    cuDoubleComplex* const V;
    const magma_int_t ldv;
    cuDoubleComplex* const TAU;
    cuDoubleComplex* const T;
    const magma_int_t ldt;
    cuDoubleComplex* const dE;
    const magma_int_t ldde;
    pthread_barrier_t barrier;
    
private:
    
    magma_zapplyQ_data(magma_zapplyQ_data& data); // disable copy
    
};

class magma_zapplyQ_id_data {
    
public:
    
    magma_zapplyQ_id_data()
    : id(-1), data(NULL)
    {}
    
    magma_zapplyQ_id_data(magma_int_t id_, magma_zapplyQ_data* data_)
    : id(id_), data(data_)
    {}
    
    magma_int_t id;
    magma_zapplyQ_data* data;
    
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////

extern "C" magma_int_t magma_zbulge_back(magma_int_t threads, char uplo, magma_int_t n, magma_int_t nb, magma_int_t ne, magma_int_t Vblksiz,
                                         cuDoubleComplex *Z, magma_int_t ldz, cuDoubleComplex *dZ, magma_int_t lddz,
                                         cuDoubleComplex *V, magma_int_t ldv, cuDoubleComplex *TAU, cuDoubleComplex *T, magma_int_t ldt, magma_int_t* info)
{
    magma_int_t mklth = threads;
    
    double timeaplQ2=0.0;
    
#if defined(USEMKL)
        mkl_set_num_threads(1);
#endif
#if defined(USEACML)
        omp_set_num_threads(1);
#endif
    
            double f= 1.;
            magma_int_t n_gpu = ne;
            
            if(threads>40){
                f = 0.5;
                n_gpu = (magma_int_t)(f*ne)/64*64;
            }
            else if(threads>10){
#if (defined(PRECISION_s) || defined(PRECISION_d))
                f = 0.68;
#else
                f = 0.72;
#endif
                n_gpu = (magma_int_t)(f*ne)/64*64;
            }
            else if(threads>5){
#if (defined(PRECISION_s) || defined(PRECISION_d))
                f = 0.82;
#else
                f = 0.86;
#endif
                n_gpu = (magma_int_t)(f*ne)/64*64;
            }            
            else if(threads>1){
#if (defined(PRECISION_s) || defined(PRECISION_d))
                f = 0.96;
#else
                f = 0.96;
#endif
                n_gpu = (magma_int_t)(f*ne)/64*64;
            }
            
            /****************************************************
             *  apply V2 from left to the eigenvectors Z. dZ = (I-V2*T2*V2')*Z
             * **************************************************/

            timeaplQ2 = magma_wtime();
            
            /*============================
             *  use GPU+CPU's
             *==========================*/  
            
            if(n_gpu < ne)
            {
                
                // define the size of Q to be done on CPU's and the size on GPU's
                // note that GPU use Q(1:N_GPU) and CPU use Q(N_GPU+1:N)

                printf("---> calling GPU + CPU(if N_CPU>0) to apply V2 to Z with NE %d     N_GPU %d   N_CPU %d\n",ne, n_gpu, ne-n_gpu); 
                
                magma_zapplyQ_data data_applyQ(threads, n, ne, n_gpu, nb, Vblksiz, Z, ldz, V, ldv, TAU, T, ldt, dZ, lddz);
                
                magma_zapplyQ_id_data* arg = new magma_zapplyQ_id_data[threads];
                pthread_t* thread_id = new pthread_t[threads];
                
                pthread_attr_t thread_attr;
                
                // ===============================
                // relaunch thread to apply Q
                // ===============================
                // Set one thread per core
                pthread_attr_init(&thread_attr);
                pthread_attr_setscope(&thread_attr, PTHREAD_SCOPE_SYSTEM);
                pthread_setconcurrency(threads);
                
                // Launch threads
                for (magma_int_t thread = 1; thread < threads; thread++)
                {
                    arg[thread] = magma_zapplyQ_id_data(thread, &data_applyQ);
                    pthread_create(&thread_id[thread], &thread_attr, magma_zapplyQ_parallel_section, &arg[thread]);
                }
                arg[0] = magma_zapplyQ_id_data(0, &data_applyQ);
                magma_zapplyQ_parallel_section(&arg[0]);
                
                // Wait for completion
                for (magma_int_t thread = 1; thread < threads; thread++)
                {
                    void *exitcodep;
                    pthread_join(thread_id[thread], &exitcodep);
                }
                
                delete[] thread_id;
                delete[] arg;
                
                magma_zsetmatrix(n, ne-n_gpu, Z + n_gpu*ldz, ldz, dZ + n_gpu*ldz, lddz);
                
                /*============================
                 *  use only GPU
                 *==========================*/  
            }else{
                magma_zsetmatrix(n, ne, Z, ldz, dZ, lddz);
                magma_zbulge_applyQ_v2('L', ne, n, nb, Vblksiz, dZ, lddz, V, ldv, T, ldt, info);
                magma_device_sync();
            }

            timeaplQ2 = magma_wtime()-timeaplQ2;
            
#if defined(USEMKL)
        mkl_set_num_threads(mklth);
#endif
#if defined(USEACML)
        omp_set_num_threads(mklth);
#endif
    
    return MAGMA_SUCCESS;
}

//##################################################################################################
static void *magma_zapplyQ_parallel_section(void *arg)
{

    magma_int_t my_core_id   = ((magma_zapplyQ_id_data*)arg) -> id;
    magma_zapplyQ_data* data = ((magma_zapplyQ_id_data*)arg) -> data;
    
    magma_int_t allcores_num   = data -> threads_num;
    magma_int_t n              = data -> n;
    magma_int_t ne             = data -> ne;
    magma_int_t n_gpu          = data -> n_gpu;
    magma_int_t nb             = data -> nb;
    magma_int_t Vblksiz        = data -> Vblksiz;
    cuDoubleComplex *E         = data -> E;
    magma_int_t lde            = data -> lde;
    cuDoubleComplex *V         = data -> V;
    magma_int_t ldv            = data -> ldv;
    cuDoubleComplex *TAU       = data -> TAU;
    cuDoubleComplex *T         = data -> T;
    magma_int_t ldt            = data -> ldt;
    cuDoubleComplex *dE        = data -> dE;
    magma_int_t ldde           = data -> ldde;
    pthread_barrier_t* barrier = &(data -> barrier);
 
    magma_int_t info;
    
    real_Double_t timeQcpu=0.0, timeQgpu=0.0;
    
    magma_int_t n_cpu = ne - n_gpu;

#if defined(SETAFFINITY)    
    cpu_set_t set;
    CPU_ZERO( &set );
    CPU_SET( my_core_id, &set );
    sched_setaffinity( 0, sizeof(set), &set) ;
#endif
   
            if(my_core_id==0)
            {
                //=============================================
                //   on GPU on thread 0:
                //    - apply V2*Z(:,1:N_GPU)
                //=============================================
                timeQgpu = magma_wtime();
                
                magma_zsetmatrix(n, n_gpu, E, lde, dE, ldde);
                magma_zbulge_applyQ_v2('L', n_gpu, n, nb, Vblksiz, dE, ldde, V, ldv, T, ldt, &info);
                
                magma_device_sync();
                timeQgpu = magma_wtime()-timeQgpu;
                printf("  Finish Q2_GPU GGG timing= %lf \n" ,timeQgpu);

            }else{
                //=============================================
                //   on CPU on threads 1:allcores_num-1:
                //    - apply V2*Z(:,N_GPU+1:NE)
                //=============================================
                if(my_core_id == 1)
                    timeQcpu = magma_wtime();
                
                magma_int_t n_loc = magma_ceildiv(n_cpu, allcores_num-1);
                cuDoubleComplex* E_loc = E + (n_gpu+ n_loc * (my_core_id-1))*lde;
                n_loc = min(n_loc,n_cpu - n_loc * (my_core_id-1));
                
                magma_ztile_bulge_applyQ('L', n_loc, n, nb, Vblksiz, E_loc, lde, V, ldv, TAU, T, ldt);
                pthread_barrier_wait(barrier);
                if(my_core_id == 1){
                    timeQcpu = magma_wtime()-timeQcpu;
                    printf("  Finish Q2_CPU CCC timing= %lf \n" ,timeQcpu);
                }
                
            } // END if my_core_id
        
    
#if defined(SETAFFINITY)    
    // unbind threads 
    magma_int_t sys_corenbr = 1;
    sys_corenbr = sysconf(_SC_NPROCESSORS_ONLN);
    CPU_ZERO( &set );
    for(magma_int_t i=0; i<sys_corenbr; i++)
        CPU_SET( i, &set );
    sched_setaffinity( 0, sizeof(set), &set) ;
#endif
    
    return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

#define E(m,n)   &(E[(m) + lde*(n)])
#define V(m)     &(V[(m)])
#define TAU(m)   &(TAU[(m)])
#define T(m)     &(T[(m)])
static void magma_ztile_bulge_applyQ(char side, magma_int_t n_loc, magma_int_t n, magma_int_t nb, magma_int_t Vblksiz,
                                     cuDoubleComplex *E, magma_int_t lde, cuDoubleComplex *V, magma_int_t ldv,
                                     cuDoubleComplex *TAU, cuDoubleComplex *T, magma_int_t ldt)//, magma_int_t* info)
{
    //%===========================
    //%   local variables
    //%===========================
    magma_int_t firstcolj;
    magma_int_t bg, rownbm;
    magma_int_t st,ed,fst,vlen,vnb,colj;
    magma_int_t vpos,tpos;
    magma_int_t cur_blksiz,avai_blksiz, ncolinvolvd;
    magma_int_t nbgr, colst, coled;
    
    if(n<=0)
        return ;
    if(n_loc<=0)
        return ;
    
    //info = 0;
    magma_int_t INFO=0;    
    
    magma_int_t nbGblk  = magma_ceildiv(n-1, Vblksiz);
    
    /*
     * version v1: for each chunck it apply all the V's then move to 
     *                    the other chunck. the locality here inside each 
     *                    chunck meaning that thread t apply V_k then move 
     *                    to V_k+1 which overlap with V_k meaning that the 
     *                    E_k+1 overlap with E_k. so here is the 
     *                    locality however thread t had to read V_k+1 and 
     *                    T_k+1 at each apply. note that all thread if they 
     *                    run at same speed they might reading the same V_k 
     *                    and T_k at the same time.
     * */
    
    magma_int_t nb_loc = 128; //$$$$$$$$
    
    magma_int_t     lwork = 2*nb_loc*max(Vblksiz,64);
    cuDoubleComplex *work, *work2;
    
    magma_zmalloc_cpu(&work, lwork);
    magma_zmalloc_cpu(&work2, lwork);
    
    magma_int_t nbchunk =  magma_ceildiv(n_loc, nb_loc);
    
    /* SIDE LEFT  meaning apply E = Q*E = (q_1*q_2*.....*q_n) * E ==> so traverse Vs in reverse order (forward) from q_n to q_1
     *            each q_i consist of applying V to a block of row E(row_i,:) and applies are overlapped meaning 
     *            that q_i+1 overlap a portion of the E(row_i, :).
     *            IN parallel E is splitten in vertical block over the threads  */
    /* SIDE RIGHT meaning apply E = E*Q = E * (q_1*q_2*.....*q_n) ==> so tarverse Vs in normal  order (forward) from q_1 to q_n 
     *            each q_i consist of applying V to a block of col E(:, col_i,:) and the applies are overlapped meaning 
     *            that q_i+1 overlap a portion of the E(:, col_i).
     *            IN parallel E is splitten in horizontal block over the threads  */
    
    //printf("  APPLY Q2   N %d  N_loc %d  nbchunk %d  NB %d  Vblksiz %d  SIDE %c \n", n, n_loc, nbchunk, nb, Vblksiz, side);
    for (magma_int_t i = 0; i<nbchunk; i++)
    {
        magma_int_t ib_loc = min(nb_loc, (n_loc - i*nb_loc));
        
        if(side=='L'){
            for (bg = nbGblk; bg>0; bg--)
            {
                firstcolj = (bg-1)*Vblksiz + 1;
                rownbm    = magma_ceildiv((n-(firstcolj+1)),nb);
                if(bg==nbGblk) rownbm    = magma_ceildiv((n-(firstcolj)),nb);  // last blk has size=1 used for complex to handle A(N,N-1)
                for (magma_int_t j = rownbm; j>0; j--)
                {
                    vlen = 0;
                    vnb  = 0;
                    colj      = (bg-1)*Vblksiz; // for k=0;I compute the fst and then can remove it from the loop
                    fst       = (rownbm -j)*nb+colj +1;
                    for (magma_int_t k=0; k<Vblksiz; k++)
                    {
                        colj     = (bg-1)*Vblksiz + k;
                        st       = (rownbm -j)*nb+colj +1;
                        ed       = min(st+nb-1,n-1);
                        if(st>ed)
                            break;
                        if((st==ed)&&(colj!=n-2))
                            break;
                        vlen=ed-fst+1;
                        vnb=k+1;
                    }        
                    colst     = (bg-1)*Vblksiz;
                    magma_bulge_findVTpos(n, nb, Vblksiz, colst, fst, ldv, ldt, &vpos, &tpos);
                    
                    if((vlen>0)&&(vnb>0)){
                        lapackf77_zlarfb( "L", "N", "F", "C", &vlen, &ib_loc, &vnb, V(vpos), &ldv, T(tpos), &ldt, E(fst,i*nb_loc), &lde, work, &ib_loc);
                    }
                    if(INFO!=0) 
                        printf("ERROR ZUNMQR INFO %d \n",INFO);
                }                
            }
        }else if (side=='R'){
            rownbm    = magma_ceildiv((n-1),nb);
            for (magma_int_t k = 1; k<=rownbm; k++)
            {
                ncolinvolvd = min(n-1, k*nb);
                avai_blksiz=min(Vblksiz,ncolinvolvd);
                nbgr = magma_ceildiv(ncolinvolvd,avai_blksiz);
                for (magma_int_t j = 1; j<=nbgr; j++)
                {
                    vlen = 0;
                    vnb  = 0;
                    cur_blksiz = min(ncolinvolvd-(j-1)*avai_blksiz, avai_blksiz);
                    colst = (j-1)*avai_blksiz;
                    coled = colst + cur_blksiz -1;
                    fst   = (rownbm -k)*nb+colst +1;
                    for (colj=colst; colj<=coled; colj++)
                    {
                        st       = (rownbm -k)*nb+colj +1;
                        ed       = min(st+nb-1,n-1);
                        if(st>ed)
                            break;
                        if((st==ed)&&(colj!=n-2))
                            break;
                        vlen=ed-fst+1;
                        vnb=vnb+1;
                    }    
                    magma_bulge_findVTpos(n, nb, Vblksiz, colst, fst, ldv, ldt, &vpos, &tpos);
                    if((vlen>0)&&(vnb>0)){
                        lapackf77_zlarfb( "R", "N", "F", "C", &ib_loc, &vlen, &vnb, V(vpos), &ldv, T(tpos), &ldt, E(i*nb_loc,fst), &lde, work, &ib_loc);       
                    }  
                }
            }
        }else{
            printf("ERROR SIDE %d \n",side);
        }
    } // END loop over the chunks
    
    magma_free_cpu(work);
    magma_free_cpu(work2);    
    
}
#undef E
#undef V
#undef TAU
#undef T
////////////////////////////////////////////////////////////////////////////////////////////////////

