/*
    -- MAGMA (version 1.2.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2010

       @generated d Tue May 15 18:17:28 2012

*/

#include "common_magma.h"

/* === Define what BLAS to use ============================================ */
#define PRECISION_d
#if (defined(PRECISION_s) || defined(PRECISION_d))
  #define magma_dgemm magmablas_dgemm
  #define magma_dtrsm magmablas_dtrsm
#endif
#if (defined(PRECISION_z))
  #define magma_dgemm magmablas_dgemm
#endif
/* === End defining what BLAS to use ======================================= */

/* to appy pivoting from the previous big-panel: need some index-adjusting */
extern "C" void
magmablas_dpermute_long3( double *dAT, int lda, int *ipiv, int nb, int ind );

extern "C" magma_int_t
magma_dgetrf1_gpu(magma_int_t m, magma_int_t n,
                  double *dAT, magma_int_t ldda,
                  magma_int_t *ipiv, magma_int_t *info);

extern "C" magma_int_t
magma_dgetrf1_mgpu(magma_int_t num_gpus,
                   magma_int_t m, magma_int_t n, magma_int_t nb, magma_int_t offset,
                   double **d_lAT, magma_int_t lddat, magma_int_t *ipiv, 
                   double **d_lAP, double *work, magma_int_t lddwork, 
                   cudaStream_t **stream, magma_int_t *info);

extern "C" void
magmablas_dgetmatrix_transpose3(magma_int_t num_gpus, cudaStream_t **stream,
                  double **dat, int ldda,
                  double   *ha, int lda,
                  double  **dB, int lddb,
                  int m, int n , int nb);

extern "C" void
magmablas_dsetmatrix_transpose3(magma_int_t num_gpus, cudaStream_t** stream,
                  double  *ha,  int lda,
                  double **dat, int ldda, int starti,
                  double **dB,  int lddb,
                  int m, int n , int nb);

extern "C" magma_int_t
magma_dgetrf3_ooc(magma_int_t num_gpus0, magma_int_t m, magma_int_t n, double *a, magma_int_t lda, 
         magma_int_t *ipiv, magma_int_t *info)
{
/*  -- MAGMA (version 1.2.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2010

    Purpose
    =======
    DGETRF_OOC computes an LU factorization of a general M-by-N matrix A
    using partial pivoting with row interchanges.  This version does not
    require work space on the GPU passed as input. GPU memory is allocated
    in the routine. The matrix may not fit entirely in the GPU memory.

    The factorization has the form
       A = P * L * U
    where P is a permutation matrix, L is lower triangular with unit
    diagonal elements (lower trapezoidal if m > n), and U is upper
    triangular (upper trapezoidal if m < n).

    This is the right-looking Level 3 BLAS version of the algorithm.

    Note: The factorization of big panel is done calling multiple-gpu-interface.
    Pivots are applied on GPU within the big panel.

    Arguments
    =========
    M       (input) INTEGER
            The number of rows of the matrix A.  M >= 0.

    N       (input) INTEGER
            The number of columns of the matrix A.  N >= 0.

    A       (input/output) DOUBLE_PRECISION array, dimension (LDA,N)
            On entry, the M-by-N matrix to be factored.
            On exit, the factors L and U from the factorization
            A = P*L*U; the unit diagonal elements of L are not stored.

            Higher performance is achieved if A is in pinned memory, e.g.
            allocated using magma_malloc_host.

    LDA     (input) INTEGER
            The leading dimension of the array A.  LDA >= max(1,M).

    IPIV    (output) INTEGER array, dimension (min(M,N))
            The pivot indices; for 1 <= i <= min(M,N), row i of the
            matrix was interchanged with row IPIV(i).

    INFO    (output) INTEGER
            = 0:  successful exit
            < 0:  if INFO = -i, the i-th argument had an illegal value
                  or another error occured, such as memory allocation failed.
            > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
                  has been completed, but the factor U is exactly
                  singular, and division by zero will occur if it is used
                  to solve a system of equations.

    =====================================================================    */

#define    A(i,j) (a   + (j)*lda + (i))
#define inAT(d,i,j) (dAT[d] + (i)*nb*ldn_local + (j)*nb)
#define inPT(d,i,j) (dPT[d] + (i)*nb*nb + (j)*nb)
/* Flops formula */
//#define PROFILE
#ifdef  PROFILE
    double flops, time_rmajor = 0, time_rmajor2 = 0, time_rmajor3 = 0, time_mem = 0;
    magma_timestr_t start, start1, start2, end1, end, start0 = get_current_time();
#define FMULS_GETRF(__m, __n) ( ((__m) < (__n)) ? (0.5 * (__m) * ((__m) * ((__n) - (1./3.) * (__m) - 1. ) + (__n)) + (2. / 3.) * (__m)) \
                                :                 (0.5 * (__n) * ((__n) * ((__m) - (1./3.) * (__n) - 1. ) + (__m)) + (2. / 3.) * (__n)) )
#define FADDS_GETRF(__m, __n) ( ((__m) < (__n)) ? (0.5 * (__m) * ((__m) * ((__n) - (1./3.) * (__m)      ) - (__n)) + (1. / 6.) * (__m)) \
                                :                 (0.5 * (__n) * ((__n) * ((__m) - (1./3.) * (__n)      ) - (__m)) + (1. / 6.) * (__n)) )
#define PRECISION_d
#if defined(PRECISION_z) || defined(PRECISION_c)
#define FLOPS(m, n) ( 6. * FMULS_GETRF(m, n) + 2. * FADDS_GETRF(m, n) )
#else
#define FLOPS(m, n) (      FMULS_GETRF(m, n) +      FADDS_GETRF(m, n) )
#endif
#endif
    double    *dAT[4], *dA[4], *dPT[4];
    double    c_one     = MAGMA_D_ONE;
    double    c_neg_one = MAGMA_D_NEG_ONE;
    magma_int_t        iinfo = 0, nb, maxm, n_local[4], ldn_local;
    magma_int_t        N, M, NB, NBk, I, d, num_gpus;
    magma_int_t        i, ii, jj, h = 2, offset, ib, rows, s;
#if CUDA_VERSION > 3010
    size_t totalMem;
#else
    unsigned int totalMem;
#endif
    CUdevice dev;
    static cudaStream_t stream[4][2];

    *info = 0;

    if (m < 0)
        *info = -1;
    else if (n < 0)
        *info = -2;
    else if (lda < max(1,m))
        *info = -4;

    if (*info != 0) {
        magma_xerbla( __func__, -(*info) );
        return *info;
    }

    /* Quick return if possible */
    if (m == 0 || n == 0)
        return *info;

    /* initialize nb */
    nb = magma_get_dgetrf_nb(m);
    maxm = ((m  + 31)/32)*32;

    /* figure out NB */
    cuDeviceGet( &dev, 0);
    cuDeviceTotalMem( &totalMem, dev );
    totalMem /= sizeof(double);
    
    /* number of columns in the big panel */
    NB = (magma_int_t)(0.8*totalMem/maxm-h*nb); 
    char * ngr_nb_char = getenv("MAGMA_NGR_NB");
    if( ngr_nb_char != NULL ) NB = max( nb, min( NB, atoi(ngr_nb_char) ) );

    if( num_gpus0 > ceil((double)NB/nb) ) {
      num_gpus = (int)ceil((double)NB/nb);
    } else {
      num_gpus = num_gpus0;
    }
    if( num_gpus*NB >= n ) {
#ifdef CHECK_DGETRF_OOC
      printf( "      * still fit in GPU memory.\n" );
#endif
      NB = n;
    } else {
#ifdef CHECK_DGETRF_OOC
      printf( "      * don't fit in GPU memory.\n" );
#endif
      NB = num_gpus*NB;
      NB = max(nb,(NB / nb) * nb); /* making sure it's devisable by nb (x64) */
    }

#ifdef CHECK_DGETRF_OOC
    if( NB != n ) printf( "      * running in out-core mode (n=%d, NB=%d, nb=%d).\n",n,NB,nb );
    else          printf( "      * running in in-core mode  (n=%d, NB=%d, nb=%d).\n",n,NB,nb );
    fflush(stdout);
#endif 

    if ( (nb <= 1) || (nb >= min(m,n)) ) {
        /* Use CPU code for scalar of one tile. */
        lapackf77_dgetrf(&m, &n, a, &lda, ipiv, info);
    } else {
        /* Use hybrid blocked code. */

    /* allocate memory on GPU to store the big panel */
#ifdef  PROFILE
    start = get_current_time();
#endif
    n_local[0] = (NB/nb)/num_gpus;
    if( NB%(nb*num_gpus) != 0 ) n_local[0] ++;
    n_local[0] *= nb;
    ldn_local = ((n_local[0]+31)/32)*32;

    for( d=0; d<num_gpus; d++ ) {
      magma_setdevice(d);
      if (MAGMA_SUCCESS != magma_dmalloc( &dA[d], (h*nb + ldn_local)*maxm )) {
        *info = MAGMA_ERR_DEVICE_ALLOC;
        return *info;
      }
      dPT[d] = dA[d] + nb*maxm; /* for storing the previous panel from CPU          */
      dAT[d] = dA[d] + h*nb*maxm;
      for( ii=0; ii<h; ii++ ) {
        magma_queue_create( &stream[d][ii] );
      }
    }
#ifdef PROFILE
    end = get_current_time();
    printf( " memory-allocation time: %e\n",GetTimerValue(start, end)/1000.0 );
    start = get_current_time();
#endif
    for( I=0; I<n; I+=NB ) {
          M = m;
          N = min( NB, n-I );       /* number of columns in this big panel             */
          s = min(max(m-I,0),N)/nb; /* number of small block-columns in this big panel */

          maxm = ((M + 31)/32)*32;
          if( num_gpus0 > ceil((double)N/nb) ) {
            num_gpus = (int)ceil((double)N/nb);
          } else {
            num_gpus = num_gpus0;
          }

          for( d=0; d<num_gpus; d++ ) {
            n_local[d] = ((N/nb)/num_gpus)*nb;
            if (d < (N/nb)%num_gpus)
              n_local[d] += nb;
            else if (d == (N/nb)%num_gpus)
              n_local[d] += N%nb;
          }
          ldn_local = ((n_local[0]+31)/32)*32;

#ifdef PROFILE
          start2 = get_current_time();
#endif
          /* upload the next big panel into GPU, transpose (A->A'), and pivot it */
          magmablas_dsetmatrix_transpose3(num_gpus, (cudaStream_t **)stream, A(0,I), lda, 
                              dAT, ldn_local, 0, dA, maxm, M, N, nb);
          //magmablas_dhtodt3(num_gpus, (cudaStream_t **)stream, A(0,I), lda, 
          //                    dAT, ldn_local, dA, maxm, M, N, nb, h);
          //magmablas_dhtodt4(num_gpus, (cudaStream_t **)stream, A(0,I), lda, 
          //                    dAT, ldn_local, dA, maxm, M, N, nb, h/2);

#ifdef PROFILE
          start1 = get_current_time();
#endif
          /* == --------------------------------------------------------------- == */
          /* == loop around the previous big-panels to update the new big-panel == */
          for( offset = 0; offset<min(m,I); offset+=NB ) 
          {
            NBk = min( m-offset, NB );
            /* start sending the first tile from the previous big-panels to gpus */
            for( d=0; d<num_gpus; d++ ) {
              magma_setdevice(d);
              magma_dsetmatrix_async( (M-offset), nb,
                                      A(offset,offset), lda,
                                      dA[d],            (maxm-offset), stream[d][0] );
            }

            /* applying the pivot from the previous big-panel */
            for( d=0; d<num_gpus; d++ ) {
              magma_setdevice(d);
              magmablas_dpermute_long3( inAT(d,0,0), ldn_local, ipiv, NBk, offset );
            }

            /* == going through each block-column of previous big-panels == */
            for( jj=0, ib=offset/nb; jj<NBk; jj+=nb, ib++ ) 
            {
              ii   = offset+jj;
              rows = maxm - ii;
              for( d=0; d<num_gpus; d++ ) {
                magma_setdevice(d);

                /* upload the previous block-column to GPU */
                magma_queue_sync( stream[d][0] );

                /* transpose the previous column */
                magmablas_dtranspose2( inPT(d,0,0), nb, dA[d], rows, M-ii, nb);

                /* start sending the next column */
                if( jj+nb < NBk )
                magma_dsetmatrix_async( (M-ii-nb), min(nb,NBk-jj-nb),
                                        A(ii+nb,ii+nb), lda,
                                        dA[d],          (rows-nb), stream[d][0] );

                /* update with the block column */
                magma_dtrsm( MagmaRight, MagmaUpper, MagmaNoTrans, MagmaUnit, 
                         n_local[d], nb, c_one, inPT(d,0,0), nb, inAT(d,ib,0), ldn_local );
                if( M > ii+nb ) {
                  magma_dgemm( MagmaNoTrans, MagmaNoTrans, 
                      n_local[d], M-(ii+nb), nb, c_neg_one, inAT(d,ib,0), ldn_local, 
                      inPT(d,1,0), nb, c_one, inAT(d,ib+1,0), ldn_local );
                }  
              } /* end of for each block-columns in a big-panel */
            }
          } /* end of for each previous big-panels */

          /* calling magma-gpu interface to panel-factorize the big panel */
          if( M > I ) {
            magma_dgetrf1_mgpu(num_gpus, M-I, N, nb, I, dAT, ldn_local, ipiv+I, dA, &a[I*lda], lda,
                               (cudaStream_t **)stream, &iinfo);
            if( iinfo < 0 ) {
              *info = iinfo;
              break;
            } else if( iinfo != 0 ) {
              *info = iinfo + I * NB;		
              //break;
            }
            /* adjust pivots */
            for( ii=I; ii<min(I+N,m); ii++ ) ipiv[ii] += I;
          }
#ifdef PROFILE
          end1 = get_current_time();
          time_rmajor  += GetTimerValue(start1, end1);
          time_rmajor3 += GetTimerValue(start2, end1);
          time_mem += (GetTimerValue(start2, end1)-GetTimerValue(start1, end1))/1000.0;
#endif      
          /* download the current big panel to CPU */
          magmablas_dgetmatrix_transpose3(num_gpus, (cudaStream_t **)stream, dAT, ldn_local, A(0,I), lda,  dA, maxm, M, N, nb);
#ifdef PROFILE
          end1 = get_current_time();
          time_rmajor2 += GetTimerValue(start1, end1);
#endif

    } /* end of for */

#ifdef PROFILE
    end = get_current_time();
    flops = FLOPS( (double)m, (double)n ) / 1000000;
    printf(" NB=%d nb=%d\n",NB,nb); 
    printf(" memcopy and transpose %e seconds\n",time_mem );
    printf(" total time %e seconds\n",GetTimerValue(start0,end)/1000.0);
    printf(" Performance %f GFlop/s, %f seconds without htod and dtoh\n",     flops / time_rmajor,  time_rmajor /1000.0);
    printf(" Performance %f GFlop/s, %f seconds with    htod\n",              flops / time_rmajor3, time_rmajor3/1000.0);
    printf(" Performance %f GFlop/s, %f seconds with    dtoh\n",              flops / time_rmajor2, time_rmajor2/1000.0);
    printf(" Performance %f GFlop/s, %f seconds without memory-allocation\n", flops / GetTimerValue(start, end), GetTimerValue(start,end)/1000.0);
#endif

    for( d=0; d<num_gpus0; d++ ) {
      magma_setdevice(d);
      magma_free( dA[d] ); 
      for( ii=0; ii<h; ii++ ) {
        magma_queue_destroy( stream[d][ii] );
      }
    }
    }
    
    return *info;
} /* magma_dgetrf_ooc */




extern "C" magma_int_t
magma_dgetrf2_piv(magma_int_t num_gpus0, magma_int_t m, magma_int_t n, double *a, magma_int_t lda, 
                  magma_int_t *ipiv, magma_int_t *info)
{
    magma_int_t nb, h = 2, num_gpus;
    magma_int_t NB, I, k1, k2, incx, minmn, maxm;

    *info = 0;

    if (m < 0)
        *info = -1;
    else if (n < 0)
        *info = -2;
    else if (lda < max(1,m))
        *info = -4;

    if (*info != 0)
        return *info;

    /* Quick return if possible */
    if (m == 0 || n == 0)
        return *info;

    /* initialize nb */
    nb = magma_get_dgetrf_nb(m);
    maxm = ((m  + 31)/32)*32;

    /* figure out NB */
#if CUDA_VERSION > 3010
    size_t totalMem;
#else
    unsigned int totalMem;
#endif
    CUdevice dev;
    cuDeviceGet( &dev, 0);
    cuDeviceTotalMem( &totalMem, dev );
    totalMem /= sizeof(double);

    /* number of columns in the big panel */
    NB = (magma_int_t)(0.8*totalMem/maxm-h*nb); 
    //NB = (magma_int_t)min(n,num_gpus*(0.8*totalMem/maxm-h*nb)); 
    //NB = (magma_int_t)min(n,(num_gpus*0.8*totalMem/(maxm))-2*nb); 
    char * ngr_nb_char = getenv("MAGMA_NGR_NB");
    if( ngr_nb_char != NULL ) NB = max( nb, min( NB, atoi(ngr_nb_char) ) );

    if( num_gpus0 > ceil((double)NB/nb) ) {
      num_gpus = (int)ceil((double)NB/nb);
    } else {
      num_gpus = num_gpus0;
    }
    if( num_gpus*NB >= n ) {
#ifdef CHECK_DGETRF_OOC
      printf( "      * still fit in GPU memory.\n" );
#endif
      NB = n;
    } else {
#ifdef CHECK_DGETRF_OOC
      printf( "      * don't fit in GPU memory.\n" );
#endif
      NB = num_gpus*NB;
      NB = max(nb,(NB / nb) * nb); /* making sure it's devisable by nb (x64) */
    }
    minmn = min(m,n);

    for( I=0; I<minmn-NB; I+=NB ) {
      k1 = 1+I+NB;
      k2 = minmn;
      incx = 1;
      lapackf77_dlaswp(&NB, &a[I*lda], &lda, &k1, &k2, ipiv, &incx);
    }

    return *info;
} /* magma_dgetrf_piv */

#undef inAT
#undef inPT
#undef A

