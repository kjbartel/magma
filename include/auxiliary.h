/*
    -- MAGMA (version 1.2.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       May 2012
*/

#ifndef _MAGMA_AUXILIARY_
#define _MAGMA_AUXILIARY_

typedef struct magma_timestr_s
{
  unsigned int sec;
  unsigned int usec;
} magma_timestr_t;

#ifdef __cplusplus
extern "C" {
#endif

int magma_get_spotrf_nb(int m);
int magma_get_sgeqrf_nb(int m);
int magma_get_sgeqlf_nb(int m);
int magma_get_sgetrf_nb(int m);
int magma_get_sgetri_nb(int m);
int magma_get_sgehrd_nb(int m);
int magma_get_ssytrd_nb(int m);
int magma_get_sgelqf_nb(int m);
int magma_get_sgebrd_nb(int m);
int magma_get_ssygst_nb(int m);

int magma_get_dpotrf_nb(int m);
int magma_get_dgeqrf_nb(int m);
int magma_get_dgeqlf_nb(int m);
int magma_get_dgetrf_nb(int m);
int magma_get_dgetri_nb(int m);
int magma_get_dgehrd_nb(int m);
int magma_get_dsytrd_nb(int m);
int magma_get_dgelqf_nb(int m);
int magma_get_dgebrd_nb(int m);
int magma_get_dsygst_nb(int m);

int magma_get_cpotrf_nb(int m);
int magma_get_cgetrf_nb(int m);
int magma_get_cgetri_nb(int m);
int magma_get_cgeqrf_nb(int m);
int magma_get_cgeqlf_nb(int m);
int magma_get_cgehrd_nb(int m);
int magma_get_chetrd_nb(int m);
int magma_get_cgelqf_nb(int m);
int magma_get_cgebrd_nb(int m);
int magma_get_chegst_nb(int m);

int magma_get_zpotrf_nb(int m);
int magma_get_zgetrf_nb(int m);
int magma_get_zgetri_nb(int m);
int magma_get_zgeqrf_nb(int m);
int magma_get_zgeqlf_nb(int m);
int magma_get_zgehrd_nb(int m);
int magma_get_zhetrd_nb(int m);
int magma_get_zgelqf_nb(int m);
int magma_get_zgebrd_nb(int m);
int magma_get_zhegst_nb(int m);

magma_timestr_t get_current_time(void);
double GetTimerValue(magma_timestr_t time_1, magma_timestr_t time_2);

double magma_wtime( void );

void printout_devices();

void swp2pswp(char trans, int n, int *ipiv, int *newipiv);

float getv(float *da);

size_t magma_strlcpy(char *dst, const char *src, size_t siz);

int magma_num_gpus( void );

#ifdef __cplusplus
}
#endif

#endif
