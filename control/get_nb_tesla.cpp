/*
    -- MAGMA (version 1.3.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2012
*/

// ==== Definition of blocking sizes for Nvidia Tesla (pre-Fermi) cards
#if (GPUSHMEM < 200)

/* ////////////////////////////////////////////////////////////////////////////
   -- Return nb for potrf based on n
*/ 
extern "C"
int magma_get_spotrf_nb(int n) {
    if (n <= 3328)
        return 128;
    else if (n <= 4256)
        return 224;
    else 
        return 288;
}

extern "C"
int magma_get_dpotrf_nb(int n) {
    if (n <= 3328)
        return 128;
    else if (n <= 4256)
        return 128;
    else 
        return 256;
}

extern "C"
int magma_get_cpotrf_nb(int n) {
    return 64;
}

extern "C"
int magma_get_zpotrf_nb(int n) {
    return 64;
}

/* ////////////////////////////////////////////////////////////////////////////
   -- Return nb for geqp3 based on m
*/
extern "C"
int magma_get_sgeqp3_nb(int m) {
    return 32;
}

extern "C"
int magma_get_dgeqp3_nb(int m) {
    return 32;
}

extern "C"
int magma_get_cgeqp3_nb(int m) {
    return 32;
}

extern "C"
int magma_get_zgeqp3_nb(int m) {
    return 32;
}

/* ////////////////////////////////////////////////////////////////////////////
   -- Return nb for geqrf based on m
*/
extern "C"
int magma_get_sgeqrf_nb(int m) {
    if (m <= 2048)
        return 32;
    else if (m <= 4032)
        return 64;
    else
        return 128;
}

extern "C"
int magma_get_dgeqrf_nb(int m) {
    if (m <= 2048)
        return 64;
    else 
        return 128;
}

extern "C"
int magma_get_cgeqrf_nb(int m) {
    if (m <= 2048)
        return 32;
    else if (m <= 4032)
        return 64;
    else
        return 128;
}

extern "C"
int magma_get_zgeqrf_nb(int m) {
    if (m <= 1024)
        return 64;
    else
        return 128;
}

/* ////////////////////////////////////////////////////////////////////////////
   -- Return nb for geqlf based on m
*/
extern "C"
int magma_get_sgeqlf_nb(int m) {
    if (m <= 1024)
        return 32;
    else if (m <= 4032)
        return 64;
    else
        return 128;
}

extern "C"
int magma_get_dgeqlf_nb(int m) {
    if (m <= 1024)
        return 32;
    else if (m <= 4032)
        return 64;
    else
        return 128;
}

extern "C"
int magma_get_cgeqlf_nb(int m) {
    if (m <= 2048)
        return 32;
    else if (m <= 4032)
        return 64;
    else
        return 128;
}

extern "C"
int magma_get_zgeqlf_nb(int m) {
    if (m <= 1024)
        return 64;
    else
        return 128;
}

/* ////////////////////////////////////////////////////////////////////////////
   -- Return nb for gelqf based on m
*/
extern "C"
int magma_get_sgelqf_nb(int m) {
    return magma_get_sgeqrf_nb(m);
}

extern "C"
int magma_get_dgelqf_nb(int m) {
    if (m <= 2048)
        return 32;
    else if (m <= 4032)
        return 64;
    else
        return 128;
}

extern "C"
int magma_get_cgelqf_nb(int m) {
    if (m <= 2048)
        return 32;
    else if (m <= 4032)
        return 64;
    else
        return 128;
}

extern "C"
int magma_get_zgelqf_nb(int m) {
    if (m <= 1024)
        return 64;
    else
        return 128;
}

/* ////////////////////////////////////////////////////////////////////////////
   -- Return nb for getrf based on m;
      the return value should be multiple of 64
*/
extern "C"
int magma_get_sgetrf_nb(int m) {
    if (m <= 2048)
        return 64;
    else
        return 128;
}

extern "C"
int magma_get_dgetrf_nb(int m) {
    if (m <= 2048)
        return 64;
    else
        return 128;
}

extern "C"
int magma_get_cgetrf_nb(int m) {
    if (m <= 2048)
        return 64;
    else
        return 128;
}

extern "C"
int magma_get_zgetrf_nb(int n) {
    return 128;
}

/* ////////////////////////////////////////////////////////////////////////////
   -- Return nb for gehrd based on m;
      the return value should be a multiple of 32
*/
extern "C"
int magma_get_sgehrd_nb(int m) {
    if (m <= 1024)
        return 32;
    else
        return 64;
}

extern "C"
int magma_get_dgehrd_nb(int m) {
    if (m <= 2048)
        return 32;
    else
        return 64;
}

extern "C"
int magma_get_cgehrd_nb(int m) {
    if (m <= 1024)
        return 32;
    else
        return 64;
}

extern "C"
int magma_get_zgehrd_nb(int m) {
    if (m <= 2048)
        return 32;
    else
        return 64;
}

/* ////////////////////////////////////////////////////////////////////////////
   -- Return nb for sytrd based on m;
      the return value should be a multiple of 32
*/
extern "C"
int magma_get_ssytrd_nb(int m) {
    return 32;
}

extern "C"
int magma_get_dsytrd_nb(int m) {
    return 32;
    //return 64;
}

extern "C"
int magma_get_chetrd_nb(int m) {
    return 32;
}

extern "C"
int magma_get_zhetrd_nb(int m) {
    return 32;
    //return 64;
}

/* ////////////////////////////////////////////////////////////////////////////
   -- Return nb for gebrd based on m;
      the return value should be a multiple of 32
*/
extern "C"
int magma_get_sgebrd_nb(int m) {
    return 32;
}

extern "C"
int magma_get_dgebrd_nb(int m) {
    return 32;
}

extern "C"
int magma_get_cgebrd_nb(int m) {
    return 32;
}

extern "C"
int magma_get_zgebrd_nb(int m) {
    return 32;
}

/* ////////////////////////////////////////////////////////////////////////////
   -- Return nb for zhegst based on m;
      the return value should be a multiple of 32
*/
extern "C"
int magma_get_ssygst_nb(int m) {
    return 64;
}

extern "C"
int magma_get_dsygst_nb(int m) {
    return 64;
}

extern "C"
int magma_get_chegst_nb(int m) {
    return 64;
}

extern "C"
int magma_get_zhegst_nb(int m) {
    return 64;
}

/* ////////////////////////////////////////////////////////////////////////////
   -- Return nb for getri based on m
*/
extern "C"
int magma_get_sgetri_nb(int m) {
    return 64;
}

extern "C"
int magma_get_dgetri_nb(int m) {
    return 64;
}

extern "C"
int magma_get_cgetri_nb(int m) {
    return 64;
}

extern "C"
int magma_get_zgetri_nb(int m) {
    return 64;
}

/* ////////////////////////////////////////////////////////////////////////////
   -- Return nb for gesvd based on m
*/

extern "C"
int magma_get_sgesvd_nb(int m){
    return magma_get_sgebrd_nb(m);
}

extern "C"
int magma_get_dgesvd_nb(int m){
    return magma_get_dgebrd_nb(m);
}

extern "C"
int magma_get_cgesvd_nb(int m){
    return magma_get_cgebrd_nb(m);
}

extern "C"
int magma_get_zgesvd_nb(int m){
    return magma_get_zgebrd_nb(m);
}

#endif  // HAVE_CUBLAS
