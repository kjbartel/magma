#ifndef _PROFILING_H_
#define _PROFILING_H_

typedef struct measure_s {
    double sum;
    double sum2;
    long   n;
} measure_t;

void profiling_display_info(const char *kernel_name, measure_t perf[STARPU_NMAXWORKERS]);
void profiling_display_efficiency(void);

void profiling_zdisplay_all(void);
void profiling_cdisplay_all(void);
void profiling_ddisplay_all(void);
void profiling_sdisplay_all(void);

void MAGMA_zload_FakeModel();
void MAGMA_cload_FakeModel();
void MAGMA_dload_FakeModel();
void MAGMA_sload_FakeModel();

void MAGMA_zrestore_Model();
void MAGMA_crestore_Model();
void MAGMA_drestore_Model();
void MAGMA_srestore_Model();

#endif
