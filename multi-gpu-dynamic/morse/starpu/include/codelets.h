/**
 *
 * @file codelets.h
 *
 *  MAGMA codelets kernel
 *  MAGMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver,
 *  and INRIA Bordeaux Sud-Ouest
 *
 * @version 2.3.1
 * @author Cedric Augonnet
 * @author Mathieu Faverge
 * @date 2011-06-01
 *
 **/

#ifndef _CODELETS_H_
#define _CODELETS_H_

#include "codelet_profile.h"

#define CODELETS_ALL(name, _nbuffers, cpu_func_name, cuda_func_name, mc_func_name)        \
struct starpu_perfmodel_t cl_##name##_fake = {                                        \
    .type   = STARPU_HISTORY_BASED,                                                \
    .symbol = "fake_"#name                                                        \
};                                                                                \
                                                                                \
static struct starpu_perfmodel_t cl_##name##_model = {                                \
    .type   = STARPU_HISTORY_BASED,                                                \
    .symbol = ""#name                                                                \
};                                                                                \
                                                                                \
static struct starpu_perfmodel_t cl_##name##_mc_model = {                        \
    .type   = STARPU_HISTORY_BASED,                                                \
    .symbol = "mc_"#name                                                        \
};                                                                                \
                                                                                \
starpu_codelet cl_##name = {                                                        \
    .where     = STARPU_CPU|STARPU_CUDA,                                        \
    .cpu_func  = ((cpu_func_name)),                                                \
    .cuda_func = ((cuda_func_name)),                                                \
    .nbuffers  = ((_nbuffers)),                                                        \
    .model     = &cl_##name##_model                                                \
};                                                                                \
                                                                                \
starpu_codelet cl_##name##_mc = {                                                \
    .where     = STARPU_CPU|STARPU_CUDA,                                        \
    .cpu_func  = ((mc_func_name)),                                                \
    .cuda_func = ((cuda_func_name)),                                                \
    .nbuffers  = ((_nbuffers)),                                                        \
    .model     = &cl_##name##_mc_model                                                \
};                                                                                \
                                                                                \
starpu_codelet cl_##name##_cpu = {                                                \
    .where     = STARPU_CPU,                                                        \
    .cpu_func  = ((cpu_func_name)),                                                \
    .nbuffers  = ((_nbuffers)),                                                        \
    .model     = &cl_##name##_model                                                \
};                                                                                \
                                                                                \
void cl_##name##_restrict_where(uint32_t where)                                        \
{                                                                                \
    cl_##name.where = where;                                                        \
}                                                                                \
                                                                                \
void cl_##name##_restore_where(void)                                                \
{                                                                                \
    cl_##name.where = STARPU_CPU|STARPU_CUDA;                                        \
}                                                                                \
                                                                                \
void cl_##name##_load_fake_model(void)                                                \
{                                                                                \
    cl_##name.model     = &cl_##name##_fake;                                        \
    cl_##name##_mc.model  = &cl_##name##_fake;                                        \
    cl_##name##_cpu.model = &cl_##name##_fake;                                        \
}                                                                                \
                                                                                \
void cl_##name##_restore_model(void)                                                \
{                                                                                \
    cl_##name.model     = &cl_##name##_model;                                        \
    cl_##name##_mc.model  = &cl_##name##_mc_model;                                \
    cl_##name##_cpu.model = &cl_##name##_model;                                        \
}                                                                                \


#define CODELETS_CPU(name, _nbuffers, cpu_func_name)                                \
struct starpu_perfmodel_t cl_##name##_fake = {                                        \
    .type   = STARPU_HISTORY_BASED,                                                \
    .symbol = "fake_"#name                                                        \
};                                                                                \
                                                                                \
static struct starpu_perfmodel_t cl_##name##_model = {                                \
    .type   = STARPU_HISTORY_BASED,                                                \
    .symbol = ""#name                                                                \
};                                                                                \
                                                                                \
starpu_codelet cl_##name = {                                                        \
    .where     = STARPU_CPU,                                                        \
    .cpu_func  = ((cpu_func_name)),                                                \
    .nbuffers  = ((_nbuffers)),                                                        \
    .model     = &cl_##name##_model                                                \
};                                                                                \
                                                                                \
void cl_##name##_load_fake_model(void)                                                \
{                                                                                \
    cl_##name.model     = &cl_##name##_fake;                                        \
}                                                                                \
                                                                                \
void cl_##name##_restore_model(void)                                                \
{                                                                                \
    cl_##name.model     = &cl_##name##_model;                                        \
}                                                                               \


#define CODELETS_ALL_HEADER(name)                                               \
CL_CALLBACK_HEADER(name)                                                        \
void cl_##name##_load_fake_model(void);                                                \
void cl_##name##_restore_model(void);                                                \
extern starpu_codelet cl_##name;                                                \
extern starpu_codelet cl_##name##_mc;                                                \
extern starpu_codelet cl_##name##_cpu;                                                \
void cl_##name##_restrict_where(uint32_t where);                                \
void cl_##name##_restore_where(void);

#define CODELETS_CPU_HEADER(name)                                                \
CL_CALLBACK_HEADER(name)                                                        \
void cl_##name##_load_fake_model(void);                                                \
void cl_##name##_restore_model(void);                                                \
extern starpu_codelet cl_##name;

#if defined(MORSE_USE_CUDA)
#define CODELETS(name, _nbuffers, cpu_func_name, cuda_func_name, mc_func_name)  \
  CODELETS_ALL(name, _nbuffers, cpu_func_name, cuda_func_name, mc_func_name)

#define CODELETS_HEADER(name)  CODELETS_ALL_HEADER(name)
#else
#define CODELETS(name, _nbuffers, cpu_func_name, cuda_func_name, mc_func_name)  \
  CODELETS_CPU(name, _nbuffers, cpu_func_name)

#define CODELETS_HEADER(name)  CODELETS_CPU_HEADER(name)
#endif

#define SCODELETS_HEADER(name)                CODELETS_HEADER(s##name)
#define DCODELETS_HEADER(name)                CODELETS_HEADER(d##name)
#define CCODELETS_HEADER(name)                CODELETS_HEADER(c##name)
#define ZCODELETS_HEADER(name)                CODELETS_HEADER(z##name)

#define SCODELETS_CPU_HEADER(name)        CODELETS_CPU_HEADER(s##name)
#define DCODELETS_CPU_HEADER(name)        CODELETS_CPU_HEADER(d##name)
#define CCODELETS_CPU_HEADER(name)        CODELETS_CPU_HEADER(c##name)
#define ZCODELETS_CPU_HEADER(name)        CODELETS_CPU_HEADER(z##name)

#include "codelet_z.h"
#include "codelet_c.h"
#include "codelet_d.h"
#include "codelet_s.h"

#endif /* _CODELETS_H_ */
