/*
    -- MAGMA (version 1.4.0-beta2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       June 2013

       @author Raffaele Solca

*/
#ifdef SETAFFINITY

#include "affinity.h"

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include <sched.h>

affinity_set::affinity_set()
{
    CPU_ZERO(&set);
}

affinity_set::affinity_set(int cpu_nr)
{
    CPU_ZERO(&set);
    CPU_SET(cpu_nr, &set);
}

void affinity_set::add(int cpu_nr)
{
    CPU_SET(cpu_nr, &set);
}

int affinity_set::get_affinity()
{
    return sched_getaffinity( 0, sizeof(set), &set);
}

int affinity_set::set_affinity()
{
    return sched_setaffinity( 0, sizeof(set), &set);
}

void affinity_set::print_affinity(int id, const char* s)
{
    if (get_affinity() == 0)
        print_set(id, s);
    else
        printf("Error in sched_getaffinity\n");

}

void affinity_set::print_set(int id, const char* s)
{
    char cpustring[1024];

    int cpu_count=CPU_COUNT(&set);
    int charcnt = 0;

    charcnt = sprintf(cpustring,"thread %d has affinity with %d CPUS: ", id, cpu_count);

    int nrcpu=0;

    for(int icpu=0; nrcpu<cpu_count && icpu<CPU_SETSIZE; ++icpu){
        if(CPU_ISSET(icpu,&set)){
            charcnt += sprintf(&(cpustring[charcnt]),"%d,",icpu);
            ++nrcpu;
        }
    }
    charcnt += sprintf(&(cpustring[charcnt-1]),"\n") - 1; // -1 is used to remove "," after last cpu.
    printf("%s: %s", s, cpustring);
    fflush(stdout);
}

#endif

