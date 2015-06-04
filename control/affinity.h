/*
    -- MAGMA (version 1.4.0-beta2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       June 2013

       @author Raffaele Solca

*/
#include "common_magma.h"

#ifndef _MAGMA_AFFINITY
#define _MAGMA_AFFINITY

#ifdef SETAFFINITY

class affinity_set
{

public:

    affinity_set();

    affinity_set(int cpu_nr);

    void add(int cpu_nr);

    int get_affinity();

    int set_affinity();

    void print_affinity(int id, const char* s);

    void print_set(int id, const char* s);

private:

    cpu_set_t set;
};

#endif

#endif

