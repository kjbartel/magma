#include <math.h>
#include "morse_starpu.h"

void profiling_display_info(const char *kernel_name, measure_t perf[STARPU_NMAXWORKERS])
{
    int header = 1;
    unsigned worker;
    for (worker = 0; worker < starpu_worker_get_count(); worker++)
    {
        if (perf[worker].n > 0)
        {
            if ( header ) {
                fprintf(stderr, "Performance for kernel %s\n", kernel_name);
                header = 0;
            }
            char workername[128];
            starpu_worker_get_name(worker, workername, 128);
            
            long   n    = perf[worker].n;
            double sum  = perf[worker].sum;
            double sum2 = perf[worker].sum2;
            
            double avg = sum / n;
            double sd  = sqrt((sum2 - (sum*sum)/n)/n);

            fprintf(stderr, "\t%s\t%.2lf\t%.2lf\t%ld\n", workername, avg, sd, n);
        }
    }
}

void profiling_display_efficiency(void)
{
    fprintf(stderr, "Efficiency\n");

    double max_total_time = 0.0;
    unsigned worker;

    for (worker = 0; worker < starpu_worker_get_count(); worker++)
    {
        char workername[128];
        starpu_worker_get_name(worker, workername, 128);
        
        struct starpu_worker_profiling_info info;
        starpu_worker_get_profiling_info(worker, &info);
        
        double executing_time = starpu_timing_timespec_to_us(&info.executing_time);
        double total_time = starpu_timing_timespec_to_us(&info.total_time);
        
        max_total_time = (total_time > max_total_time)?total_time:max_total_time;
        
        float overhead = 100.0 - (100.0*executing_time/total_time);
        fprintf(stderr, "\t%s\ttotal %.2lf s\texec %.2lf s\toverhead %.2lf%%\n", 
                workername, total_time*1e-6, executing_time*1e-6, overhead);
    }

    fprintf(stderr, "Total execution time: %.2lf us\n", max_total_time);
}

void morse_schedprofile_display(void)
{
    fprintf(stderr, "\n");
    profiling_display_efficiency();
    
    /* Display bus consumption */
    starpu_bus_profiling_helper_display_summary();
}

