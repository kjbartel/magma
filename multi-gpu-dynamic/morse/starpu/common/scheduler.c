/**
 *
 * @file pzgeqrf.c
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.1.1
 * @author Jakub Kurzak
 * @author Hatem Ltaief
 * @date 2009-11-15
 *
 **/
#include <stdarg.h>
#include "morse_starpu.h"

#define MAXCALLBACKS        8

struct callback_wrappers {
    unsigned ncallbacks;
    void (*callback_funcs[MAXCALLBACKS])(void *arg);
    void *callback_args[MAXCALLBACKS];
};

static void callback_wrappers_func(void *arg)
{
    struct callback_wrappers *wrapper = arg;

    unsigned i;
    for (i = 0; i < wrapper->ncallbacks; i++)
    {
        if (wrapper->callback_funcs[i])
                wrapper->callback_funcs[i](wrapper->callback_args[i]);
    }

    free(wrapper);
}

void starpu_unpack_cl_args(void *_cl_arg, ...)
{
    unsigned char *cl_arg = _cl_arg;

    unsigned current_arg_offset = 0;
    va_list varg_list;

    va_start(varg_list, _cl_arg);

    /* We fill the different pointers with the appropriate arguments */
    unsigned char nargs = cl_arg[0];
    current_arg_offset += sizeof(char);

    unsigned arg;
    for (arg = 0; arg < nargs; arg++)
    {
        void *argptr = va_arg(varg_list, void *);
        size_t arg_size = *(size_t *)&cl_arg[current_arg_offset];
        current_arg_offset += sizeof(size_t);

        if (argptr)
                memcpy(argptr, &cl_arg[current_arg_offset], arg_size);

        current_arg_offset += arg_size;
    }

    va_end(varg_list);

    /* this should not really be done in StarPU but well .... */
    free(_cl_arg);
}

void starpu_Insert_Task(starpu_codelet *cl, ...)
{
    struct starpu_task *task = starpu_task_create();
    int arg_type;
    va_list varg_list;

    /* The buffer will contain : nargs, {size, content} (x nargs)*/

    /* Compute the size */
    size_t arg_buffer_size = 0;
    unsigned ncallbacks = 0;

    arg_buffer_size += sizeof(char);

    va_start(varg_list, cl);

    while( (arg_type = va_arg(varg_list, int)) != 0) {
       if (arg_type==INPUT || arg_type==INOUT || arg_type==OUTPUT || arg_type == SCRATCH || arg_type == REDUX)
       {
         va_arg(varg_list, starpu_data_handle);
       }
       else if (arg_type==VALUE) {
         va_arg(varg_list, void *);
         size_t cst_size = va_arg(varg_list, size_t);

         arg_buffer_size += sizeof(size_t);
         arg_buffer_size += cst_size;
       }
       else if (arg_type==CALLBACK) {
         va_arg(varg_list, void (*)(void *));
         va_arg(varg_list, void *);
         ncallbacks++;
       }
       else if (arg_type==PRIORITY) {
         va_arg(varg_list, int);
       }
    }

    va_end(varg_list);

    char *arg_buffer = malloc(arg_buffer_size);
    unsigned current_arg_offset = 0;

    /* We will begin the buffer with the number of args (which is stored as a char) */
    current_arg_offset += sizeof(char);
    int current_buffer = 0;
    unsigned char nargs = 0;

    struct callback_wrappers * cb_wrapper = NULL;

    if (ncallbacks > 1)
    {
        cb_wrapper = malloc(sizeof(struct callback_wrappers));
        cb_wrapper->ncallbacks = 0;

        task->callback_func = callback_wrappers_func;
        task->callback_arg = cb_wrapper;
    }

    va_start(varg_list, cl);

    while( (arg_type = va_arg(varg_list, int)) != 0) {
       
        if (arg_type==INPUT || arg_type==INOUT || arg_type==OUTPUT || arg_type == SCRATCH || arg_type == REDUX)
        {
           /* We have an access mode : we expect to find a handle */
           starpu_data_handle handle = va_arg(varg_list, starpu_data_handle);

           starpu_access_mode mode = STARPU_RW;
           switch (arg_type) {
             case INPUT:
               mode = STARPU_R;
               break;
             case INOUT:
               mode = STARPU_RW;
               break;
             case OUTPUT:
               mode = STARPU_W;
               break;
             case SCRATCH:
               mode = STARPU_SCRATCH;
               break;
             case REDUX:
               mode = STARPU_REDUX;
               break;
           }

           task->buffers[current_buffer].handle = handle;
           task->buffers[current_buffer].mode = mode;

           current_buffer++;
        }
        else if (arg_type==VALUE)
        {
          /* We have a constant value: this should be followed by a pointer to the cst value and the size of the constant */
          void *ptr = va_arg(varg_list, void *);
          size_t cst_size = va_arg(varg_list, size_t);

          *(size_t *)(&arg_buffer[current_arg_offset]) = cst_size;
          current_arg_offset += sizeof(size_t);

          memcpy(&arg_buffer[current_arg_offset], ptr, cst_size);
          current_arg_offset += cst_size;

          nargs++;
          STARPU_ASSERT(current_arg_offset <= arg_buffer_size);
        }
         else if (arg_type==CALLBACK) {
           void (*callback_func)(void *);
           callback_func = va_arg(varg_list, void (*)(void *));

           void *callback_arg;
           callback_arg = va_arg(varg_list, void *);
           if (ncallbacks == 1)
           {
             task->callback_func = callback_func;
             task->callback_arg = callback_arg;
           }
           else {
             unsigned callback_id = cb_wrapper->ncallbacks;
             cb_wrapper->callback_funcs[callback_id] = callback_func;
             cb_wrapper->callback_args[callback_id] = callback_arg;
             cb_wrapper->ncallbacks++;
           }
         }
        else if (arg_type==PRIORITY)
        {
          /* Followed by a priority level */
          int prio __attribute__((unused)) = va_arg(varg_list, int); 
#ifdef PRIO
          if (prio)
             task->priority = STARPU_MAX_PRIO;
#endif
        }
    }

    va_end(varg_list);

    arg_buffer[0] = nargs;

    STARPU_ASSERT(current_buffer == cl->nbuffers);

    task->cl = cl;
    task->cl_arg = arg_buffer;

    int ret = starpu_task_submit(task);

    if (STARPU_UNLIKELY(ret == -ENODEV))
        fprintf(stderr, "No one can execute task %p wih cl %p (symbol %s)\n", task, task->cl, (task->cl->model && task->cl->model->symbol)?task->cl->model->symbol:"none");

    STARPU_ASSERT(!ret);
}
