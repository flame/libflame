/*
    Copyright (C) 2022, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "FLAME.h"

#ifdef FLA_ENABLE_MULTITHREADING

#include <omp.h>


/* To determine the sub partition range of current thread */
void FLA_Thread_get_subrange
     (
       int thread_ID,
       int num_threads,
       integer range,
       integer *sub_range,
       integer *index
     )
{
    integer sub_region, remainder;

    sub_region = range/num_threads;
    remainder = range%num_threads;

    /* divide row/column region equally among each thread*/
    if(thread_ID < remainder)
    {
        *sub_range = sub_region + 1;
        *index = thread_ID * (*sub_range);
    }
    else
    {
        *sub_range = sub_region;
        *index = remainder + thread_ID * sub_region;
    }
}

/* To determine optimum thread number for a give API */
void FLA_Thread_optimum( API_ID  family, int *actual_num_threads)
{
    int optimal_num_threads = 0;

    switch(family)
    {
        case FLA_LABRD : 
            optimal_num_threads = 8;
            break;
        case FLA_ORMQR :
            optimal_num_threads = 16;
            break;
        case FLA_ORMLQ:
            optimal_num_threads = 16;
            break;
        default :
            optimal_num_threads = 0;
            break;
    }

    *actual_num_threads = omp_get_max_threads();

    if(optimal_num_threads && *actual_num_threads > optimal_num_threads)
        *actual_num_threads = optimal_num_threads;

    return;
}

#endif