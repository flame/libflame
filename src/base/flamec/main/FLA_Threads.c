/*
    Copyright (C) 2022-2026, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "FLAME.h"

#ifdef FLA_OPENMP_MULTITHREADING

#include <omp.h>

/* To determine the sub partition range of current thread */
void FLA_Thread_get_subrange(int thread_ID, int num_threads, fla_dim_t range, fla_dim_t *sub_range,
                             fla_dim_t *index)
{
    fla_dim_t sub_region, remainder;

    if(range <= 0)
    {
        *sub_range = 0;
        *index = 0;
    }
    else
    {
        sub_region = range / num_threads;
        remainder = range % num_threads;

        /* divide row/column region equally among each thread*/
        if(thread_ID < remainder)
        {
            *sub_range = sub_region + 1;
            *index = thread_ID * (*sub_range);
        }
        else
        {
            *sub_range = sub_region;
            *index = (*sub_range == 0) ? 0 : (remainder + thread_ID * sub_region);
        }
    }
}

/* To determine the sub partition range of current thread */
void FLA_Thread_get_subrange_chunks(int thread_ID, int num_threads, size_t elem_size_in_bytes,
                                    fla_dim_t range, fla_dim_t *sub_range, fla_dim_t *index,
                                    fla_dim_t *thread_threshold)
{
    fla_dim_t thread_required, x;
    fla_dim_t granule;
    fla_dim_t granule_times_nt;

    /* Guard against invalid caller input (e.g. truncated thread count == 0). */
    if(num_threads <= 0)
    {
        *sub_range = 0;
        *index = 0;
        *thread_threshold = 0;
        return;
    }

    if(elem_size_in_bytes == 0)
    {
        granule = 1;
    }
    else
    {
        granule = (fla_dim_t)(FLA_CACHE_LINE_SIZE_BYTES / elem_size_in_bytes);
        if(granule < 1)
            granule = 1;
    }

    if(range <= 0)
    {
        *sub_range = 0;
        *index = 0;
        *thread_threshold = 0;
    }
    else
    {
        granule_times_nt = granule * (fla_dim_t)num_threads;
        x = range / granule_times_nt;
        *sub_range = (x + ((range % granule_times_nt) > 0 ? 1 : 0)) * granule;
        thread_required = range / *sub_range;
        if((range % *sub_range) > 0)
        {
            thread_required++;
        }
        *thread_threshold = thread_required;
        if(thread_ID < (thread_required - 1))
        {
            *index = thread_ID * (*sub_range);
        }
        else if(thread_ID == (thread_required - 1))
        {
            *index = (thread_ID * (*sub_range));
            *sub_range = (range - (*index));
        }
        else
        {
            *index = 0;
            *sub_range = 0;
        }
    }
}

/* To determine optimum thread number for a give API */
void FLA_Thread_optimum(API_ID family, int *actual_num_threads)
{
    int optimal_num_threads = 0;
    extern int fla_thread_get_num_threads();

    switch(family)
    {
        case FLA_LABRD:
            optimal_num_threads = 8;
            break;
        case FLA_ORMQR:
            optimal_num_threads = 16;
            break;
        case FLA_ORMLQ:
            optimal_num_threads = 16;
            break;
        default:
            optimal_num_threads = 0;
            break;
    }

    *actual_num_threads = fla_thread_get_num_threads();

    if(optimal_num_threads && (*actual_num_threads > optimal_num_threads))
        *actual_num_threads = optimal_num_threads;

    return;
}

#endif
