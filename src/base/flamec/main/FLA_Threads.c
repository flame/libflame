/*
    Copyright (C) 2022, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "FLAME.h"
#include<omp.h>


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