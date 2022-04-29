/*
    Copyright (c) 2022 Advanced Micro Devices, Inc. All rights reserved.
*/

#include "FLA_f2c.h"

aocl_fla_progress_callback aocl_fla_progress_ptr = NULL;

void aocl_fla_set_progress(aocl_fla_progress_callback func)
{
    aocl_fla_progress_ptr = func;
}