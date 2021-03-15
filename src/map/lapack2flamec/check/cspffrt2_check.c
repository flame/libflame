/*
    Copyright (c) 2020 Advanced Micro Devices, Inc.  All rights reserved.
    Nov 06, 2020
*/

#include "FLA_lapack2flame_return_defs.h"
#include "FLA_f2c.h"

int cspffrt2_check(scomplex *ap, int *n, int * ncolm, scomplex *work, scomplex *work2)
{
    int ret_val = LAPACK_SUCCESS;

    if (*n < 0)
    {
        ret_val = LAPACK_FAILURE;
    }
    else if (*ncolm < 0 || *ncolm > *n)
    {
        ret_val = LAPACK_FAILURE;
    }
    /* Quick return if possible */
    if (*n == 0 || *ncolm == 0)
    {
        ret_val =  LAPACK_QUICK_RETURN;
    }
    return ret_val;
}

