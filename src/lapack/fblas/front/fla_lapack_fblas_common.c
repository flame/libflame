/******************************************************************************
 * * Copyright (C) 2023, Advanced Micro Devices, Inc. All rights reserved.
 * *******************************************************************************/

/*! @file fla_lapack_fblas_common.c
 *  @brief Common front-end functions
 *         to choose optimized paths
 *  *  */

#include "FLAME.h"

#ifdef FLA_ENABLE_AMD_OPT

/* IDAMAX for small sizes (<= 128)*/
integer fla_idamax(integer *n, doublereal *dx, integer *incx)
{
    return fla_idamax_small(n, dx, incx);
}
#endif
