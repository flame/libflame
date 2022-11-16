/******************************************************************************
 * * Copyright (C) 2022, Advanced Micro Devices, Inc. All rights reserved.
 * *******************************************************************************/

/*! @file fla_lapack_x86_common.c
 *  @brief Common front-end functions
 *         to choose optimized paths
 *  *  */

#include "FLAME.h"

#ifdef FLA_ENABLE_AMD_OPT

/* 3x3 Householder Rotation */
int fla_dhrot3(integer *n, doublereal *a, integer *lda, doublereal *v, doublereal *tau)
{
    fla_dhrot3_avx2(n, a, lda, v, tau);
    return 0;
}

/* 2x2 Plane Rotation */
int fla_drot(integer *n, doublereal *dx, integer *incx, doublereal *dy, integer *incy, doublereal *c__, doublereal *s)
{
    fla_drot_avx2(n, dx, incx, dy, incy, c__, s);
    return 0;
}
#endif
