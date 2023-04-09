/******************************************************************************
 * * Copyright (C) 2023, Advanced Micro Devices, Inc. All rights reserved.
 * *******************************************************************************/

/*! @file fla_lapack_x86_common.c
 *  @brief Common front-end functions
 *         to choose optimized paths
 *  *  */

#include "FLAME.h"

#ifdef FLA_ENABLE_AMD_OPT

/* 3x3 Householder Rotation */
int fla_dhrot3(integer *n,
               doublereal *a, integer *lda,
               doublereal *v, doublereal *tau)
{
    fla_dhrot3_avx2(n, a, lda, v, tau);
    return 0;
}

/* 2x2 Plane Rotation */
int fla_drot(integer *n,
             doublereal *dx, integer *incx,
             doublereal *dy, integer *incy,
             doublereal *c__, doublereal *s)
{
    fla_drot_avx2(n, dx, incx, dy, incy, c__, s);
    return 0;
}

/* complex vector scaling when increment is 1 and specific threshold */
int fla_zscal(integer *n, doublecomplex *alpha,
              doublecomplex *x, integer *incx)
{
    /* Initialize global context data */
    aocl_fla_init();

    /* Take AVX path only for increment equal to 1 and particular threshold size*/
    if(global_context.is_avx2 && *incx == 1 && *n <= FLA_ZSCAL_INLINE_SMALL)
    {
        fla_zscal_ix1_avx2(n, alpha, x);
    }
    else
    {
        zscal_(n, (dcomplex *) alpha,(dcomplex *)  x, incx);
    }
    return 0;
}

/* scales a vector by a constant when threshold <= 128 */
int fla_dscal(integer *n, doublereal *da, doublereal *dx, integer *incx)
{
    /* Initialize global context data */
    aocl_fla_init();

    if(global_context.is_avx2 && *incx == 1 && *da != 0 && *n >= 1 && *n <= FLA_DSCAL_INLINE_SMALL)
    {
        fla_dscal_ix1_avx2(n, da, dx, incx);
    }
    else
    {
        dscal_(n, da, dx, incx);
    }    
    return 0;
}

/* Double QR (DGEQRF) for small sizes */
int fla_dgeqrf_small(integer *m, integer *n,
                     doublereal *a, integer *lda,
                     doublereal *tau, doublereal *work)
{
    fla_dgeqrf_small_avx2(m, n, a, lda, tau, work);
    return 0;
}
#endif
