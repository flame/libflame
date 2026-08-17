/*
 * Copyright (C) 2025-2026, Advanced Micro Devices, Inc. All rights reserved.
 */
/* zswap.f -- translated by f2c (version 19991025). You must link the resulting object file with the libraries: -lf2c -lm (in that order) */
#include "FLA_f2c.h"
/* Subroutine */
int zswap_(integer *n, dcomplex *zx, integer *incx, dcomplex *zy, integer *incy)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    /* Local variables */
    integer i__;
    dcomplex ztemp;
    integer ix, iy;
    /* interchanges two vectors. */
    /* jack dongarra, 3/11/78. */
    /* modified 12/3/93, array(1) declarations changed to array(*) */
    /* Parameter adjustments */
    --zy;
    --zx;
    /* Function Body */
    if (*n <= 0)
    {
        return 0;
    }
    if (*incx == 1 && *incy == 1)
    {
        goto L20;
    }
    /* code for unequal increments or equal increments not equal */
    /* to 1 */
    ix = 1;
    iy = 1;
    if (*incx < 0)
    {
        ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0)
    {
        iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i__ = 1;
            i__ <= i__1;
            ++i__)
    {
        i__2 = ix;
        ztemp.real = zx[i__2].real, ztemp.imag = zx[i__2].imag;
        i__2 = ix;
        i__3 = iy;
        zx[i__2].real = zy[i__3].real, zx[i__2].imag = zy[i__3].imag;
        i__2 = iy;
        zy[i__2].real = ztemp.real, zy[i__2].imag = ztemp.imag;
        ix += *incx;
        iy += *incy;
        /* L10: */
    }
    return 0;
    /* code for both increments equal to 1 */
L20:
    i__1 = *n;
    for (i__ = 1;
            i__ <= i__1;
            ++i__)
    {
        i__2 = i__;
        ztemp.real = zx[i__2].real, ztemp.imag = zx[i__2].imag;
        i__2 = i__;
        i__3 = i__;
        zx[i__2].real = zy[i__3].real, zx[i__2].imag = zy[i__3].imag;
        i__2 = i__;
        zy[i__2].real = ztemp.real, zy[i__2].imag = ztemp.imag;
        /* L30: */
    }
    return 0;
}
/* zswap_ */

