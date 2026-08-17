/*
 * Copyright (C) 2025-2026, Advanced Micro Devices, Inc. All rights reserved.
 */
/* zdscal.f -- translated by f2c (version 19991025). You must link the resulting object file with the libraries: -lf2c -lm (in that order) */
#include "FLA_f2c.h"
/* Subroutine */
int zdscal_(integer *n, doublereal *da, dcomplex *zx, integer *incx)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    dcomplex z__1, z__2;
    /* Local variables */
    integer i__, ix;
    /* scales a vector by a constant. */
    /* jack dongarra, 3/11/78. */
    /* modified 3/93 to return if incx .le. 0. */
    /* modified 12/3/93, array(1) declarations changed to array(*) */
    /* Parameter adjustments */
    --zx;
    /* Function Body */
    if (*n <= 0 || *incx <= 0)
    {
        return 0;
    }
    if (*incx == 1)
    {
        goto L20;
    }
    /* code for increment not equal to 1 */
    ix = 1;
    i__1 = *n;
    for (i__ = 1;
            i__ <= i__1;
            ++i__)
    {
        i__2 = ix;
        z__2.real = *da, z__2.imag = 0.;
        i__3 = ix;
        z__1.real = z__2.real * zx[i__3].real - z__2.imag * zx[i__3].imag, z__1.imag = z__2.real * zx[i__3].imag + z__2.imag * zx[i__3].real;
        zx[i__2].real = z__1.real, zx[i__2].imag = z__1.imag;
        ix += *incx;
        /* L10: */
    }
    return 0;
    /* code for increment equal to 1 */
L20:
    i__1 = *n;
    for (i__ = 1;
            i__ <= i__1;
            ++i__)
    {
        i__2 = i__;
        z__2.real = *da, z__2.imag = 0.;
        i__3 = i__;
        z__1.real = z__2.real * zx[i__3].real - z__2.imag * zx[i__3].imag, z__1.imag = z__2.real * zx[i__3].imag + z__2.imag * zx[i__3].real;
        zx[i__2].real = z__1.real, zx[i__2].imag = z__1.imag;
        /* L30: */
    }
    return 0;
}
/* zdscal_ */

