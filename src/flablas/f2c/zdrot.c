/*
 * Copyright (C) 2025-2026, Advanced Micro Devices, Inc. All rights reserved.
 */
/* zdrot.f -- translated by f2c (version 19991025). You must link the resulting object file with the libraries: -lf2c -lm (in that order) */
#include "FLA_f2c.h"
/* Subroutine */
int zdrot_(integer *n, dcomplex *zx, integer *incx, dcomplex *zy, integer *incy, doublereal *c__, doublereal *s)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4;
    dcomplex z__1, z__2, z__3;
    /* Local variables */
    integer i__;
    dcomplex ztemp;
    integer ix, iy;
    /* applies a plane rotation, where the cos and sin (c and s) are */
    /* double precision and the vectors zx and zy are double scomplex. */
    /* jack dongarra, linpack, 3/11/78. */
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
        z__2.real = *c__ * zx[i__2].real, z__2.imag = *c__ * zx[i__2].imag;
        i__3 = iy;
        z__3.real = *s * zy[i__3].real, z__3.imag = *s * zy[i__3].imag;
        z__1.real = z__2.real + z__3.real, z__1.imag = z__2.imag + z__3.imag;
        ztemp.real = z__1.real, ztemp.imag = z__1.imag;
        i__2 = iy;
        i__3 = iy;
        z__2.real = *c__ * zy[i__3].real, z__2.imag = *c__ * zy[i__3].imag;
        i__4 = ix;
        z__3.real = *s * zx[i__4].real, z__3.imag = *s * zx[i__4].imag;
        z__1.real = z__2.real - z__3.real, z__1.imag = z__2.imag - z__3.imag;
        zy[i__2].real = z__1.real, zy[i__2].imag = z__1.imag;
        i__2 = ix;
        zx[i__2].real = ztemp.real, zx[i__2].imag = ztemp.imag;
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
        z__2.real = *c__ * zx[i__2].real, z__2.imag = *c__ * zx[i__2].imag;
        i__3 = i__;
        z__3.real = *s * zy[i__3].real, z__3.imag = *s * zy[i__3].imag;
        z__1.real = z__2.real + z__3.real, z__1.imag = z__2.imag + z__3.imag;
        ztemp.real = z__1.real, ztemp.imag = z__1.imag;
        i__2 = i__;
        i__3 = i__;
        z__2.real = *c__ * zy[i__3].real, z__2.imag = *c__ * zy[i__3].imag;
        i__4 = i__;
        z__3.real = *s * zx[i__4].real, z__3.imag = *s * zx[i__4].imag;
        z__1.real = z__2.real - z__3.real, z__1.imag = z__2.imag - z__3.imag;
        zy[i__2].real = z__1.real, zy[i__2].imag = z__1.imag;
        i__2 = i__;
        zx[i__2].real = ztemp.real, zx[i__2].imag = ztemp.imag;
        /* L30: */
    }
    return 0;
}
/* zdrot_ */

