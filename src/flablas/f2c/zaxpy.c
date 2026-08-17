/*
 * Copyright (C) 2025-2026, Advanced Micro Devices, Inc. All rights reserved.
 */
/* zaxpy.f -- translated by f2c (version 19991025). You must link the resulting object file with the libraries: -lf2c -lm (in that order) */
#include "FLA_f2c.h"
/* Subroutine */
int zaxpy_(integer *n, dcomplex *za, dcomplex *zx, integer *incx, dcomplex *zy, integer *incy)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4;
    dcomplex z__1, z__2;
    /* Local variables */
    integer i__;
    extern doublereal dcabs1_(dcomplex *);
    integer ix, iy;
    /* constant times a vector plus a vector. */
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
    if (dcabs1_(za) == 0.)
    {
        return 0;
    }
    if (*incx == 1 && *incy == 1)
    {
        goto L20;
    }
    /* code for unequal increments or equal increments */
    /* not equal to 1 */
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
        i__2 = iy;
        i__3 = iy;
        i__4 = ix;
        z__2.real = za->real * zx[i__4].real - za->imag * zx[i__4].imag, z__2.imag = za->real * zx[ i__4].imag + za->imag * zx[i__4].real;
        z__1.real = zy[i__3].real + z__2.real, z__1.imag = zy[i__3].imag + z__2.imag;
        zy[i__2].real = z__1.real, zy[i__2].imag = z__1.imag;
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
        i__3 = i__;
        i__4 = i__;
        z__2.real = za->real * zx[i__4].real - za->imag * zx[i__4].imag, z__2.imag = za->real * zx[ i__4].imag + za->imag * zx[i__4].real;
        z__1.real = zy[i__3].real + z__2.real, z__1.imag = zy[i__3].imag + z__2.imag;
        zy[i__2].real = z__1.real, zy[i__2].imag = z__1.imag;
        /* L30: */
    }
    return 0;
}
/* zaxpy_ */

