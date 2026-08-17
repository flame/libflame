/*
 * Copyright (C) 2025-2026, Advanced Micro Devices, Inc. All rights reserved.
 */
/* cscal.f -- translated by f2c (version 19991025). You must link the resulting object file with the libraries: -lf2c -lm (in that order) */
#include "FLA_f2c.h"
/* Subroutine */
int cscal_(integer *n, scomplex *ca, scomplex *cx, integer * incx)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4;
    scomplex q__1;
    /* Local variables */
    integer i__, nincx;
    /* scales a vector by a constant. */
    /* jack dongarra, linpack, 3/11/78. */
    /* modified 3/93 to return if incx .le. 0. */
    /* modified 12/3/93, array(1) declarations changed to array(*) */
    /* Parameter adjustments */
    --cx;
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
    nincx = *n * *incx;
    i__1 = nincx;
    i__2 = *incx;
    for (i__ = 1;
            i__2 < 0 ? i__ >= i__1 : i__ <= i__1;
            i__ += i__2)
    {
        i__3 = i__;
        i__4 = i__;
        q__1.real = ca->real * cx[i__4].real - ca->imag * cx[i__4].imag, q__1.imag = ca->real * cx[ i__4].imag + ca->imag * cx[i__4].real;
        cx[i__3].real = q__1.real, cx[i__3].imag = q__1.imag;
        /* L10: */
    }
    return 0;
    /* code for increment equal to 1 */
L20:
    i__2 = *n;
    for (i__ = 1;
            i__ <= i__2;
            ++i__)
    {
        i__1 = i__;
        i__3 = i__;
        q__1.real = ca->real * cx[i__3].real - ca->imag * cx[i__3].imag, q__1.imag = ca->real * cx[ i__3].imag + ca->imag * cx[i__3].real;
        cx[i__1].real = q__1.real, cx[i__1].imag = q__1.imag;
        /* L30: */
    }
    return 0;
}
/* cscal_ */

