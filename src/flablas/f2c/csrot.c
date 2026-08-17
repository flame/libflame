/*
 * Copyright (C) 2025-2026, Advanced Micro Devices, Inc. All rights reserved.
 */
/* csrot.f -- translated by f2c (version 19991025). You must link the resulting object file with the libraries: -lf2c -lm (in that order) */
#include "FLA_f2c.h"
/* Subroutine */
int csrot_(integer *n, scomplex *cx, integer *incx, scomplex * cy, integer *incy, real *c__, real *s)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4;
    scomplex q__1, q__2, q__3;
    /* Local variables */
    integer i__;
    scomplex ctemp;
    integer ix, iy;
    /* applies a plane rotation, where the cos and sin (c and s) are real */
    /* and the vectors cx and cy are scomplex. */
    /* jack dongarra, linpack, 3/11/78. */
    /* Parameter adjustments */
    --cy;
    --cx;
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
        q__2.real = *c__ * cx[i__2].real, q__2.imag = *c__ * cx[i__2].imag;
        i__3 = iy;
        q__3.real = *s * cy[i__3].real, q__3.imag = *s * cy[i__3].imag;
        q__1.real = q__2.real + q__3.real, q__1.imag = q__2.imag + q__3.imag;
        ctemp.real = q__1.real, ctemp.imag = q__1.imag;
        i__2 = iy;
        i__3 = iy;
        q__2.real = *c__ * cy[i__3].real, q__2.imag = *c__ * cy[i__3].imag;
        i__4 = ix;
        q__3.real = *s * cx[i__4].real, q__3.imag = *s * cx[i__4].imag;
        q__1.real = q__2.real - q__3.real, q__1.imag = q__2.imag - q__3.imag;
        cy[i__2].real = q__1.real, cy[i__2].imag = q__1.imag;
        i__2 = ix;
        cx[i__2].real = ctemp.real, cx[i__2].imag = ctemp.imag;
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
        q__2.real = *c__ * cx[i__2].real, q__2.imag = *c__ * cx[i__2].imag;
        i__3 = i__;
        q__3.real = *s * cy[i__3].real, q__3.imag = *s * cy[i__3].imag;
        q__1.real = q__2.real + q__3.real, q__1.imag = q__2.imag + q__3.imag;
        ctemp.real = q__1.real, ctemp.imag = q__1.imag;
        i__2 = i__;
        i__3 = i__;
        q__2.real = *c__ * cy[i__3].real, q__2.imag = *c__ * cy[i__3].imag;
        i__4 = i__;
        q__3.real = *s * cx[i__4].real, q__3.imag = *s * cx[i__4].imag;
        q__1.real = q__2.real - q__3.real, q__1.imag = q__2.imag - q__3.imag;
        cy[i__2].real = q__1.real, cy[i__2].imag = q__1.imag;
        i__2 = i__;
        cx[i__2].real = ctemp.real, cx[i__2].imag = ctemp.imag;
        /* L30: */
    }
    return 0;
}
/* csrot_ */

