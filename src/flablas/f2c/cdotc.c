/*
 * Copyright (C) 2025-2026, Advanced Micro Devices, Inc. All rights reserved.
 */
/* cdotc.f -- translated by f2c (version 19991025). You must link the resulting object file with the libraries: -lf2c -lm (in that order) */
#include "FLA_f2c.h"
scomplex cdotc_(/*scomplex * ret_val,*/
    integer *n, scomplex *cx, integer *incx, scomplex *cy, integer *incy)
{
    scomplex ret_val;
    /* System generated locals */
    integer i__1, i__2;
    scomplex q__1, q__2, q__3;
    /* Builtin functions */
    void r_cnjg(scomplex *, scomplex *);
    /* Local variables */
    integer i__;
    scomplex ctemp;
    integer ix, iy;
    /* forms the dot product of two vectors, conjugating the first */
    /* vector. */
    /* jack dongarra, linpack, 3/11/78. */
    /* modified 12/3/93, array(1) declarations changed to array(*) */
    /* Parameter adjustments */
    --cy;
    --cx;
    /* Function Body */
    ctemp.real = 0.f, ctemp.imag = 0.f;
    ret_val.real = 0.f, ret_val.imag = 0.f;
    if (*n <= 0)
    {
        return ret_val ;
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
        r_cnjg(&q__3, &cx[ix]);
        i__2 = iy;
        q__2.real = q__3.real * cy[i__2].real - q__3.imag * cy[i__2].imag, q__2.imag = q__3.real * cy[i__2].imag + q__3.imag * cy[i__2].real;
        q__1.real = ctemp.real + q__2.real, q__1.imag = ctemp.imag + q__2.imag;
        ctemp.real = q__1.real, ctemp.imag = q__1.imag;
        ix += *incx;
        iy += *incy;
        /* L10: */
    }
    ret_val.real = ctemp.real, ret_val.imag = ctemp.imag;
    return ret_val ;
    /* code for both increments equal to 1 */
L20:
    i__1 = *n;
    for (i__ = 1;
            i__ <= i__1;
            ++i__)
    {
        r_cnjg(&q__3, &cx[i__]);
        i__2 = i__;
        q__2.real = q__3.real * cy[i__2].real - q__3.imag * cy[i__2].imag, q__2.imag = q__3.real * cy[i__2].imag + q__3.imag * cy[i__2].real;
        q__1.real = ctemp.real + q__2.real, q__1.imag = ctemp.imag + q__2.imag;
        ctemp.real = q__1.real, ctemp.imag = q__1.imag;
        /* L30: */
    }
    ret_val.real = ctemp.real, ret_val.imag = ctemp.imag;
    return ret_val ;
}
/* cdotc_ */

