/*
 * Copyright (C) 2025-2026, Advanced Micro Devices, Inc. All rights reserved.
 */
/* zdotc.f -- translated by f2c (version 19991025). You must link the resulting object file with the libraries: -lf2c -lm (in that order) */
#include "FLA_f2c.h"
dcomplex zdotc_(/*dcomplex * ret_val,*/
    integer *n, dcomplex *zx, integer *incx, dcomplex *zy, integer *incy)
{
    dcomplex ret_val;
    /* System generated locals */
    integer i__1, i__2;
    dcomplex z__1, z__2, z__3;
    /* Builtin functions */
    void d_cnjg(dcomplex *, dcomplex *);
    /* Local variables */
    integer i__;
    dcomplex ztemp;
    integer ix, iy;
    /* forms the dot product of a vector. */
    /* jack dongarra, 3/11/78. */
    /* modified 12/3/93, array(1) declarations changed to array(*) */
    /* Parameter adjustments */
    --zy;
    --zx;
    /* Function Body */
    ztemp.real = 0., ztemp.imag = 0.;
    ret_val.real = 0., ret_val.imag = 0.;
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
        d_cnjg(&z__3, &zx[ix]);
        i__2 = iy;
        z__2.real = z__3.real * zy[i__2].real - z__3.imag * zy[i__2].imag, z__2.imag = z__3.real * zy[i__2].imag + z__3.imag * zy[i__2].real;
        z__1.real = ztemp.real + z__2.real, z__1.imag = ztemp.imag + z__2.imag;
        ztemp.real = z__1.real, ztemp.imag = z__1.imag;
        ix += *incx;
        iy += *incy;
        /* L10: */
    }
    ret_val.real = ztemp.real, ret_val.imag = ztemp.imag;
    return ret_val ;
    /* code for both increments equal to 1 */
L20:
    i__1 = *n;
    for (i__ = 1;
            i__ <= i__1;
            ++i__)
    {
        d_cnjg(&z__3, &zx[i__]);
        i__2 = i__;
        z__2.real = z__3.real * zy[i__2].real - z__3.imag * zy[i__2].imag, z__2.imag = z__3.real * zy[i__2].imag + z__3.imag * zy[i__2].real;
        z__1.real = ztemp.real + z__2.real, z__1.imag = ztemp.imag + z__2.imag;
        ztemp.real = z__1.real, ztemp.imag = z__1.imag;
        /* L30: */
    }
    ret_val.real = ztemp.real, ret_val.imag = ztemp.imag;
    return ret_val ;
}
/* zdotc_ */

