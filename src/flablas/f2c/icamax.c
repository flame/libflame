/* icamax.f -- translated by f2c (version 19991025). You must link the resulting object file with the libraries: -lf2c -lm (in that order) */
#include "FLA_f2c.h"
integer icamax_(integer *n, complex *cx, integer *incx)
{
    /* System generated locals */
    integer ret_val, i__1, i__2;
    real r__1, r__2;
    /* Builtin functions */
    double r_imag(complex *);
    /* Local variables */
    real smax;
    integer i__, ix;
    /* finds the index of element having max. absolute value. */
    /* jack dongarra, linpack, 3/11/78. */
    /* modified 3/93 to return if incx .le. 0. */
    /* modified 12/3/93, array(1) declarations changed to array(*) */
    /* Parameter adjustments */
    --cx;
    /* Function Body */
    ret_val = 0;
    if (*n < 1 || *incx <= 0)
    {
        return ret_val;
    }
    ret_val = 1;
    if (*n == 1)
    {
        return ret_val;
    }
    if (*incx == 1)
    {
        goto L20;
    }
    /* code for increment not equal to 1 */
    ix = 1;
    smax = (r__1 = cx[1].r, f2c_abs(r__1)) + (r__2 = r_imag(&cx[1]), f2c_abs(r__2));
    ix += *incx;
    i__1 = *n;
    for (i__ = 2;
            i__ <= i__1;
            ++i__)
    {
        i__2 = ix;
        if ((r__1 = cx[i__2].r, f2c_abs(r__1)) + (r__2 = r_imag(&cx[ix]), f2c_abs( r__2)) <= smax)
        {
            goto L5;
        }
        ret_val = i__;
        i__2 = ix;
        smax = (r__1 = cx[i__2].r, f2c_abs(r__1)) + (r__2 = r_imag(&cx[ix]), f2c_abs( r__2));
L5:
        ix += *incx;
        /* L10: */
    }
    return ret_val;
    /* code for increment equal to 1 */
L20:
    smax = (r__1 = cx[1].r, f2c_abs(r__1)) + (r__2 = r_imag(&cx[1]), f2c_abs(r__2));
    i__1 = *n;
    for (i__ = 2;
            i__ <= i__1;
            ++i__)
    {
        i__2 = i__;
        if ((r__1 = cx[i__2].r, f2c_abs(r__1)) + (r__2 = r_imag(&cx[i__]), f2c_abs( r__2)) <= smax)
        {
            goto L30;
        }
        ret_val = i__;
        i__2 = i__;
        smax = (r__1 = cx[i__2].r, f2c_abs(r__1)) + (r__2 = r_imag(&cx[i__]), f2c_abs( r__2));
L30:
        ;
    }
    return ret_val;
}
/* icamax_ */

