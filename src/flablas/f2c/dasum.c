/* dasum.f -- translated by f2c (version 19991025). You must link the resulting object file with the libraries: -lf2c -lm (in that order) */
#include "FLA_f2c.h"
doublereal dasum_(integer *n, doublereal *dx, integer *incx)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal ret_val, d__1, d__2, d__3, d__4, d__5, d__6;
    /* Local variables */
    integer i__, m;
    doublereal dtemp;
    integer nincx, mp1;
    /* takes the sum of the absolute values. */
    /* jack dongarra, linpack, 3/11/78. */
    /* modified 3/93 to return if incx .le. 0. */
    /* modified 12/3/93, array(1) declarations changed to array(*) */
    /* Parameter adjustments */
    --dx;
    /* Function Body */
    ret_val = 0.;
    dtemp = 0.;
    if (*n <= 0 || *incx <= 0)
    {
        return ret_val;
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
        dtemp += (d__1 = dx[i__], f2c_abs(d__1));
        /* L10: */
    }
    ret_val = dtemp;
    return ret_val;
    /* code for increment equal to 1 */
    /* clean-up loop */
L20:
    m = *n % 6;
    if (m == 0)
    {
        goto L40;
    }
    i__2 = m;
    for (i__ = 1;
            i__ <= i__2;
            ++i__)
    {
        dtemp += (d__1 = dx[i__], f2c_abs(d__1));
        /* L30: */
    }
    if (*n < 6)
    {
        goto L60;
    }
L40:
    mp1 = m + 1;
    i__2 = *n;
    for (i__ = mp1;
            i__ <= i__2;
            i__ += 6)
    {
        dtemp = dtemp + (d__1 = dx[i__], f2c_abs(d__1)) + (d__2 = dx[i__ + 1], f2c_abs(d__2)) + (d__3 = dx[i__ + 2], f2c_abs(d__3)) + (d__4 = dx[i__ + 3], f2c_abs(d__4)) + (d__5 = dx[i__ + 4], f2c_abs(d__5)) + (d__6 = dx[i__ + 5], f2c_abs(d__6));
        /* L50: */
    }
L60:
    ret_val = dtemp;
    return ret_val;
}
/* dasum_ */

