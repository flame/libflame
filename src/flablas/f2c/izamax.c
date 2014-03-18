/* izamax.f -- translated by f2c (version 19991025). You must link the resulting object file with the libraries: -lf2c -lm (in that order) */
#include "FLA_f2c.h"
integer izamax_(integer *n, doublecomplex *zx, integer *incx)
{
    /* System generated locals */
    integer ret_val, i__1;
    /* Local variables */
    doublereal smax;
    integer i__;
    extern doublereal dcabs1_(doublecomplex *);
    integer ix;
    /* finds the index of element having max. absolute value. */
    /* jack dongarra, 1/15/85. */
    /* modified 3/93 to return if incx .le. 0. */
    /* modified 12/3/93, array(1) declarations changed to array(*) */
    /* Parameter adjustments */
    --zx;
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
    smax = dcabs1_(&zx[1]);
    ix += *incx;
    i__1 = *n;
    for (i__ = 2;
            i__ <= i__1;
            ++i__)
    {
        if (dcabs1_(&zx[ix]) <= smax)
        {
            goto L5;
        }
        ret_val = i__;
        smax = dcabs1_(&zx[ix]);
L5:
        ix += *incx;
        /* L10: */
    }
    return ret_val;
    /* code for increment equal to 1 */
L20:
    smax = dcabs1_(&zx[1]);
    i__1 = *n;
    for (i__ = 2;
            i__ <= i__1;
            ++i__)
    {
        if (dcabs1_(&zx[i__]) <= smax)
        {
            goto L30;
        }
        ret_val = i__;
        smax = dcabs1_(&zx[i__]);
L30:
        ;
    }
    return ret_val;
}
/* izamax_ */

