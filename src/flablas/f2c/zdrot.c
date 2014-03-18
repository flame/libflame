/* zdrot.f -- translated by f2c (version 19991025). You must link the resulting object file with the libraries: -lf2c -lm (in that order) */
#include "FLA_f2c.h"
/* Subroutine */
int zdrot_(integer *n, doublecomplex *zx, integer *incx, doublecomplex *zy, integer *incy, doublereal *c__, doublereal *s)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4;
    doublecomplex z__1, z__2, z__3;
    /* Local variables */
    integer i__;
    doublecomplex ztemp;
    integer ix, iy;
    /* applies a plane rotation, where the cos and sin (c and s) are */
    /* double precision and the vectors zx and zy are double complex. */
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
        z__2.r = *c__ * zx[i__2].r, z__2.i = *c__ * zx[i__2].i;
        i__3 = iy;
        z__3.r = *s * zy[i__3].r, z__3.i = *s * zy[i__3].i;
        z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
        ztemp.r = z__1.r, ztemp.i = z__1.i;
        i__2 = iy;
        i__3 = iy;
        z__2.r = *c__ * zy[i__3].r, z__2.i = *c__ * zy[i__3].i;
        i__4 = ix;
        z__3.r = *s * zx[i__4].r, z__3.i = *s * zx[i__4].i;
        z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - z__3.i;
        zy[i__2].r = z__1.r, zy[i__2].i = z__1.i;
        i__2 = ix;
        zx[i__2].r = ztemp.r, zx[i__2].i = ztemp.i;
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
        z__2.r = *c__ * zx[i__2].r, z__2.i = *c__ * zx[i__2].i;
        i__3 = i__;
        z__3.r = *s * zy[i__3].r, z__3.i = *s * zy[i__3].i;
        z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
        ztemp.r = z__1.r, ztemp.i = z__1.i;
        i__2 = i__;
        i__3 = i__;
        z__2.r = *c__ * zy[i__3].r, z__2.i = *c__ * zy[i__3].i;
        i__4 = i__;
        z__3.r = *s * zx[i__4].r, z__3.i = *s * zx[i__4].i;
        z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - z__3.i;
        zy[i__2].r = z__1.r, zy[i__2].i = z__1.i;
        i__2 = i__;
        zx[i__2].r = ztemp.r, zx[i__2].i = ztemp.i;
        /* L30: */
    }
    return 0;
}
/* zdrot_ */

