/* ../netlib/clacrt.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief \b CLACRT performs a linear transformation of a pair of complex vectors. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CLACRT + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clacrt. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clacrt. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clacrt. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CLACRT( N, CX, INCX, CY, INCY, C, S ) */
/* .. Scalar Arguments .. */
/* INTEGER INCX, INCY, N */
/* COMPLEX C, S */
/* .. */
/* .. Array Arguments .. */
/* COMPLEX CX( * ), CY( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CLACRT performs the operation */
/* > */
/* > ( c s )( x ) ==> ( x ) */
/* > ( -s c )( y ) ( y ) */
/* > */
/* > where c and s are complex and the vectors x and y are complex. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The number of elements in the vectors CX and CY. */
/* > \endverbatim */
/* > */
/* > \param[in,out] CX */
/* > \verbatim */
/* > CX is COMPLEX array, dimension (N) */
/* > On input, the vector x. */
/* > On output, CX is overwritten with c*x + s*y. */
/* > \endverbatim */
/* > */
/* > \param[in] INCX */
/* > \verbatim */
/* > INCX is INTEGER */
/* > The increment between successive values of CX. INCX <> 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] CY */
/* > \verbatim */
/* > CY is COMPLEX array, dimension (N) */
/* > On input, the vector y. */
/* > On output, CY is overwritten with -s*x + c*y. */
/* > \endverbatim */
/* > */
/* > \param[in] INCY */
/* > \verbatim */
/* > INCY is INTEGER */
/* > The increment between successive values of CY. INCY <> 0. */
/* > \endverbatim */
/* > */
/* > \param[in] C */
/* > \verbatim */
/* > C is COMPLEX */
/* > \endverbatim */
/* > */
/* > \param[in] S */
/* > \verbatim */
/* > S is COMPLEX */
/* > C and S define the matrix */
/* > [ C S ]. */
/* > [ -S C ] */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date September 2012 */
/* > \ingroup complexOTHERauxiliary */
/* ===================================================================== */
/* Subroutine */
int clacrt_(integer *n, complex *cx, integer *incx, complex * cy, integer *incy, complex *c__, complex *s)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4;
    complex q__1, q__2, q__3;
    /* Local variables */
    integer i__, ix, iy;
    complex ctemp;
    /* -- LAPACK auxiliary routine (version 3.4.2) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* September 2012 */
    /* .. Scalar Arguments .. */
    /* .. */
    /* .. Array Arguments .. */
    /* .. */
    /* ===================================================================== */
    /* .. Local Scalars .. */
    /* .. */
    /* .. Executable Statements .. */
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
    /* Code for unequal increments or equal increments not equal to 1 */
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
        q__2.r = c__->r * cx[i__2].r - c__->i * cx[i__2].i;
        q__2.i = c__->r * cx[i__2].i + c__->i * cx[i__2].r; // , expr subst
        i__3 = iy;
        q__3.r = s->r * cy[i__3].r - s->i * cy[i__3].i;
        q__3.i = s->r * cy[ i__3].i + s->i * cy[i__3].r; // , expr subst
        q__1.r = q__2.r + q__3.r;
        q__1.i = q__2.i + q__3.i; // , expr subst
        ctemp.r = q__1.r;
        ctemp.i = q__1.i; // , expr subst
        i__2 = iy;
        i__3 = iy;
        q__2.r = c__->r * cy[i__3].r - c__->i * cy[i__3].i;
        q__2.i = c__->r * cy[i__3].i + c__->i * cy[i__3].r; // , expr subst
        i__4 = ix;
        q__3.r = s->r * cx[i__4].r - s->i * cx[i__4].i;
        q__3.i = s->r * cx[ i__4].i + s->i * cx[i__4].r; // , expr subst
        q__1.r = q__2.r - q__3.r;
        q__1.i = q__2.i - q__3.i; // , expr subst
        cy[i__2].r = q__1.r;
        cy[i__2].i = q__1.i; // , expr subst
        i__2 = ix;
        cx[i__2].r = ctemp.r;
        cx[i__2].i = ctemp.i; // , expr subst
        ix += *incx;
        iy += *incy;
        /* L10: */
    }
    return 0;
    /* Code for both increments equal to 1 */
L20:
    i__1 = *n;
    for (i__ = 1;
            i__ <= i__1;
            ++i__)
    {
        i__2 = i__;
        q__2.r = c__->r * cx[i__2].r - c__->i * cx[i__2].i;
        q__2.i = c__->r * cx[i__2].i + c__->i * cx[i__2].r; // , expr subst
        i__3 = i__;
        q__3.r = s->r * cy[i__3].r - s->i * cy[i__3].i;
        q__3.i = s->r * cy[ i__3].i + s->i * cy[i__3].r; // , expr subst
        q__1.r = q__2.r + q__3.r;
        q__1.i = q__2.i + q__3.i; // , expr subst
        ctemp.r = q__1.r;
        ctemp.i = q__1.i; // , expr subst
        i__2 = i__;
        i__3 = i__;
        q__2.r = c__->r * cy[i__3].r - c__->i * cy[i__3].i;
        q__2.i = c__->r * cy[i__3].i + c__->i * cy[i__3].r; // , expr subst
        i__4 = i__;
        q__3.r = s->r * cx[i__4].r - s->i * cx[i__4].i;
        q__3.i = s->r * cx[ i__4].i + s->i * cx[i__4].r; // , expr subst
        q__1.r = q__2.r - q__3.r;
        q__1.i = q__2.i - q__3.i; // , expr subst
        cy[i__2].r = q__1.r;
        cy[i__2].i = q__1.i; // , expr subst
        i__2 = i__;
        cx[i__2].r = ctemp.r;
        cx[i__2].i = ctemp.i; // , expr subst
        /* L30: */
    }
    return 0;
}
/* clacrt_ */
