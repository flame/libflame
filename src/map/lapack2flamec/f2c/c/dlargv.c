/* ../netlib/dlargv.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief \b DLARGV generates a vector of plane rotations with real cosines and real sines. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download DLARGV + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlargv. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlargv. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlargv. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE DLARGV( N, X, INCX, Y, INCY, C, INCC ) */
/* .. Scalar Arguments .. */
/* INTEGER INCC, INCX, INCY, N */
/* .. */
/* .. Array Arguments .. */
/* DOUBLE PRECISION C( * ), X( * ), Y( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLARGV generates a vector of real plane rotations, determined by */
/* > elements of the real vectors x and y. For i = 1,2,...,n */
/* > */
/* > ( c(i) s(i) ) ( x(i) ) = ( a(i) ) */
/* > ( -s(i) c(i) ) ( y(i) ) = ( 0 ) */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The number of plane rotations to be generated. */
/* > \endverbatim */
/* > */
/* > \param[in,out] X */
/* > \verbatim */
/* > X is DOUBLE PRECISION array, */
/* > dimension (1+(N-1)*INCX) */
/* > On entry, the vector x. */
/* > On exit, x(i) is overwritten by a(i), for i = 1,...,n. */
/* > \endverbatim */
/* > */
/* > \param[in] INCX */
/* > \verbatim */
/* > INCX is INTEGER */
/* > The increment between elements of X. INCX > 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] Y */
/* > \verbatim */
/* > Y is DOUBLE PRECISION array, */
/* > dimension (1+(N-1)*INCY) */
/* > On entry, the vector y. */
/* > On exit, the sines of the plane rotations. */
/* > \endverbatim */
/* > */
/* > \param[in] INCY */
/* > \verbatim */
/* > INCY is INTEGER */
/* > The increment between elements of Y. INCY > 0. */
/* > \endverbatim */
/* > */
/* > \param[out] C */
/* > \verbatim */
/* > C is DOUBLE PRECISION array, dimension (1+(N-1)*INCC) */
/* > The cosines of the plane rotations. */
/* > \endverbatim */
/* > */
/* > \param[in] INCC */
/* > \verbatim */
/* > INCC is INTEGER */
/* > The increment between elements of C. INCC > 0. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date September 2012 */
/* > \ingroup doubleOTHERauxiliary */
/* ===================================================================== */
/* Subroutine */
int dlargv_(integer *n, doublereal *x, integer *incx, doublereal *y, integer *incy, doublereal *c__, integer *incc)
{
    /* System generated locals */
    integer i__1;
    /* Builtin functions */
    double sqrt(doublereal);
    /* Local variables */
    doublereal f, g;
    integer i__;
    doublereal t;
    integer ic, ix, iy;
    doublereal tt;
    /* -- LAPACK auxiliary routine (version 3.4.2) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* September 2012 */
    /* .. Scalar Arguments .. */
    /* .. */
    /* .. Array Arguments .. */
    /* .. */
    /* ===================================================================== */
    /* .. Parameters .. */
    /* .. */
    /* .. Local Scalars .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Parameter adjustments */
    --c__;
    --y;
    --x;
    /* Function Body */
    ix = 1;
    iy = 1;
    ic = 1;
    i__1 = *n;
    for (i__ = 1;
            i__ <= i__1;
            ++i__)
    {
        f = x[ix];
        g = y[iy];
        if (g == 0.)
        {
            c__[ic] = 1.;
        }
        else if (f == 0.)
        {
            c__[ic] = 0.;
            y[iy] = 1.;
            x[ix] = g;
        }
        else if (f2c_abs(f) > f2c_abs(g))
        {
            t = g / f;
            tt = sqrt(t * t + 1.);
            c__[ic] = 1. / tt;
            y[iy] = t * c__[ic];
            x[ix] = f * tt;
        }
        else
        {
            t = f / g;
            tt = sqrt(t * t + 1.);
            y[iy] = 1. / tt;
            c__[ic] = t * y[iy];
            x[ix] = g * tt;
        }
        ic += *incc;
        iy += *incy;
        ix += *incx;
        /* L10: */
    }
    return 0;
    /* End of DLARGV */
}
/* dlargv_ */
