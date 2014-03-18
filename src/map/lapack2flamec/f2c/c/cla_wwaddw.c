/* ../netlib/cla_wwaddw.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief \b CLA_WWADDW adds a vector into a doubled-single vector. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CLA_WWADDW + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cla_wwa ddw.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cla_wwa ddw.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cla_wwa ddw.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CLA_WWADDW( N, X, Y, W ) */
/* .. Scalar Arguments .. */
/* INTEGER N */
/* .. */
/* .. Array Arguments .. */
/* COMPLEX X( * ), Y( * ), W( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CLA_WWADDW adds a vector W into a doubled-single vector (X, Y). */
/* > */
/* > This works for all extant IBM's hex and binary floating point */
/* > arithmetics, but not for decimal. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The length of vectors X, Y, and W. */
/* > \endverbatim */
/* > */
/* > \param[in,out] X */
/* > \verbatim */
/* > X is COMPLEX array, dimension (N) */
/* > The first part of the doubled-single accumulation vector. */
/* > \endverbatim */
/* > */
/* > \param[in,out] Y */
/* > \verbatim */
/* > Y is COMPLEX array, dimension (N) */
/* > The second part of the doubled-single accumulation vector. */
/* > \endverbatim */
/* > */
/* > \param[in] W */
/* > \verbatim */
/* > W is COMPLEX array, dimension (N) */
/* > The vector to be added. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date September 2012 */
/* > \ingroup complexOTHERcomputational */
/* ===================================================================== */
/* Subroutine */
int cla_wwaddw_(integer *n, complex *x, complex *y, complex *w)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5;
    complex q__1, q__2, q__3;
    /* Local variables */
    integer i__;
    complex s;
    /* -- LAPACK computational routine (version 3.4.2) -- */
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
    --w;
    --y;
    --x;
    /* Function Body */
    i__1 = *n;
    for (i__ = 1;
            i__ <= i__1;
            ++i__)
    {
        i__2 = i__;
        i__3 = i__;
        q__1.r = x[i__2].r + w[i__3].r;
        q__1.i = x[i__2].i + w[i__3].i; // , expr subst
        s.r = q__1.r;
        s.i = q__1.i; // , expr subst
        q__2.r = s.r + s.r;
        q__2.i = s.i + s.i; // , expr subst
        q__1.r = q__2.r - s.r;
        q__1.i = q__2.i - s.i; // , expr subst
        s.r = q__1.r;
        s.i = q__1.i; // , expr subst
        i__2 = i__;
        i__3 = i__;
        q__3.r = x[i__3].r - s.r;
        q__3.i = x[i__3].i - s.i; // , expr subst
        i__4 = i__;
        q__2.r = q__3.r + w[i__4].r;
        q__2.i = q__3.i + w[i__4].i; // , expr subst
        i__5 = i__;
        q__1.r = q__2.r + y[i__5].r;
        q__1.i = q__2.i + y[i__5].i; // , expr subst
        y[i__2].r = q__1.r;
        y[i__2].i = q__1.i; // , expr subst
        i__2 = i__;
        x[i__2].r = s.r;
        x[i__2].i = s.i; // , expr subst
        /* L10: */
    }
    return 0;
}
/* cla_wwaddw__ */
