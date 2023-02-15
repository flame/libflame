/* ../netlib/sladiv.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief \b SLADIV performs complex division in real arithmetic, avoiding unnecessary overflow. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download SLADIV + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sladiv. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sladiv. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sladiv. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE SLADIV( A, B, C, D, P, Q ) */
/* .. Scalar Arguments .. */
/* REAL A, B, C, D, P, Q */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLADIV performs complex division in real arithmetic */
/* > */
/* > a + i*b */
/* > p + i*q = --------- */
/* > c + i*d */
/* > */
/* > The algorithm is due to Michael Baudin and Robert L. Smith */
/* > and can be found in the paper */
/* > "A Robust Complex Division in Scilab" */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] A */
/* > \verbatim */
/* > A is REAL */
/* > \endverbatim */
/* > */
/* > \param[in] B */
/* > \verbatim */
/* > B is REAL */
/* > \endverbatim */
/* > */
/* > \param[in] C */
/* > \verbatim */
/* > C is REAL */
/* > \endverbatim */
/* > */
/* > \param[in] D */
/* > \verbatim */
/* > D is REAL */
/* > The scalars a, b, c, and d in the above expression. */
/* > \endverbatim */
/* > */
/* > \param[out] P */
/* > \verbatim */
/* > P is REAL */
/* > \endverbatim */
/* > */
/* > \param[out] Q */
/* > \verbatim */
/* > Q is REAL */
/* > The scalars p and q in the above expression. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date January 2013 */
/* > \ingroup auxOTHERauxiliary */
/* ===================================================================== */
/* Subroutine */
int sladiv_(real *a, real *b, real *c__, real *d__, real *p, real *q)
{
    /* System generated locals */
    real r__1, r__2;
    /* Local variables */
    real s, aa, ab, bb, cc, cd, dd, be, un, ov, eps;
    extern real slamch_(char *);
    extern /* Subroutine */
    int sladiv1_(real *, real *, real *, real *, real *, real *);
    /* -- LAPACK auxiliary routine (version 3.5.0) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* January 2013 */
    /* .. Scalar Arguments .. */
    /* .. */
    /* ===================================================================== */
    /* .. Parameters .. */
    /* .. Local Scalars .. */
    /* .. */
    /* .. External Functions .. */
    /* .. */
    /* .. External Subroutines .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    aa = *a;
    bb = *b;
    cc = *c__;
    dd = *d__;
    /* Computing MAX */
    r__1 = f2c_abs(*a);
    r__2 = f2c_abs(*b); // , expr subst
    ab = max(r__1,r__2);
    /* Computing MAX */
    r__1 = f2c_abs(*c__);
    r__2 = f2c_abs(*d__); // , expr subst
    cd = max(r__1,r__2);
    s = 1.f;
    ov = slamch_("Overflow threshold");
    un = slamch_("Safe minimum");
    eps = slamch_("Epsilon");
    be = 2.f / (eps * eps);
    if (ab >= ov * .5f)
    {
        aa *= .5f;
        bb *= .5f;
        s *= 2.f;
    }
    if (cd >= ov * .5f)
    {
        cc *= .5f;
        dd *= .5f;
        s *= .5f;
    }
    if (ab <= un * 2.f / eps)
    {
        aa *= be;
        bb *= be;
        s /= be;
    }
    if (cd <= un * 2.f / eps)
    {
        cc *= be;
        dd *= be;
        s *= be;
    }
    if (f2c_abs(*d__) <= f2c_abs(*c__))
    {
        sladiv1_(&aa, &bb, &cc, &dd, p, q);
    }
    else
    {
        sladiv1_(&bb, &aa, &dd, &cc, p, q);
        *q = -(*q);
    }
    *p *= s;
    *q *= s;
    return 0;
    /* End of SLADIV */
}
/* sladiv_ */
/* Subroutine */
int sladiv1_(real *a, real *b, real *c__, real *d__, real *p, real *q)
{
    real r__, t;
    extern real sladiv2_(real *, real *, real *, real *, real *, real *);
    /* -- LAPACK auxiliary routine (version 3.5.0) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* January 2013 */
    /* .. Scalar Arguments .. */
    /* .. */
    /* ===================================================================== */
    /* .. Parameters .. */
    /* .. Local Scalars .. */
    /* .. */
    /* .. External Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    r__ = *d__ / *c__;
    t = 1.f / (*c__ + *d__ * r__);
    *p = sladiv2_(a, b, c__, d__, &r__, &t);
    *a = -(*a);
    *q = sladiv2_(b, a, c__, d__, &r__, &t);
    return 0;
    /* End of SLADIV1 */
}
/* sladiv1_ */
real sladiv2_(real *a, real *b, real *c__, real *d__, real *r__, real *t)
{
    /* System generated locals */
    real ret_val;
    /* Local variables */
    real br;
    /* -- LAPACK auxiliary routine (version 3.5.0) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* January 2013 */
    /* .. Scalar Arguments .. */
    /* .. */
    /* ===================================================================== */
    /* .. Parameters .. */
    /* .. Local Scalars .. */
    /* .. */
    /* .. Executable Statements .. */
    if (*r__ != 0.f)
    {
        br = *b * *r__;
        if (br != 0.f)
        {
            ret_val = (*a + br) * *t;
        }
        else
        {
            ret_val = *a * *t + *b * *t * *r__;
        }
    }
    else
    {
        ret_val = (*a + *d__ * (*b / *c__)) * *t;
    }
    return ret_val;
    /* End of SLADIV */
}
/* sladiv2_ */
