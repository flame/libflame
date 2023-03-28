/* ../netlib/dladiv.f -- translated by f2c (version 20160102). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief \b DLADIV performs complex division in real arithmetic, avoiding unnecessary overflow. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download DLADIV + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dladiv. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dladiv. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dladiv. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE DLADIV( A, B, C, D, P, Q ) */
/* .. Scalar Arguments .. */
/* DOUBLE PRECISION A, B, C, D, P, Q */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLADIV performs complex division in real arithmetic */
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
/* > A is DOUBLE PRECISION */
/* > \endverbatim */
/* > */
/* > \param[in] B */
/* > \verbatim */
/* > B is DOUBLE PRECISION */
/* > \endverbatim */
/* > */
/* > \param[in] C */
/* > \verbatim */
/* > C is DOUBLE PRECISION */
/* > \endverbatim */
/* > */
/* > \param[in] D */
/* > \verbatim */
/* > D is DOUBLE PRECISION */
/* > The scalars a, b, c, and d in the above expression. */
/* > \endverbatim */
/* > */
/* > \param[out] P */
/* > \verbatim */
/* > P is DOUBLE PRECISION */
/* > \endverbatim */
/* > */
/* > \param[out] Q */
/* > \verbatim */
/* > Q is DOUBLE PRECISION */
/* > The scalars p and q in the above expression. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date January 2013 */
/* > \ingroup doubleOTHERauxiliary */
/* ===================================================================== */
/* Subroutine */
int dladiv_(doublereal *a, doublereal *b, doublereal *c__, doublereal *d__, doublereal *p, doublereal *q)
{
    AOCL_DTL_TRACE_ENTRY_INDENT
    /* System generated locals */
    doublereal d__1, d__2;
    /* Local variables */
    static TLS_CLASS_SPEC doublereal be, un, ov, eps;
    doublereal s, aa, ab, bb, cc, cd, dd;
    extern doublereal dlamch_(char *);
    extern /* Subroutine */
    int dladiv1_(doublereal *, doublereal *, doublereal *, doublereal *, doublereal *, doublereal *);
    /* -- LAPACK auxiliary routine (version 3.7.0) -- */
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
    d__1 = f2c_abs(*a);
    d__2 = f2c_abs(*b); // , expr subst
    ab = fla_max(d__1,d__2);
    /* Computing MAX */
    d__1 = f2c_abs(*c__);
    d__2 = f2c_abs(*d__); // , expr subst
    cd = fla_max(d__1,d__2);

    if(ov == 0)
    {
        ov = dlamch_("Overflow threshold") * .5;
        eps = dlamch_("Epsilon");
        un = dlamch_("Safe minimum") * 2. / eps;
        be = 2. / (eps * eps);
    }

    s = 1.;
    if (ab >= ov )
    {
        aa *= .5;
        bb *= .5;
        s *= 2.;
    }
    if (cd >= ov )
    {
        cc *= .5;
        dd *= .5;
        s *= .5;
    }
    if (ab <= un )
    {
        aa *= be;
        bb *= be;
        s /= be;
    }
    if (cd <= un )
    {
        cc *= be;
        dd *= be;
        s *= be;
    }

    if (f2c_abs(*d__) <= f2c_abs(*c__))
    {
        dladiv1_(&aa, &bb, &cc, &dd, p, q);
    }
    else
    {
        dladiv1_(&bb, &aa, &dd, &cc, p, q);
        *q = -(*q);
    }
    *p *= s;
    *q *= s;
    AOCL_DTL_TRACE_EXIT_INDENT
    return 0;
    /* End of DLADIV */
}
/* dladiv_ */
/* > \ingroup doubleOTHERauxiliary */
/* Subroutine */
int dladiv1_(doublereal *a, doublereal *b, doublereal *c__, doublereal *d__, doublereal *p, doublereal *q)
{
    AOCL_DTL_TRACE_ENTRY_INDENT
    doublereal r__, t;
    extern doublereal dladiv2_(doublereal *, doublereal *, doublereal *, doublereal *, doublereal *, doublereal *);
    /* -- LAPACK auxiliary routine (version 3.7.0) -- */
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
    t = 1. / (*c__ + *d__ * r__);
    *p = dladiv2_(a, b, c__, d__, &r__, &t);
    *a = -(*a);
    *q = dladiv2_(b, a, c__, d__, &r__, &t);
    AOCL_DTL_TRACE_EXIT_INDENT
    return 0;
    /* End of DLADIV1 */
}
/* dladiv1_ */
/* > \ingroup doubleOTHERauxiliary */
doublereal dladiv2_(doublereal *a, doublereal *b, doublereal *c__, doublereal *d__, doublereal *r__, doublereal *t)
{
    AOCL_DTL_TRACE_ENTRY_INDENT
    /* System generated locals */
    doublereal ret_val;
    /* Local variables */
    doublereal br;
    /* -- LAPACK auxiliary routine (version 3.7.0) -- */
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
    if (*r__ != 0.)
    {
        br = *b * *r__;
        if (br != 0.)
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
    AOCL_DTL_TRACE_EXIT_INDENT
    return ret_val;
    /* End of DLADIV12 */
}
/* dladiv2_ */
