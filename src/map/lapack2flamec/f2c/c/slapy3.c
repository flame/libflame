/* slapy3.f -- translated by f2c (version 20190311). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief \b SLAPY3 returns sqrt(x2+y2+z2). */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download SLAPY3 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slapy3. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slapy3. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slapy3. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* REAL FUNCTION SLAPY3( X, Y, Z ) */
/* .. Scalar Arguments .. */
/* REAL X, Y, Z */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLAPY3 returns sqrt(x**2+y**2+z**2), taking care not to cause */
/* > unnecessary overflow and unnecessary underflow. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] X */
/* > \verbatim */
/* > X is REAL */
/* > \endverbatim */
/* > */
/* > \param[in] Y */
/* > \verbatim */
/* > Y is REAL */
/* > \endverbatim */
/* > */
/* > \param[in] Z */
/* > \verbatim */
/* > Z is REAL */
/* > X, Y and Z specify the values x, y and z. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \ingroup OTHERauxiliary */
/* ===================================================================== */
real slapy3_(real *x, real *y, real *z__)
{
    AOCL_DTL_TRACE_ENTRY_INDENT
    /* System generated locals */
    real ret_val, r__1, r__2, r__3;
    /* Builtin functions */
    double sqrt(doublereal);
    /* Local variables */
    real w, xabs, yabs, zabs;
    extern real slamch_(char *);
    real hugeval;
    /* -- LAPACK auxiliary routine -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* .. Scalar Arguments .. */
    /* .. */
    /* ===================================================================== */
    /* .. Parameters .. */
    /* .. */
    /* .. Local Scalars .. */
    /* .. */
    /* .. External Subroutines .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    hugeval = slamch_("Overflow");
    xabs = f2c_abs(*x);
    yabs = f2c_abs(*y);
    zabs = f2c_abs(*z__);
    /* Computing MAX */
    r__1 = fla_max(xabs,yabs);
    w = fla_max(r__1,zabs);
    if (w == 0.f || w > hugeval)
    {
        /* W can be zero for fla_max(0,nan,0) */
        /* adding all three entries together will make sure */
        /* NaN will not disappear. */
        ret_val = xabs + yabs + zabs;
    }
    else
    {
        /* Computing 2nd power */
        r__1 = xabs / w;
        /* Computing 2nd power */
        r__2 = yabs / w;
        /* Computing 2nd power */
        r__3 = zabs / w;
        ret_val = w * sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3);
    }
    AOCL_DTL_TRACE_EXIT_INDENT
    return ret_val;
    /* End of SLAPY3 */
}
/* slapy3_ */
