/* slapy2.f -- translated by f2c (version 20190311). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief \b SLAPY2 returns sqrt(x2+y2). */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download SLAPY2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slapy2. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slapy2. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slapy2. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* REAL FUNCTION SLAPY2( X, Y ) */
/* .. Scalar Arguments .. */
/* REAL X, Y */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLAPY2 returns sqrt(x**2+y**2), taking care not to cause unnecessary */
/* > overflow and unnecessary underflow. */
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
/* > X and Y specify the values x and y. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \ingroup OTHERauxiliary */
/* ===================================================================== */
real slapy2_(real *x, real *y)
{
    AOCL_DTL_TRACE_ENTRY_INDENT
    /* System generated locals */
    real ret_val, r__1;
    /* Builtin functions */
    double sqrt(doublereal);
    /* Local variables */
    logical x_is_nan__, y_is_nan__;
    real w, z__, xabs, yabs;
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
    /* .. External Functions .. */
    /* .. */
    /* .. External Subroutines .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    x_is_nan__ = (*x != *x);
    y_is_nan__ = (*y != *y);
    if (x_is_nan__)
    {
        ret_val = *x;
    }
    if (y_is_nan__)
    {
        ret_val = *y;
    }
    hugeval = slamch_("Overflow");
    if (! (x_is_nan__ || y_is_nan__))
    {
        xabs = f2c_abs(*x);
        yabs = f2c_abs(*y);
        w = fla_max(xabs,yabs);
        z__ = fla_min(xabs,yabs);
        if (z__ == 0.f || w > hugeval)
        {
            ret_val = w;
        }
        else
        {
            /* Computing 2nd power */
            r__1 = z__ / w;
            ret_val = w * sqrt(r__1 * r__1 + 1.f);
        }
    }
    AOCL_DTL_TRACE_EXIT_INDENT
    return ret_val;
    /* End of SLAPY2 */
}
/* slapy2_ */
