/* ../netlib/izmax1.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief \b IZMAX1 finds the index of the vector element whose real part has maximum absolute value. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download IZMAX1 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/izmax1. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/izmax1. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/izmax1. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* INTEGER FUNCTION IZMAX1( N, CX, INCX ) */
/* .. Scalar Arguments .. */
/* INTEGER INCX, N */
/* .. */
/* .. Array Arguments .. */
/* COMPLEX*16 CX( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > IZMAX1 finds the index of the element whose real part has maximum */
/* > absolute value. */
/* > */
/* > Based on IZAMAX from Level 1 BLAS. */
/* > The change is to use the 'genuine' absolute value. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The number of elements in the vector CX. */
/* > \endverbatim */
/* > */
/* > \param[in] CX */
/* > \verbatim */
/* > CX is COMPLEX*16 array, dimension (N) */
/* > The vector whose elements will be summed. */
/* > \endverbatim */
/* > */
/* > \param[in] INCX */
/* > \verbatim */
/* > INCX is INTEGER */
/* > The spacing between successive values of CX. INCX >= 1. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date September 2012 */
/* > \ingroup complex16OTHERauxiliary */
/* > \par Contributors: */
/* ================== */
/* > */
/* > Nick Higham for use with ZLACON. */
/* ===================================================================== */
integer izmax1_(integer *n, doublecomplex *cx, integer *incx)
{
    /* System generated locals */
    integer ret_val, i__1;
    /* Builtin functions */
    double z_abs(doublecomplex *);
    /* Local variables */
    integer i__, ix;
    doublereal smax;
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
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Statement Functions .. */
    /* .. */
    /* .. Statement Function definitions .. */
    /* NEXT LINE IS THE ONLY MODIFICATION. */
    /* .. */
    /* .. Executable Statements .. */
    /* Parameter adjustments */
    --cx;
    /* Function Body */
    ret_val = 0;
    if (*n < 1)
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
        goto L30;
    }
    /* CODE FOR INCREMENT NOT EQUAL TO 1 */
    ix = 1;
    smax = z_abs(&cx[1]);
    ix += *incx;
    i__1 = *n;
    for (i__ = 2;
            i__ <= i__1;
            ++i__)
    {
        if (z_abs(&cx[ix]) <= smax)
        {
            goto L10;
        }
        ret_val = i__;
        smax = z_abs(&cx[ix]);
L10:
        ix += *incx;
        /* L20: */
    }
    return ret_val;
    /* CODE FOR INCREMENT EQUAL TO 1 */
L30:
    smax = z_abs(&cx[1]);
    i__1 = *n;
    for (i__ = 2;
            i__ <= i__1;
            ++i__)
    {
        if (z_abs(&cx[i__]) <= smax)
        {
            goto L40;
        }
        ret_val = i__;
        smax = z_abs(&cx[i__]);
L40:
        ;
    }
    return ret_val;
    /* End of IZMAX1 */
}
/* izmax1_ */
