/*
    Copyright (c) 2019-2023 Advanced Micro Devices, Inc.
*/
/* ../netlib/zladiv.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief \b ZLADIV performs complex division in real arithmetic, avoiding unnecessary overflow. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download ZLADIV + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zladiv. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zladiv. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zladiv. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* COMPLEX*16 FUNCTION ZLADIV( X, Y ) */
/* .. Scalar Arguments .. */
/* COMPLEX*16 X, Y */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZLADIV := X / Y, where X and Y are complex. The computation of X / Y */
/* > will not overflow on an intermediary step unless the results */
/* > overflows. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] X */
/* > \verbatim */
/* > X is COMPLEX*16 */
/* > \endverbatim */
/* > */
/* > \param[in] Y */
/* > \verbatim */
/* > Y is COMPLEX*16 */
/* > The complex scalars X and Y. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date September 2012 */
/* > \ingroup complex16OTHERauxiliary */
/* ===================================================================== */
/* Double Complex */
#ifdef FLA_ENABLE_VOID_RETURN_COMPLEX_FUNCTION
VOID zladiv_(doublecomplex * ret_val, doublecomplex *x, doublecomplex *y)
{
    AOCL_DTL_TRACE_ENTRY_INDENT
    /* System generated locals */
    doublereal d__1, d__2, d__3, d__4;
    doublecomplex z__1;
    /* Local variables */
    doublereal zi, zr;
    extern /* Subroutine */
    int dladiv_(doublereal *, doublereal *, doublereal *, doublereal *, doublereal *, doublereal *);
    /* -- LAPACK auxiliary routine (version 3.4.2) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* September 2012 */
    /* .. Scalar Arguments .. */
    /* .. */
    /* ===================================================================== */
    /* .. Local Scalars .. */
    /* .. */
    /* .. External Subroutines .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    d__1 = x->r;
    d__2 = x->i;
    d__3 = y->r;
    d__4 = y->i;
    dladiv_(&d__1, &d__2, &d__3, &d__4, &zr, &zi);
    z__1.r = zr;
    z__1.i = zi; // , expr subst
    ret_val->r = z__1.r, ret_val->i = z__1.i;
    AOCL_DTL_TRACE_EXIT_INDENT
    return ;
    /* End of ZLADIV */
}
#else
doublecomplex zladiv_(doublecomplex *x, doublecomplex *y)
{
    AOCL_DTL_TRACE_ENTRY_INDENT
    /* System generated locals */
    doublereal d__1, d__2, d__3, d__4;
    doublecomplex z__1;
    /* Local variables */
    doublereal zi, zr;
    extern /* Subroutine */
    int dladiv_(doublereal *, doublereal *, doublereal *, doublereal *, doublereal *, doublereal *);
    /* -- LAPACK auxiliary routine (version 3.4.2) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* September 2012 */
    /* .. Scalar Arguments .. */
    /* .. */
    /* ===================================================================== */
    /* .. Local Scalars .. */
    /* .. */
    /* .. External Subroutines .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    d__1 = x->r;
    d__2 = x->i;
    d__3 = y->r;
    d__4 = y->i;
    dladiv_(&d__1, &d__2, &d__3, &d__4, &zr, &zi);
    z__1.r = zr;
    z__1.i = zi; // , expr subst
    //ret_val->r = z__1.r, ret_val->i = z__1.i;
    AOCL_DTL_TRACE_EXIT_INDENT
    return z__1;
    /* End of ZLADIV */
}

#endif

void zladiv_f2c_(doublecomplex *ret_val, doublecomplex *x, doublecomplex *y)
{

#ifdef FLA_ENABLE_VOID_RETURN_COMPLEX_FUNCTION
    zladiv_(ret_val, x, y);
#else
    *ret_val = zladiv_(x, y);
#endif

    return;
}

/* zladiv_ */
