/* sroundup_lwork.f -- translated by f2c (version 20190311). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief \b SROUNDUP_LWORK */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* Definition: */
/* =========== */
/* REAL FUNCTION SROUNDUP_LWORK( LWORK ) */
/* .. Scalar Arguments .. */
/* INTEGER LWORK */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > SROUNDUP_LWORK deals with a subtle bug with returning LWORK as a Float. */
/* > This routine guarantees it is rounded up instead of down by */
/* > multiplying LWORK by 1+eps when it is necessary, where eps is the relative machine precision. */
/* > E.g., */
/* > */
/* > float( 16777217 ) == 16777216 */
/* > float( 16777217 ) * (1.+eps) == 16777218 */
/* > */
/* > \return SROUNDUP_LWORK */
/* > \verbatim */
/* > SROUNDUP_LWORK >= LWORK. */
/* > SROUNDUP_LWORK is guaranteed to have zero decimal part. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] LWORK Workspace size. */
/* Authors: */
/* ======== */
/* > \author Weslley Pereira, University of Colorado Denver, USA */
/* > \ingroup auxOTHERauxiliary */
/* > \par Further Details: */
/* ===================== */
/* > */
/* > \verbatim */
/* > This routine was inspired in the method `magma_zmake_lwork` from MAGMA. */
/* > \see https://bitbucket.org/icl/magma/src/master/control/magma_zauxiliary.cpp */
/* > \endverbatim */
/* ===================================================================== */
real sroundup_lwork(integer *lwork)
{
    /* System generated locals */
    real ret_val, eps;
    //external function
    extern real slamch_(char *);
    /* -- LAPACK auxiliary routine -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* .. Scalar Arguments .. */
    /* .. */
    /* ===================================================================== */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* .. */
    ret_val = (real) (*lwork);
    eps = slamch_("Epsilon");
    if ((integer) ret_val < *lwork)
    {
        /* Force round up of LWORK */
        ret_val *= (1. + eps);
    }
    return ret_val;
    /* End of SROUNDUP_LWORK */
}
/* sroundup_lwork__ */
