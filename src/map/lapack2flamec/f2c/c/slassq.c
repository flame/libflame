/* ../netlib/slassq.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief \b SLASSQ updates a sum of squares represented in scaled form. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download SLASSQ + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slassq. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slassq. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slassq. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE SLASSQ( N, X, INCX, SCALE, SUMSQ ) */
/* .. Scalar Arguments .. */
/* INTEGER INCX, N */
/* REAL SCALE, SUMSQ */
/* .. */
/* .. Array Arguments .. */
/* REAL X( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLASSQ returns the values scl and smsq such that */
/* > */
/* > ( scl**2 )*smsq = x( 1 )**2 +...+ x( n )**2 + ( scale**2 )*sumsq, */
/* > */
/* > where x( i ) = X( 1 + ( i - 1 )*INCX ). The value of sumsq is */
/* > assumed to be non-negative and scl returns the value */
/* > */
/* > scl = fla_max( scale, f2c_abs( x( i ) ) ). */
/* > */
/* > scale and sumsq must be supplied in SCALE and SUMSQ and */
/* > scl and smsq are overwritten on SCALE and SUMSQ respectively. */
/* > */
/* > The routine makes only one pass through the vector x. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The number of elements to be used from the vector X. */
/* > \endverbatim */
/* > */
/* > \param[in] X */
/* > \verbatim */
/* > X is REAL array, dimension (N) */
/* > The vector for which a scaled sum of squares is computed. */
/* > x( i ) = X( 1 + ( i - 1 )*INCX ), 1 <= i <= n. */
/* > \endverbatim */
/* > */
/* > \param[in] INCX */
/* > \verbatim */
/* > INCX is INTEGER */
/* > The increment between successive values of the vector X. */
/* > INCX > 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] SCALE */
/* > \verbatim */
/* > SCALE is REAL */
/* > On entry, the value scale in the equation above. */
/* > On exit, SCALE is overwritten with scl , the scaling factor */
/* > for the sum of squares. */
/* > \endverbatim */
/* > */
/* > \param[in,out] SUMSQ */
/* > \verbatim */
/* > SUMSQ is REAL */
/* > On entry, the value sumsq in the equation above. */
/* > On exit, SUMSQ is overwritten with smsq , the basic sum of */
/* > squares from which scl has been factored out. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date September 2012 */
/* > \ingroup auxOTHERauxiliary */
/* ===================================================================== */
/* Subroutine */
int slassq_(integer *n, real *x, integer *incx, real *scl, real *sumsq)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("slassq inputs: n %" FLA_IS ", incx %" FLA_IS "",*n, *incx);
    /* System generated locals */
    integer i__1, i__2;
    real r__1, r__2;
    /* Builtin functions */
    double pow_ri(real *, real *), sqrt(doublereal);
    integer i__;
    real ax;
    integer ix;
    real abig, amed, sbig, tbig, asml, ymin, ssml, tsml, ymax;
    logical notbig;
    /* ...Translated by Pacific-Sierra Research vf90 Personal 3.4N3 09:17:33 8/30/21 */
    /* .. Scalar Arguments .. */
    /* .. */
    /* .. Array Arguments .. */
    /* .. */
    /* .. Local Scalars .. */
    /* Parameter adjustments */
    --x;
    /* Function Body */
    sbig = 1.32348898E-23;
    ssml = 3.77789319E+22;
    tsml = 1.08420217E-19;
    tbig = 4.50359963E+15;
    /* .. */
    /* Quick return if possible */
    if (scl != scl || sumsq != sumsq) {
        return 0;
    }
    if (*sumsq == 0.f) {
        *scl = 1.f;
    }
    if (*scl == 0.f) {
        *scl = 1.f;
        *sumsq = 0.f;
    }
    if (*n <= 0) {
        return 0;
    }
    /* Compute the sum of squares in 3 accumulators: */
    /* abig -- sums of squares scaled down to avoid overflow */
    /* asml -- sums of squares scaled up to avoid underflow */
    /* amed -- sums of squares that do not require scaling */
    /* The thresholds and multipliers are */
    /* dtbig -- values bigger than this are scaled down by dsbig */
    /* dtsml -- values smaller than this are scaled up by dssml */
    notbig = TRUE_;
    asml = 0.f;
    amed = 0.f;
    abig = 0.f;
    ix = 1;
    if (*incx < 0) {
        ix = 1 - (*n - 1) * *incx;
    }
    i__1 = *n;
    for (i__ = 1;
            i__ <= i__1;
            ++i__) {
        ax = f2c_dabs(x[ix]);
        if (ax > tbig) {
            /* Computing 2nd power */
            r__1 = ax * sbig;
            abig += r__1 * r__1;
            notbig = FALSE_;
        }
        else if (ax < tsml) {
            if (notbig) {
                /* Computing 2nd power */
                r__1 = ax * ssml;
                asml += r__1 * r__1;
            }
        }
        else {
            /* Computing 2nd power */
            r__1 = ax;
            amed += r__1 * r__1;
        }
        ix += *incx;
    }
    /* Put the existing sum of squares into one of the accumulators */
    if (*sumsq > 0.f) {
        ax = *scl * sqrt(*sumsq);
        if (ax > tbig) {
            /* Computing 2nd power */
            r__1 = *scl * sbig;
            abig += (r__1 * r__1)* *sumsq;
            notbig = FALSE_;
        }
        else if (ax < tsml) {
            if (notbig) {
                /* Computing 2nd power */
                r__1 = *scl * ssml;
                asml += (r__1 * r__1)* *sumsq;
            }
        }
        else {
            /* Computing 2nd power */
            r__1 = *scl;
            amed += (r__1 * r__1)* *sumsq;
        }
    }
    /* Combine abig and amed or amed and asml if more than one */
    /* accumulator was used. */
    if (abig > 0.f) {
        if (amed > 0.f || amed != amed) {
            abig += amed * sbig * sbig;
        }
        *scl = 1.f / sbig;
        *sumsq = abig;
    }
    else if (asml > 0.f) {
        /* Combine amed and asml if asml > 0. */
        if (amed > 0.f || amed != amed) {
            amed = sqrt(amed);
            asml = sqrt(asml) / ssml;
            if (asml > amed) {
                ymin = amed;
                ymax = asml;
            }
            else {
                ymin = asml;
                ymax = amed;
            }
            *scl = 1.f;
            /* Computing 2nd power */
            r__1 = ymax;
            /* Computing 2nd power */
            r__2 = ymin / ymax;
            *sumsq = r__1 * r__1 * (r__2 * r__2 + 1.f);
        }
        else {
            *scl = 1.f / ssml;
            *sumsq = asml;
        }
    }
    else {
        /* Otherwise all values are mid-range or zero */
        *scl = 1.f;
        *sumsq = amed;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return 0;
}
/* slassq_ */