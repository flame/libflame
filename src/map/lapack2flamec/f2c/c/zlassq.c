/* ../netlib/zlassq.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief \b ZLASSQ updates a sum of squares represented in scaled form. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download ZLASSQ + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlassq. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlassq. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlassq. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE ZLASSQ( N, X, INCX, SCALE, SUMSQ ) */
/* .. Scalar Arguments .. */
/* INTEGER INCX, N */
/* DOUBLE PRECISION SCALE, SUMSQ */
/* .. */
/* .. Array Arguments .. */
/* COMPLEX*16 X( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZLASSQ returns the values scl and ssq such that */
/* > */
/* > ( scl**2 )*ssq = x( 1 )**2 +...+ x( n )**2 + ( scale**2 )*sumsq, */
/* > */
/* > where x( i ) = f2c_dabs( X( 1 + ( i - 1 )*INCX ) ). The value of sumsq is */
/* > assumed to be at least unity and the value of ssq will then satisfy */
/* > */
/* > 1.0 .le. ssq .le. ( sumsq + 2*n ). */
/* > */
/* > scale is assumed to be non-negative and scl returns the value */
/* > */
/* > scl = fla_max( scale, f2c_dabs( real( x( i ) ) ), f2c_dabs( aimag( x( i ) ) ) ), */
/* > i */
/* > */
/* > scale and sumsq must be supplied in SCALE and SUMSQ respectively. */
/* > SCALE and SUMSQ are overwritten by scl and ssq respectively. */
/* > */
/* > The routine makes only one pass through the vector X. */
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
/* > X is COMPLEX*16 array, dimension (N) */
/* > The vector x as described above. */
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
/* > SCALE is DOUBLE PRECISION */
/* > On entry, the value scale in the equation above. */
/* > On exit, SCALE is overwritten with the value scl . */
/* > \endverbatim */
/* > */
/* > \param[in,out] SUMSQ */
/* > \verbatim */
/* > SUMSQ is DOUBLE PRECISION */
/* > On entry, the value sumsq in the equation above. */
/* > On exit, SUMSQ is overwritten with the value ssq . */
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
/* Subroutine */
int zlassq_(integer *n, doublecomplex *x, integer *incx, doublereal *scl, doublereal *sumsq)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("zlassq inputs: n %" FLA_IS ", incx %" FLA_IS ", scl %lf, sumsq %lf", *n, *incx, *scl, *sumsq);
    /* System generated locals */
    integer i__1, i__2;
    doublereal r__1, r__2;
    /* Builtin functions */
    double pow_ri(doublereal *, integer *), d_imag(doublecomplex *), sqrt(doublereal);
    /* Local variables */
    extern logical disnan_(doublereal *);
    integer i__;
    doublereal ax;
    integer ix;
    doublereal sbi, abig, amed, sbig, tbig, asml, ymin, ssml, tsml, ymax;
    logical notbig;
    /* ...Translated by Pacific-Sierra Research vf90 Personal 3.4N3 09:17:33 8/30/21 */
    /* ...Switches: */
    /* use LA_CONSTANTS, & */
    /* only: wp=>sp, zero=>szero, one=>sone, & */
    /* sbig=>ssbig, ssml=>sssml, tbig=>stbig, tsml=>stsml */
    /* use LA_XISNAN */
    /* .. Scalar Arguments .. */
    /* .. */
    /* .. Array Arguments .. */
    /* .. */
    /* .. Local Scalars .. */
    /* Parameter adjustments */
    --x;
    /* Function Body */
    tsml = 1.4916681462400413E-154;
    tbig = 1.9979190722022350E+146;
    ssml = 4.4989137945431964E+161;
    sbig = 1.1113793747425387E-162;
    sbi = 0.;
    /* .. */
    /* Quick return if possible */
    if (disnan_(scl) || disnan_(sumsq)) {
    AOCL_DTL_TRACE_LOG_EXIT
        return 0;
    }
    if (*sumsq == 0.) {
        *scl = 1.;
    }
    if (*scl == 0.) {
        *scl = 1.;
        *sumsq = 0.;
    }
    if (*n <= 0) {
    AOCL_DTL_TRACE_LOG_EXIT
        return 0;
    }
    /* Compute the sum of squares in 3 accumulators: */
    /* abig -- sums of squares scaled down to avoid overflow */
    /* asml -- sums of squares scaled up to avoid underflow */
    /* amed -- sums of squares that do not require scaling */
    /* The thresholds and multipliers are */
    /* tbig -- values bigger than this are scaled down by sbig */
    /* tsml -- values smaller than this are scaled up by ssml */
    notbig = TRUE_;
    asml = 0.;
    amed = 0.;
    abig = 0.;
    ix = 1;
    if (*incx < 0) {
        ix = 1 - (*n - 1) * *incx;
    }
    i__1 = *n;
    for (i__ = 1;
            i__ <= i__1;
            ++i__) {
        i__2 = ix;
        ax = (r__1 = x[i__2].r, f2c_abs(r__1));
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
        ax = (r__1 = d_imag(&x[ix]), f2c_abs(r__1));
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
    if (abig > 0.) {
        if (amed > 0. || disnan_(&amed)) {
            abig += amed * sbig * sbi;
        }
        *scl = 1. / sbig;
        *sumsq = abig;
    }
    else if (asml > 0.) {
        /* Combine amed and asml if asml > 0. */
        if (amed > 0. || disnan_(&amed)) {
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
            *scl = 1.;
            /* Computing 2nd power */
            r__1 = ymax;
            /* Computing 2nd power */
            r__2 = ymin / ymax;
            *sumsq = r__1 * r__1 * (r__2 * r__2 + 1.);
        }
        else {
            *scl = 1. / ssml;
            *sumsq = asml;
        }
    }
    else {
        /* Otherwise all values are mid-range or zero */
        *scl = 1.;
        *sumsq = amed;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return 0;
}
/* zlassq_ */