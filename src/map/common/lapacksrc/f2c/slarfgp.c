/* ../netlib/slarfgp.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief \b SLARFGP generates an elementary reflector (Householder matrix) with non-negatibe beta. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download SLARFGP + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slarfgp .f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slarfgp .f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slarfgp .f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE SLARFGP( N, ALPHA, X, INCX, TAU ) */
/* .. Scalar Arguments .. */
/* INTEGER INCX, N */
/* REAL ALPHA, TAU */
/* .. */
/* .. Array Arguments .. */
/* REAL X( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLARFGP generates a real elementary reflector H of order n, such */
/* > that */
/* > */
/* > H * ( alpha ) = ( beta ), H**T * H = I. */
/* > ( x ) ( 0 ) */
/* > */
/* > where alpha and beta are scalars, beta is non-negative, and x is */
/* > an (n-1)-element real vector. H is represented in the form */
/* > */
/* > H = I - tau * ( 1 ) * ( 1 v**T ) , */
/* > ( v ) */
/* > */
/* > where tau is a real scalar and v is a real (n-1)-element */
/* > vector. */
/* > */
/* > If the elements of x are all zero, then tau = 0 and H is taken to be */
/* > the unit matrix. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the elementary reflector. */
/* > \endverbatim */
/* > */
/* > \param[in,out] ALPHA */
/* > \verbatim */
/* > ALPHA is REAL */
/* > On entry, the value alpha. */
/* > On exit, it is overwritten with the value beta. */
/* > \endverbatim */
/* > */
/* > \param[in,out] X */
/* > \verbatim */
/* > X is REAL array, dimension */
/* > (1+(N-2)*f2c_abs(INCX)) */
/* > On entry, the vector x. */
/* > On exit, it is overwritten with the vector v. */
/* > \endverbatim */
/* > */
/* > \param[in] INCX */
/* > \verbatim */
/* > INCX is INTEGER */
/* > The increment between elements of X. INCX > 0. */
/* > \endverbatim */
/* > */
/* > \param[out] TAU */
/* > \verbatim */
/* > TAU is REAL */
/* > The value tau. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date September 2012 */
/* > \ingroup realOTHERauxiliary */
/* ===================================================================== */
/* Subroutine */
int slarfgp_(integer *n, real *alpha, real *x, integer *incx, real *tau)
{
    /* System generated locals */
    integer i__1;
    real r__1;
    /* Builtin functions */
    double r_sign(real *, real *);
    /* Local variables */
    integer j;
    real savealpha;
    integer knt;
    real beta;
    extern real snrm2_(integer *, real *, integer *);
    extern /* Subroutine */
    int sscal_(integer *, real *, real *, integer *);
    real xnorm;
    extern real slapy2_(real *, real *), slamch_(char *);
    real bignum, smlnum;
    /* -- LAPACK auxiliary routine (version 3.4.2) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* September 2012 */
    /* .. Scalar Arguments .. */
    /* .. */
    /* .. Array Arguments .. */
    /* .. */
    /* ===================================================================== */
    /* .. Parameters .. */
    /* .. */
    /* .. Local Scalars .. */
    /* .. */
    /* .. External Functions .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. External Subroutines .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Parameter adjustments */
    --x;
    /* Function Body */
    if (*n <= 0)
    {
        *tau = 0.f;
        return 0;
    }
    i__1 = *n - 1;
    xnorm = snrm2_(&i__1, &x[1], incx);
    if (xnorm == 0.f)
    {
        /* H = [+/-1, 0;
        I], sign chosen so ALPHA >= 0. */
        if (*alpha >= 0.f)
        {
            /* When TAU.eq.ZERO, the vector is special-cased to be */
            /* all zeros in the application routines. We do not need */
            /* to clear it. */
            *tau = 0.f;
        }
        else
        {
            /* However, the application routines rely on explicit */
            /* zero checks when TAU.ne.ZERO, and we must clear X. */
            *tau = 2.f;
            i__1 = *n - 1;
            for (j = 1;
                    j <= i__1;
                    ++j)
            {
                x[(j - 1) * *incx + 1] = 0.f;
            }
            *alpha = -(*alpha);
        }
    }
    else
    {
        /* general case */
        r__1 = slapy2_(alpha, &xnorm);
        beta = r_sign(&r__1, alpha);
        smlnum = slamch_("S") / slamch_("E");
        knt = 0;
        if (f2c_abs(beta) < smlnum)
        {
            /* XNORM, BETA may be inaccurate;
            scale X and recompute them */
            bignum = 1.f / smlnum;
L10:
            ++knt;
            i__1 = *n - 1;
            sscal_(&i__1, &bignum, &x[1], incx);
            beta *= bignum;
            *alpha *= bignum;
            if (f2c_abs(beta) < smlnum)
            {
                goto L10;
            }
            /* New BETA is at most 1, at least SMLNUM */
            i__1 = *n - 1;
            xnorm = snrm2_(&i__1, &x[1], incx);
            r__1 = slapy2_(alpha, &xnorm);
            beta = r_sign(&r__1, alpha);
        }
        savealpha = *alpha;
        *alpha += beta;
        if (beta < 0.f)
        {
            beta = -beta;
            *tau = -(*alpha) / beta;
        }
        else
        {
            *alpha = xnorm * (xnorm / *alpha);
            *tau = *alpha / beta;
            *alpha = -(*alpha);
        }
        if (f2c_abs(*tau) <= smlnum)
        {
            /* In the case where the computed TAU ends up being a denormalized number, */
            /* it loses relative accuracy. This is a BIG problem. Solution: flush TAU */
            /* to ZERO. This explains the next IF statement. */
            /* (Bug report provided by Pat Quillen from MathWorks on Jul 29, 2009.) */
            /* (Thanks Pat. Thanks MathWorks.) */
            if (savealpha >= 0.f)
            {
                *tau = 0.f;
            }
            else
            {
                *tau = 2.f;
                i__1 = *n - 1;
                for (j = 1;
                        j <= i__1;
                        ++j)
                {
                    x[(j - 1) * *incx + 1] = 0.f;
                }
                beta = -savealpha;
            }
        }
        else
        {
            /* This is the general case. */
            i__1 = *n - 1;
            r__1 = 1.f / *alpha;
            sscal_(&i__1, &r__1, &x[1], incx);
        }
        /* If BETA is subnormal, it may lose relative accuracy */
        i__1 = knt;
        for (j = 1;
                j <= i__1;
                ++j)
        {
            beta *= smlnum;
            /* L20: */
        }
        *alpha = beta;
    }
    return 0;
    /* End of SLARFGP */
}
/* slarfgp_ */
