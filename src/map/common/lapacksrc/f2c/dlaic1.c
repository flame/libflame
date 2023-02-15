/* ../netlib/dlaic1.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
static doublereal c_b5 = 1.;
/* > \brief \b DLAIC1 applies one step of incremental condition estimation. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download DLAIC1 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaic1. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaic1. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaic1. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE DLAIC1( JOB, J, X, SEST, W, GAMMA, SESTPR, S, C ) */
/* .. Scalar Arguments .. */
/* INTEGER J, JOB */
/* DOUBLE PRECISION C, GAMMA, S, SEST, SESTPR */
/* .. */
/* .. Array Arguments .. */
/* DOUBLE PRECISION W( J ), X( J ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLAIC1 applies one step of incremental condition estimation in */
/* > its simplest version: */
/* > */
/* > Let x, twonorm(x) = 1, be an approximate singular vector of an j-by-j */
/* > lower triangular matrix L, such that */
/* > twonorm(L*x) = sest */
/* > Then DLAIC1 computes sestpr, s, c such that */
/* > the vector */
/* > [ s*x ] */
/* > xhat = [ c ] */
/* > is an approximate singular vector of */
/* > [ L 0 ] */
/* > Lhat = [ w**T gamma ] */
/* > in the sense that */
/* > twonorm(Lhat*xhat) = sestpr. */
/* > */
/* > Depending on JOB, an estimate for the largest or smallest singular */
/* > value is computed. */
/* > */
/* > Note that [s c]**T and sestpr**2 is an eigenpair of the system */
/* > */
/* > diag(sest*sest, 0) + [alpha gamma] * [ alpha ] */
/* > [ gamma ] */
/* > */
/* > where alpha = x**T*w. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] JOB */
/* > \verbatim */
/* > JOB is INTEGER */
/* > = 1: an estimate for the largest singular value is computed. */
/* > = 2: an estimate for the smallest singular value is computed. */
/* > \endverbatim */
/* > */
/* > \param[in] J */
/* > \verbatim */
/* > J is INTEGER */
/* > Length of X and W */
/* > \endverbatim */
/* > */
/* > \param[in] X */
/* > \verbatim */
/* > X is DOUBLE PRECISION array, dimension (J) */
/* > The j-vector x. */
/* > \endverbatim */
/* > */
/* > \param[in] SEST */
/* > \verbatim */
/* > SEST is DOUBLE PRECISION */
/* > Estimated singular value of j by j matrix L */
/* > \endverbatim */
/* > */
/* > \param[in] W */
/* > \verbatim */
/* > W is DOUBLE PRECISION array, dimension (J) */
/* > The j-vector w. */
/* > \endverbatim */
/* > */
/* > \param[in] GAMMA */
/* > \verbatim */
/* > GAMMA is DOUBLE PRECISION */
/* > The diagonal element gamma. */
/* > \endverbatim */
/* > */
/* > \param[out] SESTPR */
/* > \verbatim */
/* > SESTPR is DOUBLE PRECISION */
/* > Estimated singular value of (j+1) by (j+1) matrix Lhat. */
/* > \endverbatim */
/* > */
/* > \param[out] S */
/* > \verbatim */
/* > S is DOUBLE PRECISION */
/* > Sine needed in forming xhat. */
/* > \endverbatim */
/* > */
/* > \param[out] C */
/* > \verbatim */
/* > C is DOUBLE PRECISION */
/* > Cosine needed in forming xhat. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date September 2012 */
/* > \ingroup doubleOTHERauxiliary */
/* ===================================================================== */
/* Subroutine */
int dlaic1_(integer *job, integer *j, doublereal *x, doublereal *sest, doublereal *w, doublereal *gamma, doublereal * sestpr, doublereal *s, doublereal *c__)
{
    /* System generated locals */
    doublereal d__1, d__2, d__3, d__4;
    /* Builtin functions */
    double sqrt(doublereal), d_sign(doublereal *, doublereal *);
    /* Local variables */
    doublereal b, t, s1, s2, eps, tmp;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, integer *);
    doublereal sine, test, zeta1, zeta2, alpha, norma;
    extern doublereal dlamch_(char *);
    doublereal absgam, absalp, cosine, absest;
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
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. External Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Parameter adjustments */
    --w;
    --x;
    /* Function Body */
    eps = dlamch_("Epsilon");
    alpha = ddot_(j, &x[1], &c__1, &w[1], &c__1);
    absalp = f2c_abs(alpha);
    absgam = f2c_abs(*gamma);
    absest = f2c_abs(*sest);
    if (*job == 1)
    {
        /* Estimating largest singular value */
        /* special cases */
        if (*sest == 0.)
        {
            s1 = max(absgam,absalp);
            if (s1 == 0.)
            {
                *s = 0.;
                *c__ = 1.;
                *sestpr = 0.;
            }
            else
            {
                *s = alpha / s1;
                *c__ = *gamma / s1;
                tmp = sqrt(*s * *s + *c__ * *c__);
                *s /= tmp;
                *c__ /= tmp;
                *sestpr = s1 * tmp;
            }
            return 0;
        }
        else if (absgam <= eps * absest)
        {
            *s = 1.;
            *c__ = 0.;
            tmp = max(absest,absalp);
            s1 = absest / tmp;
            s2 = absalp / tmp;
            *sestpr = tmp * sqrt(s1 * s1 + s2 * s2);
            return 0;
        }
        else if (absalp <= eps * absest)
        {
            s1 = absgam;
            s2 = absest;
            if (s1 <= s2)
            {
                *s = 1.;
                *c__ = 0.;
                *sestpr = s2;
            }
            else
            {
                *s = 0.;
                *c__ = 1.;
                *sestpr = s1;
            }
            return 0;
        }
        else if (absest <= eps * absalp || absest <= eps * absgam)
        {
            s1 = absgam;
            s2 = absalp;
            if (s1 <= s2)
            {
                tmp = s1 / s2;
                *s = sqrt(tmp * tmp + 1.);
                *sestpr = s2 * *s;
                *c__ = *gamma / s2 / *s;
                *s = d_sign(&c_b5, &alpha) / *s;
            }
            else
            {
                tmp = s2 / s1;
                *c__ = sqrt(tmp * tmp + 1.);
                *sestpr = s1 * *c__;
                *s = alpha / s1 / *c__;
                *c__ = d_sign(&c_b5, gamma) / *c__;
            }
            return 0;
        }
        else
        {
            /* normal case */
            zeta1 = alpha / absest;
            zeta2 = *gamma / absest;
            b = (1. - zeta1 * zeta1 - zeta2 * zeta2) * .5;
            *c__ = zeta1 * zeta1;
            if (b > 0.)
            {
                t = *c__ / (b + sqrt(b * b + *c__));
            }
            else
            {
                t = sqrt(b * b + *c__) - b;
            }
            sine = -zeta1 / t;
            cosine = -zeta2 / (t + 1.);
            tmp = sqrt(sine * sine + cosine * cosine);
            *s = sine / tmp;
            *c__ = cosine / tmp;
            *sestpr = sqrt(t + 1.) * absest;
            return 0;
        }
    }
    else if (*job == 2)
    {
        /* Estimating smallest singular value */
        /* special cases */
        if (*sest == 0.)
        {
            *sestpr = 0.;
            if (max(absgam,absalp) == 0.)
            {
                sine = 1.;
                cosine = 0.;
            }
            else
            {
                sine = -(*gamma);
                cosine = alpha;
            }
            /* Computing MAX */
            d__1 = f2c_abs(sine);
            d__2 = f2c_abs(cosine); // , expr subst
            s1 = max(d__1,d__2);
            *s = sine / s1;
            *c__ = cosine / s1;
            tmp = sqrt(*s * *s + *c__ * *c__);
            *s /= tmp;
            *c__ /= tmp;
            return 0;
        }
        else if (absgam <= eps * absest)
        {
            *s = 0.;
            *c__ = 1.;
            *sestpr = absgam;
            return 0;
        }
        else if (absalp <= eps * absest)
        {
            s1 = absgam;
            s2 = absest;
            if (s1 <= s2)
            {
                *s = 0.;
                *c__ = 1.;
                *sestpr = s1;
            }
            else
            {
                *s = 1.;
                *c__ = 0.;
                *sestpr = s2;
            }
            return 0;
        }
        else if (absest <= eps * absalp || absest <= eps * absgam)
        {
            s1 = absgam;
            s2 = absalp;
            if (s1 <= s2)
            {
                tmp = s1 / s2;
                *c__ = sqrt(tmp * tmp + 1.);
                *sestpr = absest * (tmp / *c__);
                *s = -(*gamma / s2) / *c__;
                *c__ = d_sign(&c_b5, &alpha) / *c__;
            }
            else
            {
                tmp = s2 / s1;
                *s = sqrt(tmp * tmp + 1.);
                *sestpr = absest / *s;
                *c__ = alpha / s1 / *s;
                *s = -d_sign(&c_b5, gamma) / *s;
            }
            return 0;
        }
        else
        {
            /* normal case */
            zeta1 = alpha / absest;
            zeta2 = *gamma / absest;
            /* Computing MAX */
            d__3 = zeta1 * zeta1 + 1. + (d__1 = zeta1 * zeta2, f2c_abs(d__1));
            d__4 = (d__2 = zeta1 * zeta2, f2c_abs(d__2)) + zeta2 * zeta2; // , expr subst
            norma = max(d__3,d__4);
            /* See if root is closer to zero or to ONE */
            test = (zeta1 - zeta2) * 2. * (zeta1 + zeta2) + 1.;
            if (test >= 0.)
            {
                /* root is close to zero, compute directly */
                b = (zeta1 * zeta1 + zeta2 * zeta2 + 1.) * .5;
                *c__ = zeta2 * zeta2;
                t = *c__ / (b + sqrt((d__1 = b * b - *c__, f2c_abs(d__1))));
                sine = zeta1 / (1. - t);
                cosine = -zeta2 / t;
                *sestpr = sqrt(t + eps * 4. * eps * norma) * absest;
            }
            else
            {
                /* root is closer to ONE, shift by that amount */
                b = (zeta2 * zeta2 + zeta1 * zeta1 - 1.) * .5;
                *c__ = zeta1 * zeta1;
                if (b >= 0.)
                {
                    t = -(*c__) / (b + sqrt(b * b + *c__));
                }
                else
                {
                    t = b - sqrt(b * b + *c__);
                }
                sine = -zeta1 / t;
                cosine = -zeta2 / (t + 1.);
                *sestpr = sqrt(t + 1. + eps * 4. * eps * norma) * absest;
            }
            tmp = sqrt(sine * sine + cosine * cosine);
            *s = sine / tmp;
            *c__ = cosine / tmp;
            return 0;
        }
    }
    return 0;
    /* End of DLAIC1 */
}
/* dlaic1_ */
