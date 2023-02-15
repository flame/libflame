/* ../netlib/zlaic1.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
/* > \brief \b ZLAIC1 applies one step of incremental condition estimation. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download ZLAIC1 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlaic1. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlaic1. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlaic1. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE ZLAIC1( JOB, J, X, SEST, W, GAMMA, SESTPR, S, C ) */
/* .. Scalar Arguments .. */
/* INTEGER J, JOB */
/* DOUBLE PRECISION SEST, SESTPR */
/* COMPLEX*16 C, GAMMA, S */
/* .. */
/* .. Array Arguments .. */
/* COMPLEX*16 W( J ), X( J ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZLAIC1 applies one step of incremental condition estimation in */
/* > its simplest version: */
/* > */
/* > Let x, twonorm(x) = 1, be an approximate singular vector of an j-by-j */
/* > lower triangular matrix L, such that */
/* > twonorm(L*x) = sest */
/* > Then ZLAIC1 computes sestpr, s, c such that */
/* > the vector */
/* > [ s*x ] */
/* > xhat = [ c ] */
/* > is an approximate singular vector of */
/* > [ L 0 ] */
/* > Lhat = [ w**H gamma ] */
/* > in the sense that */
/* > twonorm(Lhat*xhat) = sestpr. */
/* > */
/* > Depending on JOB, an estimate for the largest or smallest singular */
/* > value is computed. */
/* > */
/* > Note that [s c]**H and sestpr**2 is an eigenpair of the system */
/* > */
/* > diag(sest*sest, 0) + [alpha gamma] * [ conjg(alpha) ] */
/* > [ conjg(gamma) ] */
/* > */
/* > where alpha = x**H * w. */
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
/* > X is COMPLEX*16 array, dimension (J) */
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
/* > W is COMPLEX*16 array, dimension (J) */
/* > The j-vector w. */
/* > \endverbatim */
/* > */
/* > \param[in] GAMMA */
/* > \verbatim */
/* > GAMMA is COMPLEX*16 */
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
/* > S is COMPLEX*16 */
/* > Sine needed in forming xhat. */
/* > \endverbatim */
/* > */
/* > \param[out] C */
/* > \verbatim */
/* > C is COMPLEX*16 */
/* > Cosine needed in forming xhat. */
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
int zlaic1_(integer *job, integer *j, doublecomplex *x, doublereal *sest, doublecomplex *w, doublecomplex *gamma, doublereal * sestpr, doublecomplex *s, doublecomplex *c__)
{
    /* System generated locals */
    doublereal d__1, d__2;
    doublecomplex z__1, z__2, z__3, z__4, z__5, z__6;
    /* Builtin functions */
    double z_abs(doublecomplex *);
    void d_cnjg(doublecomplex *, doublecomplex *), z_sqrt(doublecomplex *, doublecomplex *);
    double sqrt(doublereal);
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *);
    /* Local variables */
    doublereal b, t, s1, s2, scl, eps, tmp;
    doublecomplex sine;
    doublereal test, zeta1, zeta2;
    doublecomplex alpha;
    doublereal norma;
    extern /* Double Complex */
    VOID zdotc_f2c_(doublecomplex *, integer *, doublecomplex *, integer *, doublecomplex *, integer *);
    extern doublereal dlamch_(char *);
    doublereal absgam, absalp;
    doublecomplex cosine;
    doublereal absest;
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
    zdotc_f2c_(&z__1, j, &x[1], &c__1, &w[1], &c__1);
    alpha.r = z__1.r;
    alpha.i = z__1.i; // , expr subst
    absalp = z_abs(&alpha);
    absgam = z_abs(gamma);
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
                s->r = 0., s->i = 0.;
                c__->r = 1., c__->i = 0.;
                *sestpr = 0.;
            }
            else
            {
                z__1.r = alpha.r / s1;
                z__1.i = alpha.i / s1; // , expr subst
                s->r = z__1.r, s->i = z__1.i;
                z__1.r = gamma->r / s1;
                z__1.i = gamma->i / s1; // , expr subst
                c__->r = z__1.r, c__->i = z__1.i;
                d_cnjg(&z__4, s);
                z__3.r = s->r * z__4.r - s->i * z__4.i;
                z__3.i = s->r * z__4.i + s->i * z__4.r; // , expr subst
                d_cnjg(&z__6, c__);
                z__5.r = c__->r * z__6.r - c__->i * z__6.i;
                z__5.i = c__->r * z__6.i + c__->i * z__6.r; // , expr subst
                z__2.r = z__3.r + z__5.r;
                z__2.i = z__3.i + z__5.i; // , expr subst
                z_sqrt(&z__1, &z__2);
                tmp = z__1.r;
                z__1.r = s->r / tmp;
                z__1.i = s->i / tmp; // , expr subst
                s->r = z__1.r, s->i = z__1.i;
                z__1.r = c__->r / tmp;
                z__1.i = c__->i / tmp; // , expr subst
                c__->r = z__1.r, c__->i = z__1.i;
                *sestpr = s1 * tmp;
            }
            return 0;
        }
        else if (absgam <= eps * absest)
        {
            s->r = 1., s->i = 0.;
            c__->r = 0., c__->i = 0.;
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
                s->r = 1., s->i = 0.;
                c__->r = 0., c__->i = 0.;
                *sestpr = s2;
            }
            else
            {
                s->r = 0., s->i = 0.;
                c__->r = 1., c__->i = 0.;
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
                scl = sqrt(tmp * tmp + 1.);
                *sestpr = s2 * scl;
                z__2.r = alpha.r / s2;
                z__2.i = alpha.i / s2; // , expr subst
                z__1.r = z__2.r / scl;
                z__1.i = z__2.i / scl; // , expr subst
                s->r = z__1.r, s->i = z__1.i;
                z__2.r = gamma->r / s2;
                z__2.i = gamma->i / s2; // , expr subst
                z__1.r = z__2.r / scl;
                z__1.i = z__2.i / scl; // , expr subst
                c__->r = z__1.r, c__->i = z__1.i;
            }
            else
            {
                tmp = s2 / s1;
                scl = sqrt(tmp * tmp + 1.);
                *sestpr = s1 * scl;
                z__2.r = alpha.r / s1;
                z__2.i = alpha.i / s1; // , expr subst
                z__1.r = z__2.r / scl;
                z__1.i = z__2.i / scl; // , expr subst
                s->r = z__1.r, s->i = z__1.i;
                z__2.r = gamma->r / s1;
                z__2.i = gamma->i / s1; // , expr subst
                z__1.r = z__2.r / scl;
                z__1.i = z__2.i / scl; // , expr subst
                c__->r = z__1.r, c__->i = z__1.i;
            }
            return 0;
        }
        else
        {
            /* normal case */
            zeta1 = absalp / absest;
            zeta2 = absgam / absest;
            b = (1. - zeta1 * zeta1 - zeta2 * zeta2) * .5;
            d__1 = zeta1 * zeta1;
            c__->r = d__1, c__->i = 0.;
            if (b > 0.)
            {
                d__1 = b * b;
                z__4.r = d__1 + c__->r;
                z__4.i = c__->i; // , expr subst
                z_sqrt(&z__3, &z__4);
                z__2.r = b + z__3.r;
                z__2.i = z__3.i; // , expr subst
                z_div(&z__1, c__, &z__2);
                t = z__1.r;
            }
            else
            {
                d__1 = b * b;
                z__3.r = d__1 + c__->r;
                z__3.i = c__->i; // , expr subst
                z_sqrt(&z__2, &z__3);
                z__1.r = z__2.r - b;
                z__1.i = z__2.i; // , expr subst
                t = z__1.r;
            }
            z__3.r = alpha.r / absest;
            z__3.i = alpha.i / absest; // , expr subst
            z__2.r = -z__3.r;
            z__2.i = -z__3.i; // , expr subst
            z__1.r = z__2.r / t;
            z__1.i = z__2.i / t; // , expr subst
            sine.r = z__1.r;
            sine.i = z__1.i; // , expr subst
            z__3.r = gamma->r / absest;
            z__3.i = gamma->i / absest; // , expr subst
            z__2.r = -z__3.r;
            z__2.i = -z__3.i; // , expr subst
            d__1 = t + 1.;
            z__1.r = z__2.r / d__1;
            z__1.i = z__2.i / d__1; // , expr subst
            cosine.r = z__1.r;
            cosine.i = z__1.i; // , expr subst
            d_cnjg(&z__4, &sine);
            z__3.r = sine.r * z__4.r - sine.i * z__4.i;
            z__3.i = sine.r * z__4.i + sine.i * z__4.r; // , expr subst
            d_cnjg(&z__6, &cosine);
            z__5.r = cosine.r * z__6.r - cosine.i * z__6.i;
            z__5.i = cosine.r * z__6.i + cosine.i * z__6.r; // , expr subst
            z__2.r = z__3.r + z__5.r;
            z__2.i = z__3.i + z__5.i; // , expr subst
            z_sqrt(&z__1, &z__2);
            tmp = z__1.r;
            z__1.r = sine.r / tmp;
            z__1.i = sine.i / tmp; // , expr subst
            s->r = z__1.r, s->i = z__1.i;
            z__1.r = cosine.r / tmp;
            z__1.i = cosine.i / tmp; // , expr subst
            c__->r = z__1.r, c__->i = z__1.i;
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
                sine.r = 1.;
                sine.i = 0.; // , expr subst
                cosine.r = 0.;
                cosine.i = 0.; // , expr subst
            }
            else
            {
                d_cnjg(&z__2, gamma);
                z__1.r = -z__2.r;
                z__1.i = -z__2.i; // , expr subst
                sine.r = z__1.r;
                sine.i = z__1.i; // , expr subst
                d_cnjg(&z__1, &alpha);
                cosine.r = z__1.r;
                cosine.i = z__1.i; // , expr subst
            }
            /* Computing MAX */
            d__1 = z_abs(&sine);
            d__2 = z_abs(&cosine); // , expr subst
            s1 = max(d__1,d__2);
            z__1.r = sine.r / s1;
            z__1.i = sine.i / s1; // , expr subst
            s->r = z__1.r, s->i = z__1.i;
            z__1.r = cosine.r / s1;
            z__1.i = cosine.i / s1; // , expr subst
            c__->r = z__1.r, c__->i = z__1.i;
            d_cnjg(&z__4, s);
            z__3.r = s->r * z__4.r - s->i * z__4.i;
            z__3.i = s->r * z__4.i + s->i * z__4.r; // , expr subst
            d_cnjg(&z__6, c__);
            z__5.r = c__->r * z__6.r - c__->i * z__6.i;
            z__5.i = c__->r * z__6.i + c__->i * z__6.r; // , expr subst
            z__2.r = z__3.r + z__5.r;
            z__2.i = z__3.i + z__5.i; // , expr subst
            z_sqrt(&z__1, &z__2);
            tmp = z__1.r;
            z__1.r = s->r / tmp;
            z__1.i = s->i / tmp; // , expr subst
            s->r = z__1.r, s->i = z__1.i;
            z__1.r = c__->r / tmp;
            z__1.i = c__->i / tmp; // , expr subst
            c__->r = z__1.r, c__->i = z__1.i;
            return 0;
        }
        else if (absgam <= eps * absest)
        {
            s->r = 0., s->i = 0.;
            c__->r = 1., c__->i = 0.;
            *sestpr = absgam;
            return 0;
        }
        else if (absalp <= eps * absest)
        {
            s1 = absgam;
            s2 = absest;
            if (s1 <= s2)
            {
                s->r = 0., s->i = 0.;
                c__->r = 1., c__->i = 0.;
                *sestpr = s1;
            }
            else
            {
                s->r = 1., s->i = 0.;
                c__->r = 0., c__->i = 0.;
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
                scl = sqrt(tmp * tmp + 1.);
                *sestpr = absest * (tmp / scl);
                d_cnjg(&z__4, gamma);
                z__3.r = z__4.r / s2;
                z__3.i = z__4.i / s2; // , expr subst
                z__2.r = -z__3.r;
                z__2.i = -z__3.i; // , expr subst
                z__1.r = z__2.r / scl;
                z__1.i = z__2.i / scl; // , expr subst
                s->r = z__1.r, s->i = z__1.i;
                d_cnjg(&z__3, &alpha);
                z__2.r = z__3.r / s2;
                z__2.i = z__3.i / s2; // , expr subst
                z__1.r = z__2.r / scl;
                z__1.i = z__2.i / scl; // , expr subst
                c__->r = z__1.r, c__->i = z__1.i;
            }
            else
            {
                tmp = s2 / s1;
                scl = sqrt(tmp * tmp + 1.);
                *sestpr = absest / scl;
                d_cnjg(&z__4, gamma);
                z__3.r = z__4.r / s1;
                z__3.i = z__4.i / s1; // , expr subst
                z__2.r = -z__3.r;
                z__2.i = -z__3.i; // , expr subst
                z__1.r = z__2.r / scl;
                z__1.i = z__2.i / scl; // , expr subst
                s->r = z__1.r, s->i = z__1.i;
                d_cnjg(&z__3, &alpha);
                z__2.r = z__3.r / s1;
                z__2.i = z__3.i / s1; // , expr subst
                z__1.r = z__2.r / scl;
                z__1.i = z__2.i / scl; // , expr subst
                c__->r = z__1.r, c__->i = z__1.i;
            }
            return 0;
        }
        else
        {
            /* normal case */
            zeta1 = absalp / absest;
            zeta2 = absgam / absest;
            /* Computing MAX */
            d__1 = zeta1 * zeta1 + 1. + zeta1 * zeta2;
            d__2 = zeta1 * zeta2 + zeta2 * zeta2; // , expr subst
            norma = max(d__1,d__2);
            /* See if root is closer to zero or to ONE */
            test = (zeta1 - zeta2) * 2. * (zeta1 + zeta2) + 1.;
            if (test >= 0.)
            {
                /* root is close to zero, compute directly */
                b = (zeta1 * zeta1 + zeta2 * zeta2 + 1.) * .5;
                d__1 = zeta2 * zeta2;
                c__->r = d__1, c__->i = 0.;
                d__2 = b * b;
                z__2.r = d__2 - c__->r;
                z__2.i = -c__->i; // , expr subst
                d__1 = b + sqrt(z_abs(&z__2));
                z__1.r = c__->r / d__1;
                z__1.i = c__->i / d__1; // , expr subst
                t = z__1.r;
                z__2.r = alpha.r / absest;
                z__2.i = alpha.i / absest; // , expr subst
                d__1 = 1. - t;
                z__1.r = z__2.r / d__1;
                z__1.i = z__2.i / d__1; // , expr subst
                sine.r = z__1.r;
                sine.i = z__1.i; // , expr subst
                z__3.r = gamma->r / absest;
                z__3.i = gamma->i / absest; // , expr subst
                z__2.r = -z__3.r;
                z__2.i = -z__3.i; // , expr subst
                z__1.r = z__2.r / t;
                z__1.i = z__2.i / t; // , expr subst
                cosine.r = z__1.r;
                cosine.i = z__1.i; // , expr subst
                *sestpr = sqrt(t + eps * 4. * eps * norma) * absest;
            }
            else
            {
                /* root is closer to ONE, shift by that amount */
                b = (zeta2 * zeta2 + zeta1 * zeta1 - 1.) * .5;
                d__1 = zeta1 * zeta1;
                c__->r = d__1, c__->i = 0.;
                if (b >= 0.)
                {
                    z__2.r = -c__->r;
                    z__2.i = -c__->i; // , expr subst
                    d__1 = b * b;
                    z__5.r = d__1 + c__->r;
                    z__5.i = c__->i; // , expr subst
                    z_sqrt(&z__4, &z__5);
                    z__3.r = b + z__4.r;
                    z__3.i = z__4.i; // , expr subst
                    z_div(&z__1, &z__2, &z__3);
                    t = z__1.r;
                }
                else
                {
                    d__1 = b * b;
                    z__3.r = d__1 + c__->r;
                    z__3.i = c__->i; // , expr subst
                    z_sqrt(&z__2, &z__3);
                    z__1.r = b - z__2.r;
                    z__1.i = -z__2.i; // , expr subst
                    t = z__1.r;
                }
                z__3.r = alpha.r / absest;
                z__3.i = alpha.i / absest; // , expr subst
                z__2.r = -z__3.r;
                z__2.i = -z__3.i; // , expr subst
                z__1.r = z__2.r / t;
                z__1.i = z__2.i / t; // , expr subst
                sine.r = z__1.r;
                sine.i = z__1.i; // , expr subst
                z__3.r = gamma->r / absest;
                z__3.i = gamma->i / absest; // , expr subst
                z__2.r = -z__3.r;
                z__2.i = -z__3.i; // , expr subst
                d__1 = t + 1.;
                z__1.r = z__2.r / d__1;
                z__1.i = z__2.i / d__1; // , expr subst
                cosine.r = z__1.r;
                cosine.i = z__1.i; // , expr subst
                *sestpr = sqrt(t + 1. + eps * 4. * eps * norma) * absest;
            }
            d_cnjg(&z__4, &sine);
            z__3.r = sine.r * z__4.r - sine.i * z__4.i;
            z__3.i = sine.r * z__4.i + sine.i * z__4.r; // , expr subst
            d_cnjg(&z__6, &cosine);
            z__5.r = cosine.r * z__6.r - cosine.i * z__6.i;
            z__5.i = cosine.r * z__6.i + cosine.i * z__6.r; // , expr subst
            z__2.r = z__3.r + z__5.r;
            z__2.i = z__3.i + z__5.i; // , expr subst
            z_sqrt(&z__1, &z__2);
            tmp = z__1.r;
            z__1.r = sine.r / tmp;
            z__1.i = sine.i / tmp; // , expr subst
            s->r = z__1.r, s->i = z__1.i;
            z__1.r = cosine.r / tmp;
            z__1.i = cosine.i / tmp; // , expr subst
            c__->r = z__1.r, c__->i = z__1.i;
            return 0;
        }
    }
    return 0;
    /* End of ZLAIC1 */
}
/* zlaic1_ */
