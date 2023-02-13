/* ../netlib/dlaed6.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief \b DLAED6 used by sstedc. Computes one Newton step in solution of the secular equation. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download DLAED6 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaed6. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaed6. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaed6. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE DLAED6( KNITER, ORGATI, RHO, D, Z, FINIT, TAU, INFO ) */
/* .. Scalar Arguments .. */
/* LOGICAL ORGATI */
/* INTEGER INFO, KNITER */
/* DOUBLE PRECISION FINIT, RHO, TAU */
/* .. */
/* .. Array Arguments .. */
/* DOUBLE PRECISION D( 3 ), Z( 3 ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLAED6 computes the positive or negative root (closest to the origin) */
/* > of */
/* > z(1) z(2) z(3) */
/* > f(x) = rho + --------- + ---------- + --------- */
/* > d(1)-x d(2)-x d(3)-x */
/* > */
/* > It is assumed that */
/* > */
/* > if ORGATI = .true. the root is between d(2) and d(3);
*/
/* > otherwise it is between d(1) and d(2) */
/* > */
/* > This routine will be called by DLAED4 when necessary. In most cases, */
/* > the root sought is the smallest in magnitude, though it might not be */
/* > in some extremely rare situations. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] KNITER */
/* > \verbatim */
/* > KNITER is INTEGER */
/* > Refer to DLAED4 for its significance. */
/* > \endverbatim */
/* > */
/* > \param[in] ORGATI */
/* > \verbatim */
/* > ORGATI is LOGICAL */
/* > If ORGATI is true, the needed root is between d(2) and */
/* > d(3);
otherwise it is between d(1) and d(2). See */
/* > DLAED4 for further details. */
/* > \endverbatim */
/* > */
/* > \param[in] RHO */
/* > \verbatim */
/* > RHO is DOUBLE PRECISION */
/* > Refer to the equation f(x) above. */
/* > \endverbatim */
/* > */
/* > \param[in] D */
/* > \verbatim */
/* > D is DOUBLE PRECISION array, dimension (3) */
/* > D satisfies d(1) < d(2) < d(3). */
/* > \endverbatim */
/* > */
/* > \param[in] Z */
/* > \verbatim */
/* > Z is DOUBLE PRECISION array, dimension (3) */
/* > Each of the elements in z must be positive. */
/* > \endverbatim */
/* > */
/* > \param[in] FINIT */
/* > \verbatim */
/* > FINIT is DOUBLE PRECISION */
/* > The value of f at 0. It is more accurate than the one */
/* > evaluated inside this routine (if someone wants to do */
/* > so). */
/* > \endverbatim */
/* > */
/* > \param[out] TAU */
/* > \verbatim */
/* > TAU is DOUBLE PRECISION */
/* > The root of the equation f(x). */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit */
/* > > 0: if INFO = 1, failure to converge */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date September 2012 */
/* > \ingroup auxOTHERcomputational */
/* > \par Further Details: */
/* ===================== */
/* > */
/* > \verbatim */
/* > */
/* > 10/02/03: This version has a few statements commented out for thread */
/* > safety (machine parameters are computed on each entry). SJH. */
/* > */
/* > 05/10/06: Modified from a new version of Ren-Cang Li, use */
/* > Gragg-Thornton-Warner cubic convergent scheme for better stability. */
/* > \endverbatim */
/* > \par Contributors: */
/* ================== */
/* > */
/* > Ren-Cang Li, Computer Science Division, University of California */
/* > at Berkeley, USA */
/* > */
/* ===================================================================== */
/* Subroutine */
int dlaed6_(integer *kniter, logical *orgati, doublereal * rho, doublereal *d__, doublereal *z__, doublereal *finit, doublereal * tau, integer *info)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2, d__3, d__4;
    /* Builtin functions */
    double sqrt(doublereal), log(doublereal), pow_di(doublereal *, integer *);
    /* Local variables */
    doublereal a, b, c__, f;
    integer i__;
    doublereal fc, df, ddf, lbd, eta, ubd, eps, base;
    integer iter;
    doublereal temp, temp1, temp2, temp3, temp4;
    logical scale;
    integer niter;
    doublereal small1, small2, sminv1, sminv2;
    extern doublereal dlamch_(char *);
    doublereal dscale[3], sclfac, zscale[3], erretm, sclinv;
    /* -- LAPACK computational routine (version 3.4.2) -- */
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
    /* .. External Functions .. */
    /* .. */
    /* .. Local Arrays .. */
    /* .. */
    /* .. Local Scalars .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Parameter adjustments */
    --z__;
    --d__;
    /* Function Body */
    *info = 0;
    if (*orgati)
    {
        lbd = d__[2];
        ubd = d__[3];
    }
    else
    {
        lbd = d__[1];
        ubd = d__[2];
    }
    if (*finit < 0.)
    {
        lbd = 0.;
    }
    else
    {
        ubd = 0.;
    }
    niter = 1;
    *tau = 0.;
    if (*kniter == 2)
    {
        if (*orgati)
        {
            temp = (d__[3] - d__[2]) / 2.;
            c__ = *rho + z__[1] / (d__[1] - d__[2] - temp);
            a = c__ * (d__[2] + d__[3]) + z__[2] + z__[3];
            b = c__ * d__[2] * d__[3] + z__[2] * d__[3] + z__[3] * d__[2];
        }
        else
        {
            temp = (d__[1] - d__[2]) / 2.;
            c__ = *rho + z__[3] / (d__[3] - d__[2] - temp);
            a = c__ * (d__[1] + d__[2]) + z__[1] + z__[2];
            b = c__ * d__[1] * d__[2] + z__[1] * d__[2] + z__[2] * d__[1];
        }
        /* Computing MAX */
        d__1 = f2c_abs(a), d__2 = f2c_abs(b);
        d__1 = max(d__1,d__2);
        d__2 = f2c_abs(c__); // ; expr subst
        temp = max(d__1,d__2);
        a /= temp;
        b /= temp;
        c__ /= temp;
        if (c__ == 0.)
        {
            *tau = b / a;
        }
        else if (a <= 0.)
        {
            *tau = (a - sqrt((d__1 = a * a - b * 4. * c__, f2c_abs(d__1)))) / ( c__ * 2.);
        }
        else
        {
            *tau = b * 2. / (a + sqrt((d__1 = a * a - b * 4. * c__, f2c_abs(d__1)) ));
        }
        if (*tau < lbd || *tau > ubd)
        {
            *tau = (lbd + ubd) / 2.;
        }
        if (d__[1] == *tau || d__[2] == *tau || d__[3] == *tau)
        {
            *tau = 0.;
        }
        else
        {
            temp = *finit + *tau * z__[1] / (d__[1] * (d__[1] - *tau)) + *tau * z__[2] / (d__[2] * (d__[2] - *tau)) + *tau * z__[3] / ( d__[3] * (d__[3] - *tau));
            if (temp <= 0.)
            {
                lbd = *tau;
            }
            else
            {
                ubd = *tau;
            }
            if (f2c_abs(*finit) <= f2c_abs(temp))
            {
                *tau = 0.;
            }
        }
    }
    /* get machine parameters for possible scaling to avoid overflow */
    /* modified by Sven: parameters SMALL1, SMINV1, SMALL2, */
    /* SMINV2, EPS are not SAVEd anymore between one call to the */
    /* others but recomputed at each call */
    eps = dlamch_("Epsilon");
    base = dlamch_("Base");
    i__1 = (integer) (log(dlamch_("SafMin")) / log(base) / 3.);
    small1 = pow_di(&base, &i__1);
    sminv1 = 1. / small1;
    small2 = small1 * small1;
    sminv2 = sminv1 * sminv1;
    /* Determine if scaling of inputs necessary to avoid overflow */
    /* when computing 1/TEMP**3 */
    if (*orgati)
    {
        /* Computing MIN */
        d__3 = (d__1 = d__[2] - *tau, f2c_abs(d__1));
        d__4 = (d__2 = d__[3] - * tau, f2c_abs(d__2)); // , expr subst
        temp = min(d__3,d__4);
    }
    else
    {
        /* Computing MIN */
        d__3 = (d__1 = d__[1] - *tau, f2c_abs(d__1));
        d__4 = (d__2 = d__[2] - * tau, f2c_abs(d__2)); // , expr subst
        temp = min(d__3,d__4);
    }
    scale = FALSE_;
    if (temp <= small1)
    {
        scale = TRUE_;
        if (temp <= small2)
        {
            /* Scale up by power of radix nearest 1/SAFMIN**(2/3) */
            sclfac = sminv2;
            sclinv = small2;
        }
        else
        {
            /* Scale up by power of radix nearest 1/SAFMIN**(1/3) */
            sclfac = sminv1;
            sclinv = small1;
        }
        /* Scaling up safe because D, Z, TAU scaled elsewhere to be O(1) */
        for (i__ = 1;
                i__ <= 3;
                ++i__)
        {
            dscale[i__ - 1] = d__[i__] * sclfac;
            zscale[i__ - 1] = z__[i__] * sclfac;
            /* L10: */
        }
        *tau *= sclfac;
        lbd *= sclfac;
        ubd *= sclfac;
    }
    else
    {
        /* Copy D and Z to DSCALE and ZSCALE */
        for (i__ = 1;
                i__ <= 3;
                ++i__)
        {
            dscale[i__ - 1] = d__[i__];
            zscale[i__ - 1] = z__[i__];
            /* L20: */
        }
    }
    fc = 0.;
    df = 0.;
    ddf = 0.;
    for (i__ = 1;
            i__ <= 3;
            ++i__)
    {
        temp = 1. / (dscale[i__ - 1] - *tau);
        temp1 = zscale[i__ - 1] * temp;
        temp2 = temp1 * temp;
        temp3 = temp2 * temp;
        fc += temp1 / dscale[i__ - 1];
        df += temp2;
        ddf += temp3;
        /* L30: */
    }
    f = *finit + *tau * fc;
    if (f2c_abs(f) <= 0.)
    {
        goto L60;
    }
    if (f <= 0.)
    {
        lbd = *tau;
    }
    else
    {
        ubd = *tau;
    }
    /* Iteration begins -- Use Gragg-Thornton-Warner cubic convergent */
    /* scheme */
    /* It is not hard to see that */
    /* 1) Iterations will go up monotonically */
    /* if FINIT < 0;
    */
    /* 2) Iterations will go down monotonically */
    /* if FINIT > 0. */
    iter = niter + 1;
    for (niter = iter;
            niter <= 40;
            ++niter)
    {
        if (*orgati)
        {
            temp1 = dscale[1] - *tau;
            temp2 = dscale[2] - *tau;
        }
        else
        {
            temp1 = dscale[0] - *tau;
            temp2 = dscale[1] - *tau;
        }
        a = (temp1 + temp2) * f - temp1 * temp2 * df;
        b = temp1 * temp2 * f;
        c__ = f - (temp1 + temp2) * df + temp1 * temp2 * ddf;
        /* Computing MAX */
        d__1 = f2c_abs(a), d__2 = f2c_abs(b);
        d__1 = max(d__1,d__2);
        d__2 = f2c_abs(c__); // ; expr subst
        temp = max(d__1,d__2);
        a /= temp;
        b /= temp;
        c__ /= temp;
        if (c__ == 0.)
        {
            eta = b / a;
        }
        else if (a <= 0.)
        {
            eta = (a - sqrt((d__1 = a * a - b * 4. * c__, f2c_abs(d__1)))) / (c__ * 2.);
        }
        else
        {
            eta = b * 2. / (a + sqrt((d__1 = a * a - b * 4. * c__, f2c_abs(d__1))) );
        }
        if (f * eta >= 0.)
        {
            eta = -f / df;
        }
        *tau += eta;
        if (*tau < lbd || *tau > ubd)
        {
            *tau = (lbd + ubd) / 2.;
        }
        fc = 0.;
        erretm = 0.;
        df = 0.;
        ddf = 0.;
        for (i__ = 1;
                i__ <= 3;
                ++i__)
        {
            if (dscale[i__ - 1] - *tau != 0.)
            {
                temp = 1. / (dscale[i__ - 1] - *tau);
                temp1 = zscale[i__ - 1] * temp;
                temp2 = temp1 * temp;
                temp3 = temp2 * temp;
                temp4 = temp1 / dscale[i__ - 1];
                fc += temp4;
                erretm += f2c_abs(temp4);
                df += temp2;
                ddf += temp3;
            }
            else
            {
                goto L60;
            }
            /* L40: */
        }
        f = *finit + *tau * fc;
        erretm = (f2c_abs(*finit) + f2c_abs(*tau) * erretm) * 8. + f2c_abs(*tau) * df;
        if (f2c_abs(f) <= eps * erretm)
        {
            goto L60;
        }
        if (f <= 0.)
        {
            lbd = *tau;
        }
        else
        {
            ubd = *tau;
        }
        /* L50: */
    }
    *info = 1;
L60: /* Undo scaling */
    if (scale)
    {
        *tau *= sclinv;
    }
    return 0;
    /* End of DLAED6 */
}
/* dlaed6_ */
