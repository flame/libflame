/* ../netlib/slaed5.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief \b SLAED5 used by sstedc. Solves the 2-by-2 secular equation. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download SLAED5 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slaed5. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slaed5. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slaed5. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE SLAED5( I, D, Z, DELTA, RHO, DLAM ) */
/* .. Scalar Arguments .. */
/* INTEGER I */
/* REAL DLAM, RHO */
/* .. */
/* .. Array Arguments .. */
/* REAL D( 2 ), DELTA( 2 ), Z( 2 ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > This subroutine computes the I-th eigenvalue of a symmetric rank-one */
/* > modification of a 2-by-2 diagonal matrix */
/* > */
/* > diag( D ) + RHO * Z * transpose(Z) . */
/* > */
/* > The diagonal elements in the array D are assumed to satisfy */
/* > */
/* > D(i) < D(j) for i < j . */
/* > */
/* > We also assume RHO > 0 and that the Euclidean norm of the vector */
/* > Z is one. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] I */
/* > \verbatim */
/* > I is INTEGER */
/* > The index of the eigenvalue to be computed. I = 1 or I = 2. */
/* > \endverbatim */
/* > */
/* > \param[in] D */
/* > \verbatim */
/* > D is REAL array, dimension (2) */
/* > The original eigenvalues. We assume D(1) < D(2). */
/* > \endverbatim */
/* > */
/* > \param[in] Z */
/* > \verbatim */
/* > Z is REAL array, dimension (2) */
/* > The components of the updating vector. */
/* > \endverbatim */
/* > */
/* > \param[out] DELTA */
/* > \verbatim */
/* > DELTA is REAL array, dimension (2) */
/* > The vector DELTA contains the information necessary */
/* > to construct the eigenvectors. */
/* > \endverbatim */
/* > */
/* > \param[in] RHO */
/* > \verbatim */
/* > RHO is REAL */
/* > The scalar in the symmetric updating formula. */
/* > \endverbatim */
/* > */
/* > \param[out] DLAM */
/* > \verbatim */
/* > DLAM is REAL */
/* > The computed lambda_I, the I-th updated eigenvalue. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date September 2012 */
/* > \ingroup auxOTHERcomputational */
/* > \par Contributors: */
/* ================== */
/* > */
/* > Ren-Cang Li, Computer Science Division, University of California */
/* > at Berkeley, USA */
/* > */
/* ===================================================================== */
/* Subroutine */
int slaed5_(integer *i__, real *d__, real *z__, real *delta, real *rho, real *dlam)
{
    /* System generated locals */
    real r__1;
    /* Builtin functions */
    double sqrt(doublereal);
    /* Local variables */
    real b, c__, w, del, tau, temp;
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
    /* .. Local Scalars .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Parameter adjustments */
    --delta;
    --z__;
    --d__;
    /* Function Body */
    del = d__[2] - d__[1];
    if (*i__ == 1)
    {
        w = *rho * 2.f * (z__[2] * z__[2] - z__[1] * z__[1]) / del + 1.f;
        if (w > 0.f)
        {
            b = del + *rho * (z__[1] * z__[1] + z__[2] * z__[2]);
            c__ = *rho * z__[1] * z__[1] * del;
            /* B > ZERO, always */
            tau = c__ * 2.f / (b + sqrt((r__1 = b * b - c__ * 4.f, f2c_abs(r__1))) );
            *dlam = d__[1] + tau;
            delta[1] = -z__[1] / tau;
            delta[2] = z__[2] / (del - tau);
        }
        else
        {
            b = -del + *rho * (z__[1] * z__[1] + z__[2] * z__[2]);
            c__ = *rho * z__[2] * z__[2] * del;
            if (b > 0.f)
            {
                tau = c__ * -2.f / (b + sqrt(b * b + c__ * 4.f));
            }
            else
            {
                tau = (b - sqrt(b * b + c__ * 4.f)) / 2.f;
            }
            *dlam = d__[2] + tau;
            delta[1] = -z__[1] / (del + tau);
            delta[2] = -z__[2] / tau;
        }
        temp = sqrt(delta[1] * delta[1] + delta[2] * delta[2]);
        delta[1] /= temp;
        delta[2] /= temp;
    }
    else
    {
        /* Now I=2 */
        b = -del + *rho * (z__[1] * z__[1] + z__[2] * z__[2]);
        c__ = *rho * z__[2] * z__[2] * del;
        if (b > 0.f)
        {
            tau = (b + sqrt(b * b + c__ * 4.f)) / 2.f;
        }
        else
        {
            tau = c__ * 2.f / (-b + sqrt(b * b + c__ * 4.f));
        }
        *dlam = d__[2] + tau;
        delta[1] = -z__[1] / (del + tau);
        delta[2] = -z__[2] / tau;
        temp = sqrt(delta[1] * delta[1] + delta[2] * delta[2]);
        delta[1] /= temp;
        delta[2] /= temp;
    }
    return 0;
    /* End OF SLAED5 */
}
/* slaed5_ */
