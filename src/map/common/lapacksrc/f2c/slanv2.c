/* ../netlib/slanv2.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static real c_b4 = 1.f;
/* > \brief \b SLANV2 computes the Schur factorization of a real 2-by-2 nonsymmetric matrix in standard form. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download SLANV2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slanv2. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slanv2. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slanv2. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE SLANV2( A, B, C, D, RT1R, RT1I, RT2R, RT2I, CS, SN ) */
/* .. Scalar Arguments .. */
/* REAL A, B, C, CS, D, RT1I, RT1R, RT2I, RT2R, SN */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLANV2 computes the Schur factorization of a real 2-by-2 nonsymmetric */
/* > matrix in standard form: */
/* > */
/* > [ A B ] = [ CS -SN ] [ AA BB ] [ CS SN ] */
/* > [ C D ] [ SN CS ] [ CC DD ] [-SN CS ] */
/* > */
/* > where either */
/* > 1) CC = 0 so that AA and DD are real eigenvalues of the matrix, or */
/* > 2) AA = DD and BB*CC < 0, so that AA + or - sqrt(BB*CC) are complex */
/* > conjugate eigenvalues. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in,out] A */
/* > \verbatim */
/* > A is REAL */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* > B is REAL */
/* > \endverbatim */
/* > */
/* > \param[in,out] C */
/* > \verbatim */
/* > C is REAL */
/* > \endverbatim */
/* > */
/* > \param[in,out] D */
/* > \verbatim */
/* > D is REAL */
/* > On entry, the elements of the input matrix. */
/* > On exit, they are overwritten by the elements of the */
/* > standardised Schur form. */
/* > \endverbatim */
/* > */
/* > \param[out] RT1R */
/* > \verbatim */
/* > RT1R is REAL */
/* > \endverbatim */
/* > */
/* > \param[out] RT1I */
/* > \verbatim */
/* > RT1I is REAL */
/* > \endverbatim */
/* > */
/* > \param[out] RT2R */
/* > \verbatim */
/* > RT2R is REAL */
/* > \endverbatim */
/* > */
/* > \param[out] RT2I */
/* > \verbatim */
/* > RT2I is REAL */
/* > The real and imaginary parts of the eigenvalues. If the */
/* > eigenvalues are a complex conjugate pair, RT1I > 0. */
/* > \endverbatim */
/* > */
/* > \param[out] CS */
/* > \verbatim */
/* > CS is REAL */
/* > \endverbatim */
/* > */
/* > \param[out] SN */
/* > \verbatim */
/* > SN is REAL */
/* > Parameters of the rotation matrix. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date September 2012 */
/* > \ingroup realOTHERauxiliary */
/* > \par Further Details: */
/* ===================== */
/* > */
/* > \verbatim */
/* > */
/* > Modified by V. Sima, Research Institute for Informatics, Bucharest, */
/* > Romania, to reduce the risk of cancellation errors, */
/* > when computing real eigenvalues, and to ensure, if possible, that */
/* > f2c_abs(RT1R) >= f2c_abs(RT2R). */
/* > \endverbatim */
/* > */
/* ===================================================================== */
/* Subroutine */
int slanv2_(real *a, real *b, real *c__, real *d__, real * rt1r, real *rt1i, real *rt2r, real *rt2i, real *cs, real *sn)
{
    /* System generated locals */
    real r__1, r__2;
    /* Builtin functions */
    double r_sign(real *, real *), sqrt(doublereal);
    /* Local variables */
    real p, z__, aa, bb, cc, dd, cs1, sn1, sab, sac, eps, tau, temp, scale, bcmax, bcmis, sigma;
    extern real slapy2_(real *, real *), slamch_(char *);
    /* -- LAPACK auxiliary routine (version 3.4.2) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* September 2012 */
    /* .. Scalar Arguments .. */
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
    /* .. Executable Statements .. */
    eps = slamch_("P");
    if (*c__ == 0.f)
    {
        *cs = 1.f;
        *sn = 0.f;
        goto L10;
    }
    else if (*b == 0.f)
    {
        /* Swap rows and columns */
        *cs = 0.f;
        *sn = 1.f;
        temp = *d__;
        *d__ = *a;
        *a = temp;
        *b = -(*c__);
        *c__ = 0.f;
        goto L10;
    }
    else if (*a - *d__ == 0.f && r_sign(&c_b4, b) != r_sign(&c_b4, c__))
    {
        *cs = 1.f;
        *sn = 0.f;
        goto L10;
    }
    else
    {
        temp = *a - *d__;
        p = temp * .5f;
        /* Computing MAX */
        r__1 = f2c_abs(*b);
        r__2 = f2c_abs(*c__); // , expr subst
        bcmax = max(r__1,r__2);
        /* Computing MIN */
        r__1 = f2c_abs(*b);
        r__2 = f2c_abs(*c__); // , expr subst
        bcmis = min(r__1,r__2) * r_sign(&c_b4, b) * r_sign(&c_b4, c__);
        /* Computing MAX */
        r__1 = f2c_abs(p);
        scale = max(r__1,bcmax);
        z__ = p / scale * p + bcmax / scale * bcmis;
        /* If Z is of the order of the machine accuracy, postpone the */
        /* decision on the nature of eigenvalues */
        if (z__ >= eps * 4.f)
        {
            /* Real eigenvalues. Compute A and D. */
            r__1 = sqrt(scale) * sqrt(z__);
            z__ = p + r_sign(&r__1, &p);
            *a = *d__ + z__;
            *d__ -= bcmax / z__ * bcmis;
            /* Compute B and the rotation matrix */
            tau = slapy2_(c__, &z__);
            *cs = z__ / tau;
            *sn = *c__ / tau;
            *b -= *c__;
            *c__ = 0.f;
        }
        else
        {
            /* Complex eigenvalues, or real (almost) equal eigenvalues. */
            /* Make diagonal elements equal. */
            sigma = *b + *c__;
            tau = slapy2_(&sigma, &temp);
            *cs = sqrt((f2c_abs(sigma) / tau + 1.f) * .5f);
            *sn = -(p / (tau * *cs)) * r_sign(&c_b4, &sigma);
            /* Compute [ AA BB ] = [ A B ] [ CS -SN ] */
            /* [ CC DD ] [ C D ] [ SN CS ] */
            aa = *a * *cs + *b * *sn;
            bb = -(*a) * *sn + *b * *cs;
            cc = *c__ * *cs + *d__ * *sn;
            dd = -(*c__) * *sn + *d__ * *cs;
            /* Compute [ A B ] = [ CS SN ] [ AA BB ] */
            /* [ C D ] [-SN CS ] [ CC DD ] */
            *a = aa * *cs + cc * *sn;
            *b = bb * *cs + dd * *sn;
            *c__ = -aa * *sn + cc * *cs;
            *d__ = -bb * *sn + dd * *cs;
            temp = (*a + *d__) * .5f;
            *a = temp;
            *d__ = temp;
            if (*c__ != 0.f)
            {
                if (*b != 0.f)
                {
                    if (r_sign(&c_b4, b) == r_sign(&c_b4, c__))
                    {
                        /* Real eigenvalues: reduce to upper triangular form */
                        sab = sqrt((f2c_abs(*b)));
                        sac = sqrt((f2c_abs(*c__)));
                        r__1 = sab * sac;
                        p = r_sign(&r__1, c__);
                        tau = 1.f / sqrt((r__1 = *b + *c__, f2c_abs(r__1)));
                        *a = temp + p;
                        *d__ = temp - p;
                        *b -= *c__;
                        *c__ = 0.f;
                        cs1 = sab * tau;
                        sn1 = sac * tau;
                        temp = *cs * cs1 - *sn * sn1;
                        *sn = *cs * sn1 + *sn * cs1;
                        *cs = temp;
                    }
                }
                else
                {
                    *b = -(*c__);
                    *c__ = 0.f;
                    temp = *cs;
                    *cs = -(*sn);
                    *sn = temp;
                }
            }
        }
    }
L10: /* Store eigenvalues in (RT1R,RT1I) and (RT2R,RT2I). */
    *rt1r = *a;
    *rt2r = *d__;
    if (*c__ == 0.f)
    {
        *rt1i = 0.f;
        *rt2i = 0.f;
    }
    else
    {
        *rt1i = sqrt((f2c_abs(*b))) * sqrt((f2c_abs(*c__)));
        *rt2i = -(*rt1i);
    }
    return 0;
    /* End of SLANV2 */
}
/* slanv2_ */
