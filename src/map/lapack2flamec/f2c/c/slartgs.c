/* ../netlib/slartgs.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief \b SLARTGS generates a plane rotation designed to introduce a bulge in implicit QR iteration for t he bidiagonal SVD problem. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download SLARTGS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slartgs .f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slartgs .f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slartgs .f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE SLARTGS( X, Y, SIGMA, CS, SN ) */
/* .. Scalar Arguments .. */
/* REAL CS, SIGMA, SN, X, Y */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLARTGS generates a plane rotation designed to introduce a bulge in */
/* > Golub-Reinsch-style implicit QR iteration for the bidiagonal SVD */
/* > problem. X and Y are the top-row entries, and SIGMA is the shift. */
/* > The computed CS and SN define a plane rotation satisfying */
/* > */
/* > [ CS SN ] . [ X^2 - SIGMA ] = [ R ], */
/* > [ -SN CS ] [ X * Y ] [ 0 ] */
/* > */
/* > with R nonnegative. If X^2 - SIGMA and X * Y are 0, then the */
/* > rotation is by PI/2. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] X */
/* > \verbatim */
/* > X is REAL */
/* > The (1,1) entry of an upper bidiagonal matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] Y */
/* > \verbatim */
/* > Y is REAL */
/* > The (1,2) entry of an upper bidiagonal matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] SIGMA */
/* > \verbatim */
/* > SIGMA is REAL */
/* > The shift. */
/* > \endverbatim */
/* > */
/* > \param[out] CS */
/* > \verbatim */
/* > CS is REAL */
/* > The cosine of the rotation. */
/* > \endverbatim */
/* > */
/* > \param[out] SN */
/* > \verbatim */
/* > SN is REAL */
/* > The sine of the rotation. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date September 2012 */
/* > \ingroup auxOTHERcomputational */
/* ===================================================================== */
/* Subroutine */
int slartgs_(real *x, real *y, real *sigma, real *cs, real * sn)
{
    real r__, s, w, z__;
    extern real slamch_(char *);
    real thresh;
    extern /* Subroutine */
    int slartgp_(real *, real *, real *, real *, real *);
    /* -- LAPACK computational routine (version 3.4.2) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* September 2012 */
    /* .. Scalar Arguments .. */
    /* .. */
    /* =================================================================== */
    /* .. Parameters .. */
    /* .. */
    /* .. Local Scalars .. */
    /* .. */
    /* .. External Functions .. */
    /* .. Executable Statements .. */
    thresh = slamch_("E");
    /* Compute the first column of B**T*B - SIGMA^2*I, up to a scale */
    /* factor. */
    if (*sigma == 0.f && f2c_abs(*x) < thresh || f2c_abs(*x) == *sigma && *y == 0.f)
    {
        z__ = 0.f;
        w = 0.f;
    }
    else if (*sigma == 0.f)
    {
        if (*x >= 0.f)
        {
            z__ = *x;
            w = *y;
        }
        else
        {
            z__ = -(*x);
            w = -(*y);
        }
    }
    else if (f2c_abs(*x) < thresh)
    {
        z__ = -(*sigma) * *sigma;
        w = 0.f;
    }
    else
    {
        if (*x >= 0.f)
        {
            s = 1.f;
        }
        else
        {
            s = -1.f;
        }
        z__ = s * (f2c_abs(*x) - *sigma) * (s + *sigma / *x);
        w = s * *y;
    }
    /* Generate the rotation. */
    /* CALL SLARTGP( Z, W, CS, SN, R ) might seem more natural;
    */
    /* reordering the arguments ensures that if Z = 0 then the rotation */
    /* is by PI/2. */
    slartgp_(&w, &z__, sn, cs, &r__);
    return 0;
    /* End SLARTGS */
}
/* slartgs_ */
