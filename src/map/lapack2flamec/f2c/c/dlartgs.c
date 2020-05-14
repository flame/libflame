/* ../netlib/dlartgs.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief \b DLARTGS generates a plane rotation designed to introduce a bulge in implicit QR iteration for t he bidiagonal SVD problem. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download DLARTGS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlartgs .f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlartgs .f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlartgs .f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE DLARTGS( X, Y, SIGMA, CS, SN ) */
/* .. Scalar Arguments .. */
/* DOUBLE PRECISION CS, SIGMA, SN, X, Y */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLARTGS generates a plane rotation designed to introduce a bulge in */
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
/* > X is DOUBLE PRECISION */
/* > The (1,1) entry of an upper bidiagonal matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] Y */
/* > \verbatim */
/* > Y is DOUBLE PRECISION */
/* > The (1,2) entry of an upper bidiagonal matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] SIGMA */
/* > \verbatim */
/* > SIGMA is DOUBLE PRECISION */
/* > The shift. */
/* > \endverbatim */
/* > */
/* > \param[out] CS */
/* > \verbatim */
/* > CS is DOUBLE PRECISION */
/* > The cosine of the rotation. */
/* > \endverbatim */
/* > */
/* > \param[out] SN */
/* > \verbatim */
/* > SN is DOUBLE PRECISION */
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
int dlartgs_(doublereal *x, doublereal *y, doublereal *sigma, doublereal *cs, doublereal *sn)
{
    doublereal r__, s, w, z__;
    extern doublereal dlamch_(char *);
    doublereal thresh;
    extern /* Subroutine */
    int dlartgp_(doublereal *, doublereal *, doublereal *, doublereal *, doublereal *);
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
    thresh = dlamch_("E");
    /* Compute the first column of B**T*B - SIGMA^2*I, up to a scale */
    /* factor. */
    if (*sigma == 0. && f2c_abs(*x) < thresh || f2c_abs(*x) == *sigma && *y == 0.)
    {
        z__ = 0.;
        w = 0.;
    }
    else if (*sigma == 0.)
    {
        if (*x >= 0.)
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
        w = 0.;
    }
    else
    {
        if (*x >= 0.)
        {
            s = 1.;
        }
        else
        {
            s = -1.;
        }
        z__ = s * (f2c_abs(*x) - *sigma) * (s + *sigma / *x);
        w = s * *y;
    }
    /* Generate the rotation. */
    /* CALL DLARTGP( Z, W, CS, SN, R ) might seem more natural;
    */
    /* reordering the arguments ensures that if Z = 0 then the rotation */
    /* is by PI/2. */
    dlartgp_(&w, &z__, sn, cs, &r__);
    return 0;
    /* End DLARTGS */
}
/* dlartgs_ */
