/* ../netlib/csrscl.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief \b CSRSCL multiplies a vector by the reciprocal of a real scalar. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CSRSCL + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/csrscl. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/csrscl. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/csrscl. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CSRSCL( N, SA, SX, INCX ) */
/* .. Scalar Arguments .. */
/* INTEGER INCX, N */
/* REAL SA */
/* .. */
/* .. Array Arguments .. */
/* COMPLEX SX( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CSRSCL multiplies an n-element complex vector x by the real scalar */
/* > 1/a. This is done without overflow or underflow as long as */
/* > the final result x/a does not overflow or underflow. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The number of components of the vector x. */
/* > \endverbatim */
/* > */
/* > \param[in] SA */
/* > \verbatim */
/* > SA is REAL */
/* > The scalar a which is used to divide each component of x. */
/* > SA must be >= 0, or the subroutine will divide by zero. */
/* > \endverbatim */
/* > */
/* > \param[in,out] SX */
/* > \verbatim */
/* > SX is COMPLEX array, dimension */
/* > (1+(N-1)*f2c_abs(INCX)) */
/* > The n-element vector x. */
/* > \endverbatim */
/* > */
/* > \param[in] INCX */
/* > \verbatim */
/* > INCX is INTEGER */
/* > The increment between successive values of the vector SX. */
/* > > 0: SX(1) = X(1) and SX(1+(i-1)*INCX) = x(i), 1< i<= n */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date September 2012 */
/* > \ingroup complexOTHERauxiliary */
/* ===================================================================== */
/* Subroutine */
int csrscl_(integer *n, real *sa, complex *sx, integer *incx)
{
    real mul, cden;
    logical done;
    real cnum, cden1, cnum1;
    extern /* Subroutine */
    int slabad_(real *, real *);
    extern real slamch_(char *);
    extern /* Subroutine */
    int csscal_(integer *, real *, complex *, integer *);
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
    /* .. External Subroutines .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Quick return if possible */
    /* Parameter adjustments */
    --sx;
    /* Function Body */
    if (*n <= 0)
    {
        return 0;
    }
    /* Get machine parameters */
    smlnum = slamch_("S");
    bignum = 1.f / smlnum;
    slabad_(&smlnum, &bignum);
    /* Initialize the denominator to SA and the numerator to 1. */
    cden = *sa;
    cnum = 1.f;
L10:
    cden1 = cden * smlnum;
    cnum1 = cnum / bignum;
    if (f2c_abs(cden1) > f2c_abs(cnum) && cnum != 0.f)
    {
        /* Pre-multiply X by SMLNUM if CDEN is large compared to CNUM. */
        mul = smlnum;
        done = FALSE_;
        cden = cden1;
    }
    else if (f2c_abs(cnum1) > f2c_abs(cden))
    {
        /* Pre-multiply X by BIGNUM if CDEN is small compared to CNUM. */
        mul = bignum;
        done = FALSE_;
        cnum = cnum1;
    }
    else
    {
        /* Multiply X by CNUM / CDEN and return. */
        mul = cnum / cden;
        done = TRUE_;
    }
    /* Scale the vector X by MUL */
    csscal_(n, &mul, &sx[1], incx);
    if (! done)
    {
        goto L10;
    }
    return 0;
    /* End of CSRSCL */
}
/* csrscl_ */
