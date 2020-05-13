/* ../netlib/slarrk.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief \b SLARRK computes one eigenvalue of a symmetric tridiagonal matrix T to suitable accuracy. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download SLARRK + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slarrk. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slarrk. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slarrk. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE SLARRK( N, IW, GL, GU, */
/* D, E2, PIVMIN, RELTOL, W, WERR, INFO) */
/* .. Scalar Arguments .. */
/* INTEGER INFO, IW, N */
/* REAL PIVMIN, RELTOL, GL, GU, W, WERR */
/* .. */
/* .. Array Arguments .. */
/* REAL D( * ), E2( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLARRK computes one eigenvalue of a symmetric tridiagonal */
/* > matrix T to suitable accuracy. This is an auxiliary code to be */
/* > called from SSTEMR. */
/* > */
/* > To avoid overflow, the matrix must be scaled so that its */
/* > largest element is no greater than overflow**(1/2) * underflow**(1/4) in absolute value, and for greatest */
/* > accuracy, it should not be much smaller than that. */
/* > */
/* > See W. Kahan "Accurate Eigenvalues of a Symmetric Tridiagonal */
/* > Matrix", Report CS41, Computer Science Dept., Stanford */
/* > University, July 21, 1966. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the tridiagonal matrix T. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] IW */
/* > \verbatim */
/* > IW is INTEGER */
/* > The index of the eigenvalues to be returned. */
/* > \endverbatim */
/* > */
/* > \param[in] GL */
/* > \verbatim */
/* > GL is REAL */
/* > \endverbatim */
/* > */
/* > \param[in] GU */
/* > \verbatim */
/* > GU is REAL */
/* > An upper and a lower bound on the eigenvalue. */
/* > \endverbatim */
/* > */
/* > \param[in] D */
/* > \verbatim */
/* > D is REAL array, dimension (N) */
/* > The n diagonal elements of the tridiagonal matrix T. */
/* > \endverbatim */
/* > */
/* > \param[in] E2 */
/* > \verbatim */
/* > E2 is REAL array, dimension (N-1) */
/* > The (n-1) squared off-diagonal elements of the tridiagonal matrix T. */
/* > \endverbatim */
/* > */
/* > \param[in] PIVMIN */
/* > \verbatim */
/* > PIVMIN is REAL */
/* > The minimum pivot allowed in the Sturm sequence for T. */
/* > \endverbatim */
/* > */
/* > \param[in] RELTOL */
/* > \verbatim */
/* > RELTOL is REAL */
/* > The minimum relative width of an interval. When an interval */
/* > is narrower than RELTOL times the larger (in */
/* > magnitude) endpoint, then it is considered to be */
/* > sufficiently small, i.e., converged. Note: this should */
/* > always be at least radix*machine epsilon. */
/* > \endverbatim */
/* > */
/* > \param[out] W */
/* > \verbatim */
/* > W is REAL */
/* > \endverbatim */
/* > */
/* > \param[out] WERR */
/* > \verbatim */
/* > WERR is REAL */
/* > The error bound on the corresponding eigenvalue approximation */
/* > in W. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: Eigenvalue converged */
/* > = -1: Eigenvalue did NOT converge */
/* > \endverbatim */
/* > \par Internal Parameters: */
/* ========================= */
/* > */
/* > \verbatim */
/* > FUDGE REAL , default = 2 */
/* > A "fudge factor" to widen the Gershgorin intervals. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date September 2012 */
/* > \ingroup auxOTHERauxiliary */
/* ===================================================================== */
/* Subroutine */
int slarrk_(integer *n, integer *iw, real *gl, real *gu, real *d__, real *e2, real *pivmin, real *reltol, real *w, real *werr, integer *info)
{
    /* System generated locals */
    integer i__1;
    real r__1, r__2;
    /* Builtin functions */
    double log(doublereal);
    /* Local variables */
    integer i__, it;
    real mid, eps, tmp1, tmp2, left, atoli, right;
    integer itmax;
    real rtoli, tnorm;
    extern real slamch_(char *);
    integer negcnt;
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
    /* .. Executable Statements .. */
    /* Get machine constants */
    /* Parameter adjustments */
    --e2;
    --d__;
    /* Function Body */
    eps = slamch_("P");
    /* Computing MAX */
    r__1 = f2c_abs(*gl);
    r__2 = f2c_abs(*gu); // , expr subst
    tnorm = max(r__1,r__2);
    rtoli = *reltol;
    atoli = *pivmin * 4.f;
    itmax = (integer) ((log(tnorm + *pivmin) - log(*pivmin)) / log(2.f)) + 2;
    *info = -1;
    left = *gl - tnorm * 2.f * eps * *n - *pivmin * 4.f;
    right = *gu + tnorm * 2.f * eps * *n + *pivmin * 4.f;
    it = 0;
L10: /* Check if interval converged or maximum number of iterations reached */
    tmp1 = (r__1 = right - left, f2c_abs(r__1));
    /* Computing MAX */
    r__1 = f2c_abs(right);
    r__2 = f2c_abs(left); // , expr subst
    tmp2 = max(r__1,r__2);
    /* Computing MAX */
    r__1 = max(atoli,*pivmin);
    r__2 = rtoli * tmp2; // , expr subst
    if (tmp1 < max(r__1,r__2))
    {
        *info = 0;
        goto L30;
    }
    if (it > itmax)
    {
        goto L30;
    }
    /* Count number of negative pivots for mid-point */
    ++it;
    mid = (left + right) * .5f;
    negcnt = 0;
    tmp1 = d__[1] - mid;
    if (f2c_abs(tmp1) < *pivmin)
    {
        tmp1 = -(*pivmin);
    }
    if (tmp1 <= 0.f)
    {
        ++negcnt;
    }
    i__1 = *n;
    for (i__ = 2;
            i__ <= i__1;
            ++i__)
    {
        tmp1 = d__[i__] - e2[i__ - 1] / tmp1 - mid;
        if (f2c_abs(tmp1) < *pivmin)
        {
            tmp1 = -(*pivmin);
        }
        if (tmp1 <= 0.f)
        {
            ++negcnt;
        }
        /* L20: */
    }
    if (negcnt >= *iw)
    {
        right = mid;
    }
    else
    {
        left = mid;
    }
    goto L10;
L30: /* Converged or maximum number of iterations reached */
    *w = (left + right) * .5f;
    *werr = (r__1 = right - left, f2c_abs(r__1)) * .5f;
    return 0;
    /* End of SLARRK */
}
/* slarrk_ */
