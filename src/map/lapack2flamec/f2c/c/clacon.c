/* ../netlib/clacon.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
/* > \brief \b CLACON estimates the 1-norm of a square matrix, using reverse communication for evaluating matr ix-vector products. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CLACON + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clacon. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clacon. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clacon. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CLACON( N, V, X, EST, KASE ) */
/* .. Scalar Arguments .. */
/* INTEGER KASE, N */
/* REAL EST */
/* .. */
/* .. Array Arguments .. */
/* COMPLEX V( N ), X( N ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CLACON estimates the 1-norm of a square, complex matrix A. */
/* > Reverse communication is used for evaluating matrix-vector products. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the matrix. N >= 1. */
/* > \endverbatim */
/* > */
/* > \param[out] V */
/* > \verbatim */
/* > V is COMPLEX array, dimension (N) */
/* > On the final return, V = A*W, where EST = norm(V)/norm(W) */
/* > (W is not returned). */
/* > \endverbatim */
/* > */
/* > \param[in,out] X */
/* > \verbatim */
/* > X is COMPLEX array, dimension (N) */
/* > On an intermediate return, X should be overwritten by */
/* > A * X, if KASE=1, */
/* > A**H * X, if KASE=2, */
/* > where A**H is the conjugate transpose of A, and CLACON must be */
/* > re-called with all the other parameters unchanged. */
/* > \endverbatim */
/* > */
/* > \param[in,out] EST */
/* > \verbatim */
/* > EST is REAL */
/* > On entry with KASE = 1 or 2 and JUMP = 3, EST should be */
/* > unchanged from the previous call to CLACON. */
/* > On exit, EST is an estimate (a lower bound) for norm(A). */
/* > \endverbatim */
/* > */
/* > \param[in,out] KASE */
/* > \verbatim */
/* > KASE is INTEGER */
/* > On the initial call to CLACON, KASE should be 0. */
/* > On an intermediate return, KASE will be 1 or 2, indicating */
/* > whether X should be overwritten by A * X or A**H * X. */
/* > On the final return from CLACON, KASE will again be 0. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date September 2012 */
/* > \ingroup complexOTHERauxiliary */
/* > \par Further Details: */
/* ===================== */
/* > */
/* > Originally named CONEST, dated March 16, 1988. \n */
/* > Last modified: April, 1999 */
/* > \par Contributors: */
/* ================== */
/* > */
/* > Nick Higham, University of Manchester */
/* > \par References: */
/* ================ */
/* > */
/* > N.J. Higham, "FORTRAN codes for estimating the one-norm of */
/* > a real or complex matrix, with applications to condition estimation", */
/* > ACM Trans. Math. Soft., vol. 14, no. 4, pp. 381-396, December 1988. */
/* > */
/* ===================================================================== */
/* Subroutine */
int clacon_(integer *n, complex *v, complex *x, real *est, integer *kase)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    real r__1, r__2;
    complex q__1;
    /* Builtin functions */
    double c_abs(complex *), r_imag(complex *);
    /* Local variables */
    static integer i__, j, iter;
    static real temp;
    static integer jump;
    static real absxi;
    static integer jlast;
    extern /* Subroutine */
    int ccopy_(integer *, complex *, integer *, complex *, integer *);
    extern integer icmax1_(integer *, complex *, integer *);
    extern real scsum1_(integer *, complex *, integer *), slamch_(char *);
    static real safmin, altsgn, estold;
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
    /* .. Save statement .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Parameter adjustments */
    --x;
    --v;
    /* Function Body */
    safmin = slamch_("Safe minimum");
    if (*kase == 0)
    {
        i__1 = *n;
        for (i__ = 1;
                i__ <= i__1;
                ++i__)
        {
            i__2 = i__;
            r__1 = 1.f / (real) (*n);
            q__1.r = r__1;
            q__1.i = 0.f; // , expr subst
            x[i__2].r = q__1.r;
            x[i__2].i = q__1.i; // , expr subst
            /* L10: */
        }
        *kase = 1;
        jump = 1;
        return 0;
    }
    switch (jump)
    {
    case 1:
        goto L20;
    case 2:
        goto L40;
    case 3:
        goto L70;
    case 4:
        goto L90;
    case 5:
        goto L120;
    }
    /* ................ ENTRY (JUMP = 1) */
    /* FIRST ITERATION. X HAS BEEN OVERWRITTEN BY A*X. */
L20:
    if (*n == 1)
    {
        v[1].r = x[1].r;
        v[1].i = x[1].i; // , expr subst
        *est = c_abs(&v[1]);
        /* ... QUIT */
        goto L130;
    }
    *est = scsum1_(n, &x[1], &c__1);
    i__1 = *n;
    for (i__ = 1;
            i__ <= i__1;
            ++i__)
    {
        absxi = c_abs(&x[i__]);
        if (absxi > safmin)
        {
            i__2 = i__;
            i__3 = i__;
            r__1 = x[i__3].r / absxi;
            r__2 = r_imag(&x[i__]) / absxi;
            q__1.r = r__1;
            q__1.i = r__2; // , expr subst
            x[i__2].r = q__1.r;
            x[i__2].i = q__1.i; // , expr subst
        }
        else
        {
            i__2 = i__;
            x[i__2].r = 1.f;
            x[i__2].i = 0.f; // , expr subst
        }
        /* L30: */
    }
    *kase = 2;
    jump = 2;
    return 0;
    /* ................ ENTRY (JUMP = 2) */
    /* FIRST ITERATION. X HAS BEEN OVERWRITTEN BY CTRANS(A)*X. */
L40:
    j = icmax1_(n, &x[1], &c__1);
    iter = 2;
    /* MAIN LOOP - ITERATIONS 2,3,...,ITMAX. */
L50:
    i__1 = *n;
    for (i__ = 1;
            i__ <= i__1;
            ++i__)
    {
        i__2 = i__;
        x[i__2].r = 0.f;
        x[i__2].i = 0.f; // , expr subst
        /* L60: */
    }
    i__1 = j;
    x[i__1].r = 1.f;
    x[i__1].i = 0.f; // , expr subst
    *kase = 1;
    jump = 3;
    return 0;
    /* ................ ENTRY (JUMP = 3) */
    /* X HAS BEEN OVERWRITTEN BY A*X. */
L70:
    ccopy_(n, &x[1], &c__1, &v[1], &c__1);
    estold = *est;
    *est = scsum1_(n, &v[1], &c__1);
    /* TEST FOR CYCLING. */
    if (*est <= estold)
    {
        goto L100;
    }
    i__1 = *n;
    for (i__ = 1;
            i__ <= i__1;
            ++i__)
    {
        absxi = c_abs(&x[i__]);
        if (absxi > safmin)
        {
            i__2 = i__;
            i__3 = i__;
            r__1 = x[i__3].r / absxi;
            r__2 = r_imag(&x[i__]) / absxi;
            q__1.r = r__1;
            q__1.i = r__2; // , expr subst
            x[i__2].r = q__1.r;
            x[i__2].i = q__1.i; // , expr subst
        }
        else
        {
            i__2 = i__;
            x[i__2].r = 1.f;
            x[i__2].i = 0.f; // , expr subst
        }
        /* L80: */
    }
    *kase = 2;
    jump = 4;
    return 0;
    /* ................ ENTRY (JUMP = 4) */
    /* X HAS BEEN OVERWRITTEN BY CTRANS(A)*X. */
L90:
    jlast = j;
    j = icmax1_(n, &x[1], &c__1);
    if (c_abs(&x[jlast]) != c_abs(&x[j]) && iter < 5)
    {
        ++iter;
        goto L50;
    }
    /* ITERATION COMPLETE. FINAL STAGE. */
L100:
    altsgn = 1.f;
    i__1 = *n;
    for (i__ = 1;
            i__ <= i__1;
            ++i__)
    {
        i__2 = i__;
        r__1 = altsgn * ((real) (i__ - 1) / (real) (*n - 1) + 1.f);
        q__1.r = r__1;
        q__1.i = 0.f; // , expr subst
        x[i__2].r = q__1.r;
        x[i__2].i = q__1.i; // , expr subst
        altsgn = -altsgn;
        /* L110: */
    }
    *kase = 1;
    jump = 5;
    return 0;
    /* ................ ENTRY (JUMP = 5) */
    /* X HAS BEEN OVERWRITTEN BY A*X. */
L120:
    temp = scsum1_(n, &x[1], &c__1) / (real) (*n * 3) * 2.f;
    if (temp > *est)
    {
        ccopy_(n, &x[1], &c__1, &v[1], &c__1);
        *est = temp;
    }
L130:
    *kase = 0;
    return 0;
    /* End of CLACON */
}
/* clacon_ */
