/* ../netlib/dlacon.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
static doublereal c_b11 = 1.;
/* > \brief \b DLACON estimates the 1-norm of a square matrix, using reverse communication for evaluating matr ix-vector products. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download DLACON + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlacon. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlacon. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlacon. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE DLACON( N, V, X, ISGN, EST, KASE ) */
/* .. Scalar Arguments .. */
/* INTEGER KASE, N */
/* DOUBLE PRECISION EST */
/* .. */
/* .. Array Arguments .. */
/* INTEGER ISGN( * ) */
/* DOUBLE PRECISION V( * ), X( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLACON estimates the 1-norm of a square, real matrix A. */
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
/* > V is DOUBLE PRECISION array, dimension (N) */
/* > On the final return, V = A*W, where EST = norm(V)/norm(W) */
/* > (W is not returned). */
/* > \endverbatim */
/* > */
/* > \param[in,out] X */
/* > \verbatim */
/* > X is DOUBLE PRECISION array, dimension (N) */
/* > On an intermediate return, X should be overwritten by */
/* > A * X, if KASE=1, */
/* > A**T * X, if KASE=2, */
/* > and DLACON must be re-called with all the other parameters */
/* > unchanged. */
/* > \endverbatim */
/* > */
/* > \param[out] ISGN */
/* > \verbatim */
/* > ISGN is INTEGER array, dimension (N) */
/* > \endverbatim */
/* > */
/* > \param[in,out] EST */
/* > \verbatim */
/* > EST is DOUBLE PRECISION */
/* > On entry with KASE = 1 or 2 and JUMP = 3, EST should be */
/* > unchanged from the previous call to DLACON. */
/* > On exit, EST is an estimate (a lower bound) for norm(A). */
/* > \endverbatim */
/* > */
/* > \param[in,out] KASE */
/* > \verbatim */
/* > KASE is INTEGER */
/* > On the initial call to DLACON, KASE should be 0. */
/* > On an intermediate return, KASE will be 1 or 2, indicating */
/* > whether X should be overwritten by A * X or A**T * X. */
/* > On the final return from DLACON, KASE will again be 0. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date September 2012 */
/* > \ingroup doubleOTHERauxiliary */
/* > \par Contributors: */
/* ================== */
/* > */
/* > Nick Higham, University of Manchester. \n */
/* > Originally named SONEST, dated March 16, 1988. */
/* > \par References: */
/* ================ */
/* > */
/* > N.J. Higham, "FORTRAN codes for estimating the one-norm of */
/* > a real or complex matrix, with applications to condition estimation", */
/* > ACM Trans. Math. Soft., vol. 14, no. 4, pp. 381-396, December 1988. */
/* > */
/* ===================================================================== */
/* Subroutine */
int dlacon_(integer *n, doublereal *v, doublereal *x, integer *isgn, doublereal *est, integer *kase)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;
    /* Builtin functions */
    double d_sign(doublereal *, doublereal *);
    integer i_dnnt(doublereal *);
    /* Local variables */
    static integer i__, j, iter;
    static doublereal temp;
    static integer jump;
    extern doublereal dasum_(integer *, doublereal *, integer *);
    static integer jlast;
    extern /* Subroutine */
    int dcopy_(integer *, doublereal *, integer *, doublereal *, integer *);
    extern integer idamax_(integer *, doublereal *, integer *);
    static doublereal altsgn, estold;
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
    --isgn;
    --x;
    --v;
    /* Function Body */
    if (*kase == 0)
    {
        i__1 = *n;
        for (i__ = 1;
                i__ <= i__1;
                ++i__)
        {
            x[i__] = 1. / (doublereal) (*n);
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
        goto L110;
    case 5:
        goto L140;
    }
    /* ................ ENTRY (JUMP = 1) */
    /* FIRST ITERATION. X HAS BEEN OVERWRITTEN BY A*X. */
L20:
    if (*n == 1)
    {
        v[1] = x[1];
        *est = f2c_abs(v[1]);
        /* ... QUIT */
        goto L150;
    }
    *est = dasum_(n, &x[1], &c__1);
    i__1 = *n;
    for (i__ = 1;
            i__ <= i__1;
            ++i__)
    {
        x[i__] = d_sign(&c_b11, &x[i__]);
        isgn[i__] = i_dnnt(&x[i__]);
        /* L30: */
    }
    *kase = 2;
    jump = 2;
    return 0;
    /* ................ ENTRY (JUMP = 2) */
    /* FIRST ITERATION. X HAS BEEN OVERWRITTEN BY TRANSPOSE(A)*X. */
L40:
    j = idamax_(n, &x[1], &c__1);
    iter = 2;
    /* MAIN LOOP - ITERATIONS 2,3,...,ITMAX. */
L50:
    i__1 = *n;
    for (i__ = 1;
            i__ <= i__1;
            ++i__)
    {
        x[i__] = 0.;
        /* L60: */
    }
    x[j] = 1.;
    *kase = 1;
    jump = 3;
    return 0;
    /* ................ ENTRY (JUMP = 3) */
    /* X HAS BEEN OVERWRITTEN BY A*X. */
L70:
    dcopy_(n, &x[1], &c__1, &v[1], &c__1);
    estold = *est;
    *est = dasum_(n, &v[1], &c__1);
    i__1 = *n;
    for (i__ = 1;
            i__ <= i__1;
            ++i__)
    {
        d__1 = d_sign(&c_b11, &x[i__]);
        if (i_dnnt(&d__1) != isgn[i__])
        {
            goto L90;
        }
        /* L80: */
    }
    /* REPEATED SIGN VECTOR DETECTED, HENCE ALGORITHM HAS CONVERGED. */
    goto L120;
L90: /* TEST FOR CYCLING. */
    if (*est <= estold)
    {
        goto L120;
    }
    i__1 = *n;
    for (i__ = 1;
            i__ <= i__1;
            ++i__)
    {
        x[i__] = d_sign(&c_b11, &x[i__]);
        isgn[i__] = i_dnnt(&x[i__]);
        /* L100: */
    }
    *kase = 2;
    jump = 4;
    return 0;
    /* ................ ENTRY (JUMP = 4) */
    /* X HAS BEEN OVERWRITTEN BY TRANSPOSE(A)*X. */
L110:
    jlast = j;
    j = idamax_(n, &x[1], &c__1);
    if (x[jlast] != (d__1 = x[j], f2c_abs(d__1)) && iter < 5)
    {
        ++iter;
        goto L50;
    }
    /* ITERATION COMPLETE. FINAL STAGE. */
L120:
    altsgn = 1.;
    i__1 = *n;
    for (i__ = 1;
            i__ <= i__1;
            ++i__)
    {
        x[i__] = altsgn * ((doublereal) (i__ - 1) / (doublereal) (*n - 1) + 1.);
        altsgn = -altsgn;
        /* L130: */
    }
    *kase = 1;
    jump = 5;
    return 0;
    /* ................ ENTRY (JUMP = 5) */
    /* X HAS BEEN OVERWRITTEN BY A*X. */
L140:
    temp = dasum_(n, &x[1], &c__1) / (doublereal) (*n * 3) * 2.;
    if (temp > *est)
    {
        dcopy_(n, &x[1], &c__1, &v[1], &c__1);
        *est = temp;
    }
L150:
    *kase = 0;
    return 0;
    /* End of DLACON */
}
/* dlacon_ */
