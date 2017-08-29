/* ../netlib/slacn2.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
static real c_b11 = 1.f;
/* > \brief \b SLACN2 estimates the 1-norm of a square matrix, using reverse communication for evaluating matr ix-vector products. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download SLACN2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slacn2. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slacn2. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slacn2. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE SLACN2( N, V, X, ISGN, EST, KASE, ISAVE ) */
/* .. Scalar Arguments .. */
/* INTEGER KASE, N */
/* REAL EST */
/* .. */
/* .. Array Arguments .. */
/* INTEGER ISGN( * ), ISAVE( 3 ) */
/* REAL V( * ), X( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLACN2 estimates the 1-norm of a square, real matrix A. */
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
/* > V is REAL array, dimension (N) */
/* > On the final return, V = A*W, where EST = norm(V)/norm(W) */
/* > (W is not returned). */
/* > \endverbatim */
/* > */
/* > \param[in,out] X */
/* > \verbatim */
/* > X is REAL array, dimension (N) */
/* > On an intermediate return, X should be overwritten by */
/* > A * X, if KASE=1, */
/* > A**T * X, if KASE=2, */
/* > and SLACN2 must be re-called with all the other parameters */
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
/* > EST is REAL */
/* > On entry with KASE = 1 or 2 and ISAVE(1) = 3, EST should be */
/* > unchanged from the previous call to SLACN2. */
/* > On exit, EST is an estimate (a lower bound) for norm(A). */
/* > \endverbatim */
/* > */
/* > \param[in,out] KASE */
/* > \verbatim */
/* > KASE is INTEGER */
/* > On the initial call to SLACN2, KASE should be 0. */
/* > On an intermediate return, KASE will be 1 or 2, indicating */
/* > whether X should be overwritten by A * X or A**T * X. */
/* > On the final return from SLACN2, KASE will again be 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] ISAVE */
/* > \verbatim */
/* > ISAVE is INTEGER array, dimension (3) */
/* > ISAVE is used to save variables between calls to SLACN2 */
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
/* > Originally named SONEST, dated March 16, 1988. */
/* > */
/* > This is a thread safe version of SLACON, which uses the array ISAVE */
/* > in place of a SAVE statement, as follows: */
/* > */
/* > SLACON SLACN2 */
/* > JUMP ISAVE(1) */
/* > J ISAVE(2) */
/* > ITER ISAVE(3) */
/* > \endverbatim */
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
int slacn2_(integer *n, real *v, real *x, integer *isgn, real *est, integer *kase, integer *isave)
{
    /* System generated locals */
    integer i__1;
    real r__1;
    /* Builtin functions */
    double r_sign(real *, real *);
    integer i_nint(real *);
    /* Local variables */
    integer i__;
    real temp;
    integer jlast;
    extern real sasum_(integer *, real *, integer *);
    extern /* Subroutine */
    int scopy_(integer *, real *, integer *, real *, integer *);
    extern integer isamax_(integer *, real *, integer *);
    real altsgn, estold;
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
    /* Parameter adjustments */
    --isave;
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
            x[i__] = 1.f / (real) (*n);
            /* L10: */
        }
        *kase = 1;
        isave[1] = 1;
        return 0;
    }
    switch (isave[1])
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
    /* ................ ENTRY (ISAVE( 1 ) = 1) */
    /* FIRST ITERATION. X HAS BEEN OVERWRITTEN BY A*X. */
L20:
    if (*n == 1)
    {
        v[1] = x[1];
        *est = f2c_abs(v[1]);
        /* ... QUIT */
        goto L150;
    }
    *est = sasum_(n, &x[1], &c__1);
    i__1 = *n;
    for (i__ = 1;
            i__ <= i__1;
            ++i__)
    {
        x[i__] = r_sign(&c_b11, &x[i__]);
        isgn[i__] = i_nint(&x[i__]);
        /* L30: */
    }
    *kase = 2;
    isave[1] = 2;
    return 0;
    /* ................ ENTRY (ISAVE( 1 ) = 2) */
    /* FIRST ITERATION. X HAS BEEN OVERWRITTEN BY TRANSPOSE(A)*X. */
L40:
    isave[2] = isamax_(n, &x[1], &c__1);
    isave[3] = 2;
    /* MAIN LOOP - ITERATIONS 2,3,...,ITMAX. */
L50:
    i__1 = *n;
    for (i__ = 1;
            i__ <= i__1;
            ++i__)
    {
        x[i__] = 0.f;
        /* L60: */
    }
    x[isave[2]] = 1.f;
    *kase = 1;
    isave[1] = 3;
    return 0;
    /* ................ ENTRY (ISAVE( 1 ) = 3) */
    /* X HAS BEEN OVERWRITTEN BY A*X. */
L70:
    scopy_(n, &x[1], &c__1, &v[1], &c__1);
    estold = *est;
    *est = sasum_(n, &v[1], &c__1);
    i__1 = *n;
    for (i__ = 1;
            i__ <= i__1;
            ++i__)
    {
        r__1 = r_sign(&c_b11, &x[i__]);
        if (i_nint(&r__1) != isgn[i__])
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
        x[i__] = r_sign(&c_b11, &x[i__]);
        isgn[i__] = i_nint(&x[i__]);
        /* L100: */
    }
    *kase = 2;
    isave[1] = 4;
    return 0;
    /* ................ ENTRY (ISAVE( 1 ) = 4) */
    /* X HAS BEEN OVERWRITTEN BY TRANSPOSE(A)*X. */
L110:
    jlast = isave[2];
    isave[2] = isamax_(n, &x[1], &c__1);
    if (x[jlast] != (r__1 = x[isave[2]], f2c_abs(r__1)) && isave[3] < 5)
    {
        ++isave[3];
        goto L50;
    }
    /* ITERATION COMPLETE. FINAL STAGE. */
L120:
    altsgn = 1.f;
    i__1 = *n;
    for (i__ = 1;
            i__ <= i__1;
            ++i__)
    {
        x[i__] = altsgn * ((real) (i__ - 1) / (real) (*n - 1) + 1.f);
        altsgn = -altsgn;
        /* L130: */
    }
    *kase = 1;
    isave[1] = 5;
    return 0;
    /* ................ ENTRY (ISAVE( 1 ) = 5) */
    /* X HAS BEEN OVERWRITTEN BY A*X. */
L140:
    temp = sasum_(n, &x[1], &c__1) / (real) (*n * 3) * 2.f;
    if (temp > *est)
    {
        scopy_(n, &x[1], &c__1, &v[1], &c__1);
        *est = temp;
    }
L150:
    *kase = 0;
    return 0;
    /* End of SLACN2 */
}
/* slacn2_ */
