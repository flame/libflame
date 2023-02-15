/* ../netlib/zlacn2.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
/* > \brief \b ZLACN2 estimates the 1-norm of a square matrix, using reverse communication for evaluating matr ix-vector products. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download ZLACN2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlacn2. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlacn2. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlacn2. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE ZLACN2( N, V, X, EST, KASE, ISAVE ) */
/* .. Scalar Arguments .. */
/* INTEGER KASE, N */
/* DOUBLE PRECISION EST */
/* .. */
/* .. Array Arguments .. */
/* INTEGER ISAVE( 3 ) */
/* COMPLEX*16 V( * ), X( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZLACN2 estimates the 1-norm of a square, complex matrix A. */
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
/* > V is COMPLEX*16 array, dimension (N) */
/* > On the final return, V = A*W, where EST = norm(V)/norm(W) */
/* > (W is not returned). */
/* > \endverbatim */
/* > */
/* > \param[in,out] X */
/* > \verbatim */
/* > X is COMPLEX*16 array, dimension (N) */
/* > On an intermediate return, X should be overwritten by */
/* > A * X, if KASE=1, */
/* > A**H * X, if KASE=2, */
/* > where A**H is the conjugate transpose of A, and ZLACN2 must be */
/* > re-called with all the other parameters unchanged. */
/* > \endverbatim */
/* > */
/* > \param[in,out] EST */
/* > \verbatim */
/* > EST is DOUBLE PRECISION */
/* > On entry with KASE = 1 or 2 and ISAVE(1) = 3, EST should be */
/* > unchanged from the previous call to ZLACN2. */
/* > On exit, EST is an estimate (a lower bound) for norm(A). */
/* > \endverbatim */
/* > */
/* > \param[in,out] KASE */
/* > \verbatim */
/* > KASE is INTEGER */
/* > On the initial call to ZLACN2, KASE should be 0. */
/* > On an intermediate return, KASE will be 1 or 2, indicating */
/* > whether X should be overwritten by A * X or A**H * X. */
/* > On the final return from ZLACN2, KASE will again be 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] ISAVE */
/* > \verbatim */
/* > ISAVE is INTEGER array, dimension (3) */
/* > ISAVE is used to save variables between calls to ZLACN2 */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date September 2012 */
/* > \ingroup complex16OTHERauxiliary */
/* > \par Further Details: */
/* ===================== */
/* > */
/* > \verbatim */
/* > */
/* > Originally named CONEST, dated March 16, 1988. */
/* > */
/* > Last modified: April, 1999 */
/* > */
/* > This is a thread safe version of ZLACON, which uses the array ISAVE */
/* > in place of a SAVE statement, as follows: */
/* > */
/* > ZLACON ZLACN2 */
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
int zlacn2_(integer *n, doublecomplex *v, doublecomplex *x, doublereal *est, integer *kase, integer *isave)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1, d__2;
    doublecomplex z__1;
    /* Builtin functions */
    double z_abs(doublecomplex *), d_imag(doublecomplex *);
    /* Local variables */
    integer i__;
    doublereal temp, absxi;
    integer jlast;
    extern /* Subroutine */
    int zcopy_(integer *, doublecomplex *, integer *, doublecomplex *, integer *);
    extern integer izmax1_(integer *, doublecomplex *, integer *);
    extern doublereal dzsum1_(integer *, doublecomplex *, integer *), dlamch_( char *);
    doublereal safmin, altsgn, estold;
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
    --x;
    --v;
    /* Function Body */
    safmin = dlamch_("Safe minimum");
    if (*kase == 0)
    {
        i__1 = *n;
        for (i__ = 1;
                i__ <= i__1;
                ++i__)
        {
            i__2 = i__;
            d__1 = 1. / (doublereal) (*n);
            z__1.r = d__1;
            z__1.i = 0.; // , expr subst
            x[i__2].r = z__1.r;
            x[i__2].i = z__1.i; // , expr subst
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
        goto L90;
    case 5:
        goto L120;
    }
    /* ................ ENTRY (ISAVE( 1 ) = 1) */
    /* FIRST ITERATION. X HAS BEEN OVERWRITTEN BY A*X. */
L20:
    if (*n == 1)
    {
        v[1].r = x[1].r;
        v[1].i = x[1].i; // , expr subst
        *est = z_abs(&v[1]);
        /* ... QUIT */
        goto L130;
    }
    *est = dzsum1_(n, &x[1], &c__1);
    i__1 = *n;
    for (i__ = 1;
            i__ <= i__1;
            ++i__)
    {
        absxi = z_abs(&x[i__]);
        if (absxi > safmin)
        {
            i__2 = i__;
            i__3 = i__;
            d__1 = x[i__3].r / absxi;
            d__2 = d_imag(&x[i__]) / absxi;
            z__1.r = d__1;
            z__1.i = d__2; // , expr subst
            x[i__2].r = z__1.r;
            x[i__2].i = z__1.i; // , expr subst
        }
        else
        {
            i__2 = i__;
            x[i__2].r = 1.;
            x[i__2].i = 0.; // , expr subst
        }
        /* L30: */
    }
    *kase = 2;
    isave[1] = 2;
    return 0;
    /* ................ ENTRY (ISAVE( 1 ) = 2) */
    /* FIRST ITERATION. X HAS BEEN OVERWRITTEN BY CTRANS(A)*X. */
L40:
    isave[2] = izmax1_(n, &x[1], &c__1);
    isave[3] = 2;
    /* MAIN LOOP - ITERATIONS 2,3,...,ITMAX. */
L50:
    i__1 = *n;
    for (i__ = 1;
            i__ <= i__1;
            ++i__)
    {
        i__2 = i__;
        x[i__2].r = 0.;
        x[i__2].i = 0.; // , expr subst
        /* L60: */
    }
    i__1 = isave[2];
    x[i__1].r = 1.;
    x[i__1].i = 0.; // , expr subst
    *kase = 1;
    isave[1] = 3;
    return 0;
    /* ................ ENTRY (ISAVE( 1 ) = 3) */
    /* X HAS BEEN OVERWRITTEN BY A*X. */
L70:
    zcopy_(n, &x[1], &c__1, &v[1], &c__1);
    estold = *est;
    *est = dzsum1_(n, &v[1], &c__1);
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
        absxi = z_abs(&x[i__]);
        if (absxi > safmin)
        {
            i__2 = i__;
            i__3 = i__;
            d__1 = x[i__3].r / absxi;
            d__2 = d_imag(&x[i__]) / absxi;
            z__1.r = d__1;
            z__1.i = d__2; // , expr subst
            x[i__2].r = z__1.r;
            x[i__2].i = z__1.i; // , expr subst
        }
        else
        {
            i__2 = i__;
            x[i__2].r = 1.;
            x[i__2].i = 0.; // , expr subst
        }
        /* L80: */
    }
    *kase = 2;
    isave[1] = 4;
    return 0;
    /* ................ ENTRY (ISAVE( 1 ) = 4) */
    /* X HAS BEEN OVERWRITTEN BY CTRANS(A)*X. */
L90:
    jlast = isave[2];
    isave[2] = izmax1_(n, &x[1], &c__1);
    if (z_abs(&x[jlast]) != z_abs(&x[isave[2]]) && isave[3] < 5)
    {
        ++isave[3];
        goto L50;
    }
    /* ITERATION COMPLETE. FINAL STAGE. */
L100:
    altsgn = 1.;
    i__1 = *n;
    for (i__ = 1;
            i__ <= i__1;
            ++i__)
    {
        i__2 = i__;
        d__1 = altsgn * ((doublereal) (i__ - 1) / (doublereal) (*n - 1) + 1.);
        z__1.r = d__1;
        z__1.i = 0.; // , expr subst
        x[i__2].r = z__1.r;
        x[i__2].i = z__1.i; // , expr subst
        altsgn = -altsgn;
        /* L110: */
    }
    *kase = 1;
    isave[1] = 5;
    return 0;
    /* ................ ENTRY (ISAVE( 1 ) = 5) */
    /* X HAS BEEN OVERWRITTEN BY A*X. */
L120:
    temp = dzsum1_(n, &x[1], &c__1) / (doublereal) (*n * 3) * 2.;
    if (temp > *est)
    {
        zcopy_(n, &x[1], &c__1, &v[1], &c__1);
        *est = temp;
    }
L130:
    *kase = 0;
    return 0;
    /* End of ZLACN2 */
}
/* zlacn2_ */
