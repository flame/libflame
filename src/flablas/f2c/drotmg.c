/* drotmg.f -- translated by f2c (version 19991025). You must link the resulting object file with the libraries: -lf2c -lm (in that order) */
#include "FLA_f2c.h"
/* Subroutine */
int drotmg_(doublereal *dd1, doublereal *dd2, doublereal * dx1, doublereal *dy1, doublereal *dparam)
{
    /* Initialized data */
    static const doublereal zero = 0.;
    static const doublereal one = 1.;
    static const doublereal two = 2.;
    static const doublereal gam = 4096.;
    static const doublereal gamsq = 16777216.;
    static const doublereal rgamsq = 5.9604645e-8;
    /* Format strings */
    char fmt_120[] = "";
    char fmt_150[] = "";
    char fmt_180[] = "";
    char fmt_210[] = "";
    /* System generated locals */
    doublereal d__1;
    /* Local variables */
    doublereal dflag, dtemp, du, dp1, dp2, dq2, dq1, dh11, dh21, dh12, dh22;
    integer igo;
    /* Assigned format variables */
    char *igo_fmt;
    /* CONSTRUCT THE MODIFIED GIVENS TRANSFORMATION MATRIX H WHICH ZEROS */
    /* THE SECOND COMPONENT OF THE 2-VECTOR (DSQRT(DD1)*DX1,DSQRT(DD2)* */
    /* DY2)**T. */
    /* WITH DPARAM(1)=DFLAG, H HAS ONE OF THE FOLLOWING FORMS.. */
    /* DFLAG=-1.D0 DFLAG=0.D0 DFLAG=1.D0 DFLAG=-2.D0 */
    /* (DH11 DH12) (1.D0 DH12) (DH11 1.D0) (1.D0 0.D0) */
    /* H=( ) ( ) ( ) ( ) */
    /* (DH21 DH22), (DH21 1.D0), (-1.D0 DH22), (0.D0 1.D0). */
    /* LOCATIONS 2-4 OF DPARAM CONTAIN DH11, DH21, DH12, AND DH22 */
    /* RESPECTIVELY. (VALUES OF 1.D0, -1.D0, OR 0.D0 IMPLIED BY THE */
    /* VALUE OF DPARAM(1) ARE NOT STORED IN DPARAM.) */
    /* THE VALUES OF GAMSQ AND RGAMSQ SET IN THE DATA STATEMENT MAY BE */
    /* INEXACT. THIS IS OK AS THEY ARE ONLY USED FOR TESTING THE SIZE */
    /* OF DD1 AND DD2. ALL ACTUAL SCALING OF DATA IS DONE USING GAM. */
    /* Parameter adjustments */
    --dparam;
    /* Function Body */
    if (! (*dd1 < zero))
    {
        goto L10;
    }
    /* GO ZERO-H-D-AND-DX1.. */
    goto L60;
L10: /* CASE-DD1-NONNEGATIVE */
    dp2 = *dd2 * *dy1;
    if (! (dp2 == zero))
    {
        goto L20;
    }
    dflag = -two;
    goto L260;
    /* REGULAR-CASE.. */
L20:
    dp1 = *dd1 * *dx1;
    dq2 = dp2 * *dy1;
    dq1 = dp1 * *dx1;
    if (! (f2c_abs(dq1) > f2c_abs(dq2)))
    {
        goto L40;
    }
    dh21 = -(*dy1) / *dx1;
    dh12 = dp2 / dp1;
    du = one - dh12 * dh21;
    if (! (du <= zero))
    {
        goto L30;
    }
    /* GO ZERO-H-D-AND-DX1.. */
    goto L60;
L30:
    dflag = zero;
    *dd1 /= du;
    *dd2 /= du;
    *dx1 *= du;
    /* GO SCALE-CHECK.. */
    goto L100;
L40:
    if (! (dq2 < zero))
    {
        goto L50;
    }
    /* GO ZERO-H-D-AND-DX1.. */
    goto L60;
L50:
    dflag = one;
    dh11 = dp1 / dp2;
    dh22 = *dx1 / *dy1;
    du = one + dh11 * dh22;
    dtemp = *dd2 / du;
    *dd2 = *dd1 / du;
    *dd1 = dtemp;
    *dx1 = *dy1 * du;
    /* GO SCALE-CHECK */
    goto L100;
    /* PROCEDURE..ZERO-H-D-AND-DX1.. */
L60:
    dflag = -one;
    dh11 = zero;
    dh12 = zero;
    dh21 = zero;
    dh22 = zero;
    *dd1 = zero;
    *dd2 = zero;
    *dx1 = zero;
    /* RETURN.. */
    goto L220;
    /* PROCEDURE..FIX-H.. */
L70:
    if (! (dflag >= zero))
    {
        goto L90;
    }
    if (! (dflag == zero))
    {
        goto L80;
    }
    dh11 = one;
    dh22 = one;
    dflag = -one;
    goto L90;
L80:
    dh21 = -one;
    dh12 = one;
    dflag = -one;
L90:
    switch (igo)
    {
    case 0:
        goto L120;
    case 1:
        goto L150;
    case 2:
        goto L180;
    case 3:
        goto L210;
    }
    /* PROCEDURE..SCALE-CHECK */
L100:
L110:
    if (! (*dd1 <= rgamsq))
    {
        goto L130;
    }
    if (*dd1 == zero)
    {
        goto L160;
    }
    igo = 0;
    igo_fmt = fmt_120;
    /* FIX-H.. */
    goto L70;
L120: /* Computing 2nd power */
    d__1 = gam;
    *dd1 *= d__1 * d__1;
    *dx1 /= gam;
    dh11 /= gam;
    dh12 /= gam;
    goto L110;
L130:
L140:
    if (! (*dd1 >= gamsq))
    {
        goto L160;
    }
    igo = 1;
    igo_fmt = fmt_150;
    /* FIX-H.. */
    goto L70;
L150: /* Computing 2nd power */
    d__1 = gam;
    *dd1 /= d__1 * d__1;
    *dx1 *= gam;
    dh11 *= gam;
    dh12 *= gam;
    goto L140;
L160:
L170:
    if (! (f2c_abs(*dd2) <= rgamsq))
    {
        goto L190;
    }
    if (*dd2 == zero)
    {
        goto L220;
    }
    igo = 2;
    igo_fmt = fmt_180;
    /* FIX-H.. */
    goto L70;
L180: /* Computing 2nd power */
    d__1 = gam;
    *dd2 *= d__1 * d__1;
    dh21 /= gam;
    dh22 /= gam;
    goto L170;
L190:
L200:
    if (! (f2c_abs(*dd2) >= gamsq))
    {
        goto L220;
    }
    igo = 3;
    igo_fmt = fmt_210;
    /* FIX-H.. */
    goto L70;
L210: /* Computing 2nd power */
    d__1 = gam;
    *dd2 /= d__1 * d__1;
    dh21 *= gam;
    dh22 *= gam;
    goto L200;
L220:
    if (dflag < 0.)
    {
        goto L250;
    }
    else if (dflag == 0)
    {
        goto L230;
    }
    else
    {
        goto L240;
    }
L230:
    dparam[3] = dh21;
    dparam[4] = dh12;
    goto L260;
L240:
    dparam[2] = dh11;
    dparam[5] = dh22;
    goto L260;
L250:
    dparam[2] = dh11;
    dparam[3] = dh21;
    dparam[4] = dh12;
    dparam[5] = dh22;
L260:
    dparam[1] = dflag;
    return 0;
}
/* drotmg_ */

