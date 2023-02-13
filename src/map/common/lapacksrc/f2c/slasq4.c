/* ../netlib/slasq4.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief \b SLASQ4 computes an approximation to the smallest eigenvalue using values of d from the previous transform. Used by sbdsqr. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download SLASQ4 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slasq4. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slasq4. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slasq4. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE SLASQ4( I0, N0, Z, PP, N0IN, DMIN, DMIN1, DMIN2, DN, */
/* DN1, DN2, TAU, TTYPE, G ) */
/* .. Scalar Arguments .. */
/* INTEGER I0, N0, N0IN, PP, TTYPE */
/* REAL DMIN, DMIN1, DMIN2, DN, DN1, DN2, G, TAU */
/* .. */
/* .. Array Arguments .. */
/* REAL Z( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLASQ4 computes an approximation TAU to the smallest eigenvalue */
/* > using values of d from the previous transform. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] I0 */
/* > \verbatim */
/* > I0 is INTEGER */
/* > First index. */
/* > \endverbatim */
/* > */
/* > \param[in] N0 */
/* > \verbatim */
/* > N0 is INTEGER */
/* > Last index. */
/* > \endverbatim */
/* > */
/* > \param[in] Z */
/* > \verbatim */
/* > Z is REAL array, dimension ( 4*N ) */
/* > Z holds the qd array. */
/* > \endverbatim */
/* > */
/* > \param[in] PP */
/* > \verbatim */
/* > PP is INTEGER */
/* > PP=0 for ping, PP=1 for pong. */
/* > \endverbatim */
/* > */
/* > \param[in] N0IN */
/* > \verbatim */
/* > N0IN is INTEGER */
/* > The value of N0 at start of EIGTEST. */
/* > \endverbatim */
/* > */
/* > \param[in] DMIN */
/* > \verbatim */
/* > DMIN is REAL */
/* > Minimum value of d. */
/* > \endverbatim */
/* > */
/* > \param[in] DMIN1 */
/* > \verbatim */
/* > DMIN1 is REAL */
/* > Minimum value of d, excluding D( N0 ). */
/* > \endverbatim */
/* > */
/* > \param[in] DMIN2 */
/* > \verbatim */
/* > DMIN2 is REAL */
/* > Minimum value of d, excluding D( N0 ) and D( N0-1 ). */
/* > \endverbatim */
/* > */
/* > \param[in] DN */
/* > \verbatim */
/* > DN is REAL */
/* > d(N) */
/* > \endverbatim */
/* > */
/* > \param[in] DN1 */
/* > \verbatim */
/* > DN1 is REAL */
/* > d(N-1) */
/* > \endverbatim */
/* > */
/* > \param[in] DN2 */
/* > \verbatim */
/* > DN2 is REAL */
/* > d(N-2) */
/* > \endverbatim */
/* > */
/* > \param[out] TAU */
/* > \verbatim */
/* > TAU is REAL */
/* > This is the shift. */
/* > \endverbatim */
/* > */
/* > \param[out] TTYPE */
/* > \verbatim */
/* > TTYPE is INTEGER */
/* > Shift type. */
/* > \endverbatim */
/* > */
/* > \param[in,out] G */
/* > \verbatim */
/* > G is REAL */
/* > G is passed as an argument in order to save its value between */
/* > calls to SLASQ4. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date September 2012 */
/* > \ingroup auxOTHERcomputational */
/* > \par Further Details: */
/* ===================== */
/* > */
/* > \verbatim */
/* > */
/* > CNST1 = 9/16 */
/* > \endverbatim */
/* > */
/* ===================================================================== */
/* Subroutine */
int slasq4_(integer *i0, integer *n0, real *z__, integer *pp, integer *n0in, real *dmin__, real *dmin1, real *dmin2, real *dn, real *dn1, real *dn2, real *tau, integer *ttype, real *g)
{
    /* System generated locals */
    integer i__1;
    real r__1, r__2;
    /* Builtin functions */
    double sqrt(doublereal);
    /* Local variables */
    real s, a2, b1, b2;
    integer i4, nn, np;
    real gam, gap1, gap2;
    /* -- LAPACK computational routine (version 3.4.2) -- */
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
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* A negative DMIN forces the shift to take that absolute value */
    /* TTYPE records the type of shift. */
    /* Parameter adjustments */
    --z__;
    /* Function Body */
    if (*dmin__ <= 0.f)
    {
        *tau = -(*dmin__);
        *ttype = -1;
        return 0;
    }
    nn = (*n0 << 2) + *pp;
    if (*n0in == *n0)
    {
        /* No eigenvalues deflated. */
        if (*dmin__ == *dn || *dmin__ == *dn1)
        {
            b1 = sqrt(z__[nn - 3]) * sqrt(z__[nn - 5]);
            b2 = sqrt(z__[nn - 7]) * sqrt(z__[nn - 9]);
            a2 = z__[nn - 7] + z__[nn - 5];
            /* Cases 2 and 3. */
            if (*dmin__ == *dn && *dmin1 == *dn1)
            {
                gap2 = *dmin2 - a2 - *dmin2 * .25f;
                if (gap2 > 0.f && gap2 > b2)
                {
                    gap1 = a2 - *dn - b2 / gap2 * b2;
                }
                else
                {
                    gap1 = a2 - *dn - (b1 + b2);
                }
                if (gap1 > 0.f && gap1 > b1)
                {
                    /* Computing MAX */
                    r__1 = *dn - b1 / gap1 * b1;
                    r__2 = *dmin__ * .5f; // , expr subst
                    s = max(r__1,r__2);
                    *ttype = -2;
                }
                else
                {
                    s = 0.f;
                    if (*dn > b1)
                    {
                        s = *dn - b1;
                    }
                    if (a2 > b1 + b2)
                    {
                        /* Computing MIN */
                        r__1 = s;
                        r__2 = a2 - (b1 + b2); // , expr subst
                        s = min(r__1,r__2);
                    }
                    /* Computing MAX */
                    r__1 = s;
                    r__2 = *dmin__ * .333f; // , expr subst
                    s = max(r__1,r__2);
                    *ttype = -3;
                }
            }
            else
            {
                /* Case 4. */
                *ttype = -4;
                s = *dmin__ * .25f;
                if (*dmin__ == *dn)
                {
                    gam = *dn;
                    a2 = 0.f;
                    if (z__[nn - 5] > z__[nn - 7])
                    {
                        return 0;
                    }
                    b2 = z__[nn - 5] / z__[nn - 7];
                    np = nn - 9;
                }
                else
                {
                    np = nn - (*pp << 1);
                    b2 = z__[np - 2];
                    gam = *dn1;
                    if (z__[np - 4] > z__[np - 2])
                    {
                        return 0;
                    }
                    a2 = z__[np - 4] / z__[np - 2];
                    if (z__[nn - 9] > z__[nn - 11])
                    {
                        return 0;
                    }
                    b2 = z__[nn - 9] / z__[nn - 11];
                    np = nn - 13;
                }
                /* Approximate contribution to norm squared from I < NN-1. */
                a2 += b2;
                i__1 = (*i0 << 2) - 1 + *pp;
                for (i4 = np;
                        i4 >= i__1;
                        i4 += -4)
                {
                    if (b2 == 0.f)
                    {
                        goto L20;
                    }
                    b1 = b2;
                    if (z__[i4] > z__[i4 - 2])
                    {
                        return 0;
                    }
                    b2 *= z__[i4] / z__[i4 - 2];
                    a2 += b2;
                    if (max(b2,b1) * 100.f < a2 || .563f < a2)
                    {
                        goto L20;
                    }
                    /* L10: */
                }
L20:
                a2 *= 1.05f;
                /* Rayleigh quotient residual bound. */
                if (a2 < .563f)
                {
                    s = gam * (1.f - sqrt(a2)) / (a2 + 1.f);
                }
            }
        }
        else if (*dmin__ == *dn2)
        {
            /* Case 5. */
            *ttype = -5;
            s = *dmin__ * .25f;
            /* Compute contribution to norm squared from I > NN-2. */
            np = nn - (*pp << 1);
            b1 = z__[np - 2];
            b2 = z__[np - 6];
            gam = *dn2;
            if (z__[np - 8] > b2 || z__[np - 4] > b1)
            {
                return 0;
            }
            a2 = z__[np - 8] / b2 * (z__[np - 4] / b1 + 1.f);
            /* Approximate contribution to norm squared from I < NN-2. */
            if (*n0 - *i0 > 2)
            {
                b2 = z__[nn - 13] / z__[nn - 15];
                a2 += b2;
                i__1 = (*i0 << 2) - 1 + *pp;
                for (i4 = nn - 17;
                        i4 >= i__1;
                        i4 += -4)
                {
                    if (b2 == 0.f)
                    {
                        goto L40;
                    }
                    b1 = b2;
                    if (z__[i4] > z__[i4 - 2])
                    {
                        return 0;
                    }
                    b2 *= z__[i4] / z__[i4 - 2];
                    a2 += b2;
                    if (max(b2,b1) * 100.f < a2 || .563f < a2)
                    {
                        goto L40;
                    }
                    /* L30: */
                }
L40:
                a2 *= 1.05f;
            }
            if (a2 < .563f)
            {
                s = gam * (1.f - sqrt(a2)) / (a2 + 1.f);
            }
        }
        else
        {
            /* Case 6, no information to guide us. */
            if (*ttype == -6)
            {
                *g += (1.f - *g) * .333f;
            }
            else if (*ttype == -18)
            {
                *g = .083250000000000005f;
            }
            else
            {
                *g = .25f;
            }
            s = *g * *dmin__;
            *ttype = -6;
        }
    }
    else if (*n0in == *n0 + 1)
    {
        /* One eigenvalue just deflated. Use DMIN1, DN1 for DMIN and DN. */
        if (*dmin1 == *dn1 && *dmin2 == *dn2)
        {
            /* Cases 7 and 8. */
            *ttype = -7;
            s = *dmin1 * .333f;
            if (z__[nn - 5] > z__[nn - 7])
            {
                return 0;
            }
            b1 = z__[nn - 5] / z__[nn - 7];
            b2 = b1;
            if (b2 == 0.f)
            {
                goto L60;
            }
            i__1 = (*i0 << 2) - 1 + *pp;
            for (i4 = (*n0 << 2) - 9 + *pp;
                    i4 >= i__1;
                    i4 += -4)
            {
                a2 = b1;
                if (z__[i4] > z__[i4 - 2])
                {
                    return 0;
                }
                b1 *= z__[i4] / z__[i4 - 2];
                b2 += b1;
                if (max(b1,a2) * 100.f < b2)
                {
                    goto L60;
                }
                /* L50: */
            }
L60:
            b2 = sqrt(b2 * 1.05f);
            /* Computing 2nd power */
            r__1 = b2;
            a2 = *dmin1 / (r__1 * r__1 + 1.f);
            gap2 = *dmin2 * .5f - a2;
            if (gap2 > 0.f && gap2 > b2 * a2)
            {
                /* Computing MAX */
                r__1 = s;
                r__2 = a2 * (1.f - a2 * 1.01f * (b2 / gap2) * b2); // , expr subst
                s = max(r__1,r__2);
            }
            else
            {
                /* Computing MAX */
                r__1 = s;
                r__2 = a2 * (1.f - b2 * 1.01f); // , expr subst
                s = max(r__1,r__2);
                *ttype = -8;
            }
        }
        else
        {
            /* Case 9. */
            s = *dmin1 * .25f;
            if (*dmin1 == *dn1)
            {
                s = *dmin1 * .5f;
            }
            *ttype = -9;
        }
    }
    else if (*n0in == *n0 + 2)
    {
        /* Two eigenvalues deflated. Use DMIN2, DN2 for DMIN and DN. */
        /* Cases 10 and 11. */
        if (*dmin2 == *dn2 && z__[nn - 5] * 2.f < z__[nn - 7])
        {
            *ttype = -10;
            s = *dmin2 * .333f;
            if (z__[nn - 5] > z__[nn - 7])
            {
                return 0;
            }
            b1 = z__[nn - 5] / z__[nn - 7];
            b2 = b1;
            if (b2 == 0.f)
            {
                goto L80;
            }
            i__1 = (*i0 << 2) - 1 + *pp;
            for (i4 = (*n0 << 2) - 9 + *pp;
                    i4 >= i__1;
                    i4 += -4)
            {
                if (z__[i4] > z__[i4 - 2])
                {
                    return 0;
                }
                b1 *= z__[i4] / z__[i4 - 2];
                b2 += b1;
                if (b1 * 100.f < b2)
                {
                    goto L80;
                }
                /* L70: */
            }
L80:
            b2 = sqrt(b2 * 1.05f);
            /* Computing 2nd power */
            r__1 = b2;
            a2 = *dmin2 / (r__1 * r__1 + 1.f);
            gap2 = z__[nn - 7] + z__[nn - 9] - sqrt(z__[nn - 11]) * sqrt(z__[ nn - 9]) - a2;
            if (gap2 > 0.f && gap2 > b2 * a2)
            {
                /* Computing MAX */
                r__1 = s;
                r__2 = a2 * (1.f - a2 * 1.01f * (b2 / gap2) * b2); // , expr subst
                s = max(r__1,r__2);
            }
            else
            {
                /* Computing MAX */
                r__1 = s;
                r__2 = a2 * (1.f - b2 * 1.01f); // , expr subst
                s = max(r__1,r__2);
            }
        }
        else
        {
            s = *dmin2 * .25f;
            *ttype = -11;
        }
    }
    else if (*n0in > *n0 + 2)
    {
        /* Case 12, more than two eigenvalues deflated. No information. */
        s = 0.f;
        *ttype = -12;
    }
    *tau = s;
    return 0;
    /* End of SLASQ4 */
}
/* slasq4_ */
