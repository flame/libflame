/* ../netlib/zlarnv.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief \b ZLARNV returns a vector of random numbers from a uniform or normal distribution. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download ZLARNV + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlarnv. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlarnv. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlarnv. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE ZLARNV( IDIST, ISEED, N, X ) */
/* .. Scalar Arguments .. */
/* INTEGER IDIST, N */
/* .. */
/* .. Array Arguments .. */
/* INTEGER ISEED( 4 ) */
/* COMPLEX*16 X( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZLARNV returns a vector of n random complex numbers from a uniform or */
/* > normal distribution. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] IDIST */
/* > \verbatim */
/* > IDIST is INTEGER */
/* > Specifies the distribution of the random numbers: */
/* > = 1: real and imaginary parts each uniform (0,1) */
/* > = 2: real and imaginary parts each uniform (-1,1) */
/* > = 3: real and imaginary parts each normal (0,1) */
/* > = 4: uniformly distributed on the disc f2c_abs(z) < 1 */
/* > = 5: uniformly distributed on the circle f2c_abs(z) = 1 */
/* > \endverbatim */
/* > */
/* > \param[in,out] ISEED */
/* > \verbatim */
/* > ISEED is INTEGER array, dimension (4) */
/* > On entry, the seed of the random number generator;
the array */
/* > elements must be between 0 and 4095, and ISEED(4) must be */
/* > odd. */
/* > On exit, the seed is updated. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The number of random numbers to be generated. */
/* > \endverbatim */
/* > */
/* > \param[out] X */
/* > \verbatim */
/* > X is COMPLEX*16 array, dimension (N) */
/* > The generated random numbers. */
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
/* > This routine calls the auxiliary routine DLARUV to generate random */
/* > real numbers from a uniform (0,1) distribution, in batches of up to */
/* > 128 using vectorisable code. The Box-Muller method is used to */
/* > transform numbers from a uniform to a normal distribution. */
/* > \endverbatim */
/* > */
/* ===================================================================== */
/* Subroutine */
int zlarnv_(integer *idist, integer *iseed, integer *n, doublecomplex *x)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5;
    doublereal d__1, d__2;
    doublecomplex z__1, z__2, z__3;
    /* Builtin functions */
    double log(doublereal), sqrt(doublereal);
    void z_exp(doublecomplex *, doublecomplex *);
    /* Local variables */
    integer i__;
    doublereal u[128];
    integer il, iv;
    extern /* Subroutine */
    int dlaruv_(integer *, integer *, doublereal *);
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
    /* .. Local Arrays .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. External Subroutines .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Parameter adjustments */
    --x;
    --iseed;
    /* Function Body */
    i__1 = *n;
    for (iv = 1;
            iv <= i__1;
            iv += 64)
    {
        /* Computing MIN */
        i__2 = 64;
        i__3 = *n - iv + 1; // , expr subst
        il = min(i__2,i__3);
        /* Call DLARUV to generate 2*IL real numbers from a uniform (0,1) */
        /* distribution (2*IL <= LV) */
        i__2 = il << 1;
        dlaruv_(&iseed[1], &i__2, u);
        if (*idist == 1)
        {
            /* Copy generated numbers */
            i__2 = il;
            for (i__ = 1;
                    i__ <= i__2;
                    ++i__)
            {
                i__3 = iv + i__ - 1;
                i__4 = (i__ << 1) - 2;
                i__5 = (i__ << 1) - 1;
                z__1.r = u[i__4];
                z__1.i = u[i__5]; // , expr subst
                x[i__3].r = z__1.r;
                x[i__3].i = z__1.i; // , expr subst
                /* L10: */
            }
        }
        else if (*idist == 2)
        {
            /* Convert generated numbers to uniform (-1,1) distribution */
            i__2 = il;
            for (i__ = 1;
                    i__ <= i__2;
                    ++i__)
            {
                i__3 = iv + i__ - 1;
                d__1 = u[(i__ << 1) - 2] * 2. - 1.;
                d__2 = u[(i__ << 1) - 1] * 2. - 1.;
                z__1.r = d__1;
                z__1.i = d__2; // , expr subst
                x[i__3].r = z__1.r;
                x[i__3].i = z__1.i; // , expr subst
                /* L20: */
            }
        }
        else if (*idist == 3)
        {
            /* Convert generated numbers to normal (0,1) distribution */
            i__2 = il;
            for (i__ = 1;
                    i__ <= i__2;
                    ++i__)
            {
                i__3 = iv + i__ - 1;
                d__1 = sqrt(log(u[(i__ << 1) - 2]) * -2.);
                d__2 = u[(i__ << 1) - 1] * 6.2831853071795864769252867663;
                z__3.r = 0.;
                z__3.i = d__2; // , expr subst
                z_exp(&z__2, &z__3);
                z__1.r = d__1 * z__2.r;
                z__1.i = d__1 * z__2.i; // , expr subst
                x[i__3].r = z__1.r;
                x[i__3].i = z__1.i; // , expr subst
                /* L30: */
            }
        }
        else if (*idist == 4)
        {
            /* Convert generated numbers to complex numbers uniformly */
            /* distributed on the unit disk */
            i__2 = il;
            for (i__ = 1;
                    i__ <= i__2;
                    ++i__)
            {
                i__3 = iv + i__ - 1;
                d__1 = sqrt(u[(i__ << 1) - 2]);
                d__2 = u[(i__ << 1) - 1] * 6.2831853071795864769252867663;
                z__3.r = 0.;
                z__3.i = d__2; // , expr subst
                z_exp(&z__2, &z__3);
                z__1.r = d__1 * z__2.r;
                z__1.i = d__1 * z__2.i; // , expr subst
                x[i__3].r = z__1.r;
                x[i__3].i = z__1.i; // , expr subst
                /* L40: */
            }
        }
        else if (*idist == 5)
        {
            /* Convert generated numbers to complex numbers uniformly */
            /* distributed on the unit circle */
            i__2 = il;
            for (i__ = 1;
                    i__ <= i__2;
                    ++i__)
            {
                i__3 = iv + i__ - 1;
                d__1 = u[(i__ << 1) - 1] * 6.2831853071795864769252867663;
                z__2.r = 0.;
                z__2.i = d__1; // , expr subst
                z_exp(&z__1, &z__2);
                x[i__3].r = z__1.r;
                x[i__3].i = z__1.i; // , expr subst
                /* L50: */
            }
        }
        /* L60: */
    }
    return 0;
    /* End of ZLARNV */
}
/* zlarnv_ */
