/* ../netlib/clarfg.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static complex c_b5 =
{
    1.f,0.f
}
;
/* > \brief \b CLARFG generates an elementary reflector (Householder matrix). */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CLARFG + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clarfg. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clarfg. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clarfg. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CLARFG( N, ALPHA, X, INCX, TAU ) */
/* .. Scalar Arguments .. */
/* INTEGER INCX, N */
/* COMPLEX ALPHA, TAU */
/* .. */
/* .. Array Arguments .. */
/* COMPLEX X( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CLARFG generates a complex elementary reflector H of order n, such */
/* > that */
/* > */
/* > H**H * ( alpha ) = ( beta ), H**H * H = I. */
/* > ( x ) ( 0 ) */
/* > */
/* > where alpha and beta are scalars, with beta real, and x is an */
/* > (n-1)-element complex vector. H is represented in the form */
/* > */
/* > H = I - tau * ( 1 ) * ( 1 v**H ) , */
/* > ( v ) */
/* > */
/* > where tau is a complex scalar and v is a complex (n-1)-element */
/* > vector. Note that H is not hermitian. */
/* > */
/* > If the elements of x are all zero and alpha is real, then tau = 0 */
/* > and H is taken to be the unit matrix. */
/* > */
/* > Otherwise 1 <= real(tau) <= 2 and f2c_abs(tau-1) <= 1 . */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the elementary reflector. */
/* > \endverbatim */
/* > */
/* > \param[in,out] ALPHA */
/* > \verbatim */
/* > ALPHA is COMPLEX */
/* > On entry, the value alpha. */
/* > On exit, it is overwritten with the value beta. */
/* > \endverbatim */
/* > */
/* > \param[in,out] X */
/* > \verbatim */
/* > X is COMPLEX array, dimension */
/* > (1+(N-2)*f2c_abs(INCX)) */
/* > On entry, the vector x. */
/* > On exit, it is overwritten with the vector v. */
/* > \endverbatim */
/* > */
/* > \param[in] INCX */
/* > \verbatim */
/* > INCX is INTEGER */
/* > The increment between elements of X. INCX > 0. */
/* > \endverbatim */
/* > */
/* > \param[out] TAU */
/* > \verbatim */
/* > TAU is COMPLEX */
/* > The value tau. */
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
int clarfg_(integer *n, complex *alpha, complex *x, integer * incx, complex *tau)
{
    /* System generated locals */
    integer i__1;
    real r__1, r__2;
    complex q__1, q__2;
    /* Builtin functions */
    double r_imag(complex *), r_sign(real *, real *);
    /* Local variables */
    integer j, knt;
    real beta;
    extern /* Subroutine */
    int cscal_(integer *, complex *, complex *, integer *);
    real alphi, alphr, xnorm;
    extern real scnrm2_(integer *, complex *, integer *), slapy3_(real *, real *, real *);
    extern /* Complex */
    VOID cladiv_(complex *, complex *, complex *);
    extern real slamch_(char *);
    extern /* Subroutine */
    int csscal_(integer *, real *, complex *, integer *);
    real safmin, rsafmn;
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
    /* .. External Subroutines .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Parameter adjustments */
    --x;
    /* Function Body */
    if (*n <= 0)
    {
        tau->r = 0.f, tau->i = 0.f;
        return 0;
    }
    i__1 = *n - 1;
    xnorm = scnrm2_(&i__1, &x[1], incx);
    alphr = alpha->r;
    alphi = r_imag(alpha);
    if (xnorm == 0.f && alphi == 0.f)
    {
        /* H = I */
        tau->r = 0.f, tau->i = 0.f;
    }
    else
    {
        /* general case */
        r__1 = slapy3_(&alphr, &alphi, &xnorm);
        beta = -r_sign(&r__1, &alphr);
        safmin = slamch_("S") / slamch_("E");
        rsafmn = 1.f / safmin;
        knt = 0;
        if (f2c_abs(beta) < safmin)
        {
            /* XNORM, BETA may be inaccurate;
            scale X and recompute them */
L10:
            ++knt;
            i__1 = *n - 1;
            csscal_(&i__1, &rsafmn, &x[1], incx);
            beta *= rsafmn;
            alphi *= rsafmn;
            alphr *= rsafmn;
            if (f2c_abs(beta) < safmin)
            {
                goto L10;
            }
            /* New BETA is at most 1, at least SAFMIN */
            i__1 = *n - 1;
            xnorm = scnrm2_(&i__1, &x[1], incx);
            q__1.r = alphr;
            q__1.i = alphi; // , expr subst
            alpha->r = q__1.r, alpha->i = q__1.i;
            r__1 = slapy3_(&alphr, &alphi, &xnorm);
            beta = -r_sign(&r__1, &alphr);
        }
        r__1 = (beta - alphr) / beta;
        r__2 = -alphi / beta;
        q__1.r = r__1;
        q__1.i = r__2; // , expr subst
        tau->r = q__1.r, tau->i = q__1.i;
        q__2.r = alpha->r - beta;
        q__2.i = alpha->i; // , expr subst
        cladiv_(&q__1, &c_b5, &q__2);
        alpha->r = q__1.r, alpha->i = q__1.i;
        i__1 = *n - 1;
        cscal_(&i__1, alpha, &x[1], incx);
        /* If ALPHA is subnormal, it may lose relative accuracy */
        i__1 = knt;
        for (j = 1;
                j <= i__1;
                ++j)
        {
            beta *= safmin;
            /* L20: */
        }
        alpha->r = beta, alpha->i = 0.f;
    }
    return 0;
    /* End of CLARFG */
}
/* clarfg_ */
