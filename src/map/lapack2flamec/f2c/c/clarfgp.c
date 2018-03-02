/* ../netlib/clarfgp.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static complex c_b5 =
{
    1.f,0.f
}
;
/* > \brief \b CLARFGP generates an elementary reflector (Householder matrix) with non-negatibe beta. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CLARFGP + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clarfgp .f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clarfgp .f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clarfgp .f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CLARFGP( N, ALPHA, X, INCX, TAU ) */
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
/* > CLARFGP generates a complex elementary reflector H of order n, such */
/* > that */
/* > */
/* > H**H * ( alpha ) = ( beta ), H**H * H = I. */
/* > ( x ) ( 0 ) */
/* > */
/* > where alpha and beta are scalars, beta is real and non-negative, and */
/* > x is an (n-1)-element complex vector. H is represented in the form */
/* > */
/* > H = I - tau * ( 1 ) * ( 1 v**H ) , */
/* > ( v ) */
/* > */
/* > where tau is a complex scalar and v is a complex (n-1)-element */
/* > vector. Note that H is not hermitian. */
/* > */
/* > If the elements of x are all zero and alpha is real, then tau = 0 */
/* > and H is taken to be the unit matrix. */
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
int clarfgp_(integer *n, complex *alpha, complex *x, integer *incx, complex *tau)
{
    /* System generated locals */
    integer i__1, i__2;
    real r__1, r__2;
    complex q__1, q__2;
    /* Builtin functions */
    double r_imag(complex *), r_sign(real *, real *), c_abs(complex *);
    /* Local variables */
    integer j;
    complex savealpha;
    integer knt;
    real beta;
    extern /* Subroutine */
    int cscal_(integer *, complex *, complex *, integer *);
    real alphi, alphr, xnorm;
    extern real scnrm2_(integer *, complex *, integer *), slapy2_(real *, real *), slapy3_(real *, real *, real *);
    extern /* Complex */
    VOID cladiv_(complex *, complex *, complex *);
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
    if (xnorm == 0.f)
    {
        /* H = [1-alpha/f2c_abs(alpha) 0;
        0 I], sign chosen so ALPHA >= 0. */
        if (alphi == 0.f)
        {
            if (alphr >= 0.f)
            {
                /* When TAU.eq.ZERO, the vector is special-cased to be */
                /* all zeros in the application routines. We do not need */
                /* to clear it. */
                tau->r = 0.f, tau->i = 0.f;
            }
            else
            {
                /* However, the application routines rely on explicit */
                /* zero checks when TAU.ne.ZERO, and we must clear X. */
                tau->r = 2.f, tau->i = 0.f;
                i__1 = *n - 1;
                for (j = 1;
                        j <= i__1;
                        ++j)
                {
                    i__2 = (j - 1) * *incx + 1;
                    x[i__2].r = 0.f;
                    x[i__2].i = 0.f; // , expr subst
                }
                q__1.r = -alpha->r;
                q__1.i = -alpha->i; // , expr subst
                alpha->r = q__1.r, alpha->i = q__1.i;
            }
        }
        else
        {
            /* Only "reflecting" the diagonal entry to be real and non-negative. */
            xnorm = slapy2_(&alphr, &alphi);
            r__1 = 1.f - alphr / xnorm;
            r__2 = -alphi / xnorm;
            q__1.r = r__1;
            q__1.i = r__2; // , expr subst
            tau->r = q__1.r, tau->i = q__1.i;
            i__1 = *n - 1;
            for (j = 1;
                    j <= i__1;
                    ++j)
            {
                i__2 = (j - 1) * *incx + 1;
                x[i__2].r = 0.f;
                x[i__2].i = 0.f; // , expr subst
            }
            alpha->r = xnorm, alpha->i = 0.f;
        }
    }
    else
    {
        /* general case */
        r__1 = slapy3_(&alphr, &alphi, &xnorm);
        beta = r_sign(&r__1, &alphr);
        smlnum = slamch_("S") / slamch_("E");
        bignum = 1.f / smlnum;
        knt = 0;
        if (f2c_abs(beta) < smlnum)
        {
            /* XNORM, BETA may be inaccurate;
            scale X and recompute them */
L10:
            ++knt;
            i__1 = *n - 1;
            csscal_(&i__1, &bignum, &x[1], incx);
            beta *= bignum;
            alphi *= bignum;
            alphr *= bignum;
            if (f2c_abs(beta) < smlnum)
            {
                goto L10;
            }
            /* New BETA is at most 1, at least SMLNUM */
            i__1 = *n - 1;
            xnorm = scnrm2_(&i__1, &x[1], incx);
            q__1.r = alphr;
            q__1.i = alphi; // , expr subst
            alpha->r = q__1.r, alpha->i = q__1.i;
            r__1 = slapy3_(&alphr, &alphi, &xnorm);
            beta = r_sign(&r__1, &alphr);
        }
        savealpha.r = alpha->r;
        savealpha.i = alpha->i; // , expr subst
        q__1.r = alpha->r + beta;
        q__1.i = alpha->i; // , expr subst
        alpha->r = q__1.r, alpha->i = q__1.i;
        if (beta < 0.f)
        {
            beta = -beta;
            q__2.r = -alpha->r;
            q__2.i = -alpha->i; // , expr subst
            q__1.r = q__2.r / beta;
            q__1.i = q__2.i / beta; // , expr subst
            tau->r = q__1.r, tau->i = q__1.i;
        }
        else
        {
            alphr = alphi * (alphi / alpha->r);
            alphr += xnorm * (xnorm / alpha->r);
            r__1 = alphr / beta;
            r__2 = -alphi / beta;
            q__1.r = r__1;
            q__1.i = r__2; // , expr subst
            tau->r = q__1.r, tau->i = q__1.i;
            r__1 = -alphr;
            q__1.r = r__1;
            q__1.i = alphi; // , expr subst
            alpha->r = q__1.r, alpha->i = q__1.i;
        }
        cladiv_(&q__1, &c_b5, alpha);
        alpha->r = q__1.r, alpha->i = q__1.i;
        if (c_abs(tau) <= smlnum)
        {
            /* In the case where the computed TAU ends up being a denormalized number, */
            /* it loses relative accuracy. This is a BIG problem. Solution: flush TAU */
            /* to ZERO (or TWO or whatever makes a nonnegative real number for BETA). */
            /* (Bug report provided by Pat Quillen from MathWorks on Jul 29, 2009.) */
            /* (Thanks Pat. Thanks MathWorks.) */
            alphr = savealpha.r;
            alphi = r_imag(&savealpha);
            if (alphi == 0.f)
            {
                if (alphr >= 0.f)
                {
                    tau->r = 0.f, tau->i = 0.f;
                }
                else
                {
                    tau->r = 2.f, tau->i = 0.f;
                    i__1 = *n - 1;
                    for (j = 1;
                            j <= i__1;
                            ++j)
                    {
                        i__2 = (j - 1) * *incx + 1;
                        x[i__2].r = 0.f;
                        x[i__2].i = 0.f; // , expr subst
                    }
                    q__1.r = -savealpha.r;
                    q__1.i = -savealpha.i; // , expr subst
                    beta = q__1.r;
                }
            }
            else
            {
                xnorm = slapy2_(&alphr, &alphi);
                r__1 = 1.f - alphr / xnorm;
                r__2 = -alphi / xnorm;
                q__1.r = r__1;
                q__1.i = r__2; // , expr subst
                tau->r = q__1.r, tau->i = q__1.i;
                i__1 = *n - 1;
                for (j = 1;
                        j <= i__1;
                        ++j)
                {
                    i__2 = (j - 1) * *incx + 1;
                    x[i__2].r = 0.f;
                    x[i__2].i = 0.f; // , expr subst
                }
                beta = xnorm;
            }
        }
        else
        {
            /* This is the general case. */
            i__1 = *n - 1;
            cscal_(&i__1, alpha, &x[1], incx);
        }
        /* If BETA is subnormal, it may lose relative accuracy */
        i__1 = knt;
        for (j = 1;
                j <= i__1;
                ++j)
        {
            beta *= smlnum;
            /* L20: */
        }
        alpha->r = beta, alpha->i = 0.f;
    }
    return 0;
    /* End of CLARFGP */
}
/* clarfgp_ */
