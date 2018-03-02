/* ../netlib/clartg.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief \b CLARTG generates a plane rotation with real cosine and complex sine. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CLARTG + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clartg. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clartg. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clartg. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CLARTG( F, G, CS, SN, R ) */
/* .. Scalar Arguments .. */
/* REAL CS */
/* COMPLEX F, G, R, SN */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CLARTG generates a plane rotation so that */
/* > */
/* > [ CS SN ] [ F ] [ R ] */
/* > [ __ ] . [ ] = [ ] where CS**2 + |SN|**2 = 1. */
/* > [ -SN CS ] [ G ] [ 0 ] */
/* > */
/* > This is a faster version of the BLAS1 routine CROTG, except for */
/* > the following differences: */
/* > F and G are unchanged on return. */
/* > If G=0, then CS=1 and SN=0. */
/* > If F=0, then CS=0 and SN is chosen so that R is real. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] F */
/* > \verbatim */
/* > F is COMPLEX */
/* > The first component of vector to be rotated. */
/* > \endverbatim */
/* > */
/* > \param[in] G */
/* > \verbatim */
/* > G is COMPLEX */
/* > The second component of vector to be rotated. */
/* > \endverbatim */
/* > */
/* > \param[out] CS */
/* > \verbatim */
/* > CS is REAL */
/* > The cosine of the rotation. */
/* > \endverbatim */
/* > */
/* > \param[out] SN */
/* > \verbatim */
/* > SN is COMPLEX */
/* > The sine of the rotation. */
/* > \endverbatim */
/* > */
/* > \param[out] R */
/* > \verbatim */
/* > R is COMPLEX */
/* > The nonzero component of the rotated vector. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date November 2013 */
/* > \ingroup complexOTHERauxiliary */
/* > \par Further Details: */
/* ===================== */
/* > */
/* > \verbatim */
/* > */
/* > 3-5-96 - Modified with a new algorithm by W. Kahan and J. Demmel */
/* > */
/* > This version has a few statements commented out for thread safety */
/* > (machine parameters are computed on each entry). 10 feb 03, SJH. */
/* > \endverbatim */
/* > */
/* ===================================================================== */
/* Subroutine */
int clartg_(complex *f, complex *g, real *cs, complex *sn, complex *r__)
{
    /* System generated locals */
    integer i__1;
    real r__1, r__2, r__3, r__4, r__5, r__6, r__7, r__8, r__9, r__10;
    complex q__1, q__2, q__3;
    /* Builtin functions */
    double log(doublereal), pow_ri(real *, integer *), r_imag(complex *), c_abs(complex *), sqrt(doublereal);
    void r_cnjg(complex *, complex *);
    /* Local variables */
    real d__;
    integer i__;
    real f2, g2;
    complex ff;
    real di, dr;
    complex fs, gs;
    real f2s, g2s, eps, scale;
    integer count;
    real safmn2, safmx2;
    extern real slapy2_(real *, real *), slamch_(char *);
    real safmin;
    extern logical sisnan_(real *);
    /* -- LAPACK auxiliary routine (version 3.5.0) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* November 2013 */
    /* .. Scalar Arguments .. */
    /* .. */
    /* ===================================================================== */
    /* .. Parameters .. */
    /* .. */
    /* .. Local Scalars .. */
    /* LOGICAL FIRST */
    /* .. */
    /* .. External Functions .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Statement Functions .. */
    /* .. */
    /* .. Statement Function definitions .. */
    /* .. */
    /* .. Executable Statements .. */
    safmin = slamch_("S");
    eps = slamch_("E");
    r__1 = slamch_("B");
    i__1 = (integer) (log(safmin / eps) / log(slamch_("B")) / 2.f);
    safmn2 = pow_ri(&r__1, &i__1);
    safmx2 = 1.f / safmn2;
    /* Computing MAX */
    /* Computing MAX */
    r__7 = (r__1 = f->r, f2c_abs(r__1));
    r__8 = (r__2 = r_imag(f), f2c_abs(r__2)); // , expr subst
    /* Computing MAX */
    r__9 = (r__3 = g->r, f2c_abs(r__3));
    r__10 = (r__4 = r_imag(g), f2c_abs(r__4)); // , expr subst
    r__5 = max(r__7,r__8);
    r__6 = max(r__9,r__10); // , expr subst
    scale = max(r__5,r__6);
    fs.r = f->r;
    fs.i = f->i; // , expr subst
    gs.r = g->r;
    gs.i = g->i; // , expr subst
    count = 0;
    if (scale >= safmx2)
    {
L10:
        ++count;
        q__1.r = safmn2 * fs.r;
        q__1.i = safmn2 * fs.i; // , expr subst
        fs.r = q__1.r;
        fs.i = q__1.i; // , expr subst
        q__1.r = safmn2 * gs.r;
        q__1.i = safmn2 * gs.i; // , expr subst
        gs.r = q__1.r;
        gs.i = q__1.i; // , expr subst
        scale *= safmn2;
        if (scale >= safmx2)
        {
            goto L10;
        }
    }
    else if (scale <= safmn2)
    {
        r__1 = c_abs(g);
        if (g->r == 0.f && g->i == 0.f || sisnan_(&r__1))
        {
            *cs = 1.f;
            sn->r = 0.f, sn->i = 0.f;
            r__->r = f->r, r__->i = f->i;
            return 0;
        }
L20:
        --count;
        q__1.r = safmx2 * fs.r;
        q__1.i = safmx2 * fs.i; // , expr subst
        fs.r = q__1.r;
        fs.i = q__1.i; // , expr subst
        q__1.r = safmx2 * gs.r;
        q__1.i = safmx2 * gs.i; // , expr subst
        gs.r = q__1.r;
        gs.i = q__1.i; // , expr subst
        scale *= safmx2;
        if (scale <= safmn2)
        {
            goto L20;
        }
    }
    /* Computing 2nd power */
    r__1 = fs.r;
    /* Computing 2nd power */
    r__2 = r_imag(&fs);
    f2 = r__1 * r__1 + r__2 * r__2;
    /* Computing 2nd power */
    r__1 = gs.r;
    /* Computing 2nd power */
    r__2 = r_imag(&gs);
    g2 = r__1 * r__1 + r__2 * r__2;
    if (f2 <= max(g2,1.f) * safmin)
    {
        /* This is a rare case: F is very small. */
        if (f->r == 0.f && f->i == 0.f)
        {
            *cs = 0.f;
            r__2 = g->r;
            r__3 = r_imag(g);
            r__1 = slapy2_(&r__2, &r__3);
            r__->r = r__1, r__->i = 0.f;
            /* Do complex/real division explicitly with two real divisions */
            r__1 = gs.r;
            r__2 = r_imag(&gs);
            d__ = slapy2_(&r__1, &r__2);
            r__1 = gs.r / d__;
            r__2 = -r_imag(&gs) / d__;
            q__1.r = r__1;
            q__1.i = r__2; // , expr subst
            sn->r = q__1.r, sn->i = q__1.i;
            return 0;
        }
        r__1 = fs.r;
        r__2 = r_imag(&fs);
        f2s = slapy2_(&r__1, &r__2);
        /* G2 and G2S are accurate */
        /* G2 is at least SAFMIN, and G2S is at least SAFMN2 */
        g2s = sqrt(g2);
        /* Error in CS from underflow in F2S is at most */
        /* UNFL / SAFMN2 .lt. sqrt(UNFL*EPS) .lt. EPS */
        /* If MAX(G2,ONE)=G2, then F2 .lt. G2*SAFMIN, */
        /* and so CS .lt. sqrt(SAFMIN) */
        /* If MAX(G2,ONE)=ONE, then F2 .lt. SAFMIN */
        /* and so CS .lt. sqrt(SAFMIN)/SAFMN2 = sqrt(EPS) */
        /* Therefore, CS = F2S/G2S / sqrt( 1 + (F2S/G2S)**2 ) = F2S/G2S */
        *cs = f2s / g2s;
        /* Make sure f2c_abs(FF) = 1 */
        /* Do complex/real division explicitly with 2 real divisions */
        /* Computing MAX */
        r__3 = (r__1 = f->r, f2c_abs(r__1));
        r__4 = (r__2 = r_imag(f), f2c_abs(r__2)); // , expr subst
        if (max(r__3,r__4) > 1.f)
        {
            r__1 = f->r;
            r__2 = r_imag(f);
            d__ = slapy2_(&r__1, &r__2);
            r__1 = f->r / d__;
            r__2 = r_imag(f) / d__;
            q__1.r = r__1;
            q__1.i = r__2; // , expr subst
            ff.r = q__1.r;
            ff.i = q__1.i; // , expr subst
        }
        else
        {
            dr = safmx2 * f->r;
            di = safmx2 * r_imag(f);
            d__ = slapy2_(&dr, &di);
            r__1 = dr / d__;
            r__2 = di / d__;
            q__1.r = r__1;
            q__1.i = r__2; // , expr subst
            ff.r = q__1.r;
            ff.i = q__1.i; // , expr subst
        }
        r__1 = gs.r / g2s;
        r__2 = -r_imag(&gs) / g2s;
        q__2.r = r__1;
        q__2.i = r__2; // , expr subst
        q__1.r = ff.r * q__2.r - ff.i * q__2.i;
        q__1.i = ff.r * q__2.i + ff.i * q__2.r; // , expr subst
        sn->r = q__1.r, sn->i = q__1.i;
        q__2.r = *cs * f->r;
        q__2.i = *cs * f->i; // , expr subst
        q__3.r = sn->r * g->r - sn->i * g->i;
        q__3.i = sn->r * g->i + sn->i * g->r; // , expr subst
        q__1.r = q__2.r + q__3.r;
        q__1.i = q__2.i + q__3.i; // , expr subst
        r__->r = q__1.r, r__->i = q__1.i;
    }
    else
    {
        /* This is the most common case. */
        /* Neither F2 nor F2/G2 are less than SAFMIN */
        /* F2S cannot overflow, and it is accurate */
        f2s = sqrt(g2 / f2 + 1.f);
        /* Do the F2S(real)*FS(complex) multiply with two real multiplies */
        r__1 = f2s * fs.r;
        r__2 = f2s * r_imag(&fs);
        q__1.r = r__1;
        q__1.i = r__2; // , expr subst
        r__->r = q__1.r, r__->i = q__1.i;
        *cs = 1.f / f2s;
        d__ = f2 + g2;
        /* Do complex/real division explicitly with two real divisions */
        r__1 = r__->r / d__;
        r__2 = r_imag(r__) / d__;
        q__1.r = r__1;
        q__1.i = r__2; // , expr subst
        sn->r = q__1.r, sn->i = q__1.i;
        r_cnjg(&q__2, &gs);
        q__1.r = sn->r * q__2.r - sn->i * q__2.i;
        q__1.i = sn->r * q__2.i + sn->i * q__2.r; // , expr subst
        sn->r = q__1.r, sn->i = q__1.i;
        if (count != 0)
        {
            if (count > 0)
            {
                i__1 = count;
                for (i__ = 1;
                        i__ <= i__1;
                        ++i__)
                {
                    q__1.r = safmx2 * r__->r;
                    q__1.i = safmx2 * r__->i; // , expr subst
                    r__->r = q__1.r, r__->i = q__1.i;
                    /* L30: */
                }
            }
            else
            {
                i__1 = -count;
                for (i__ = 1;
                        i__ <= i__1;
                        ++i__)
                {
                    q__1.r = safmn2 * r__->r;
                    q__1.i = safmn2 * r__->i; // , expr subst
                    r__->r = q__1.r, r__->i = q__1.i;
                    /* L40: */
                }
            }
        }
    }
    return 0;
    /* End of CLARTG */
}
/* clartg_ */
