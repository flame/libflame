/* ../netlib/zlartg.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief \b ZLARTG generates a plane rotation with real cosine and complex sine. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download ZLARTG + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlartg. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlartg. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlartg. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE ZLARTG( F, G, CS, SN, R ) */
/* .. Scalar Arguments .. */
/* DOUBLE PRECISION CS */
/* COMPLEX*16 F, G, R, SN */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZLARTG generates a plane rotation so that */
/* > */
/* > [ CS SN ] [ F ] [ R ] */
/* > [ __ ] . [ ] = [ ] where CS**2 + |SN|**2 = 1. */
/* > [ -SN CS ] [ G ] [ 0 ] */
/* > */
/* > This is a faster version of the BLAS1 routine ZROTG, except for */
/* > the following differences: */
/* > F and G are unchanged on return. */
/* > If G=0, then CS=1 and SN=0. */
/* > If F=0, then CS=0 and SN is chosen so that R is real. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] F */
/* > \verbatim */
/* > F is COMPLEX*16 */
/* > The first component of vector to be rotated. */
/* > \endverbatim */
/* > */
/* > \param[in] G */
/* > \verbatim */
/* > G is COMPLEX*16 */
/* > The second component of vector to be rotated. */
/* > \endverbatim */
/* > */
/* > \param[out] CS */
/* > \verbatim */
/* > CS is DOUBLE PRECISION */
/* > The cosine of the rotation. */
/* > \endverbatim */
/* > */
/* > \param[out] SN */
/* > \verbatim */
/* > SN is COMPLEX*16 */
/* > The sine of the rotation. */
/* > \endverbatim */
/* > */
/* > \param[out] R */
/* > \verbatim */
/* > R is COMPLEX*16 */
/* > The nonzero component of the rotated vector. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date November 2013 */
/* > \ingroup complex16OTHERauxiliary */
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
int zlartg_(doublecomplex *f, doublecomplex *g, doublereal * cs, doublecomplex *sn, doublecomplex *r__)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2, d__3, d__4, d__5, d__6, d__7, d__8, d__9, d__10;
    doublecomplex z__1, z__2, z__3;
    /* Builtin functions */
    double log(doublereal), pow_di(doublereal *, integer *), d_imag( doublecomplex *), z_abs(doublecomplex *), sqrt(doublereal);
    void d_cnjg(doublecomplex *, doublecomplex *);
    /* Local variables */
    doublereal d__;
    integer i__;
    doublereal f2, g2;
    doublecomplex ff;
    doublereal di, dr;
    doublecomplex fs, gs;
    doublereal f2s, g2s, eps, scale;
    integer count;
    doublereal safmn2;
    extern doublereal dlapy2_(doublereal *, doublereal *);
    doublereal safmx2;
    extern doublereal dlamch_(char *);
    extern logical disnan_(doublereal *);
    doublereal safmin;
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
    safmin = dlamch_("S");
    eps = dlamch_("E");
    d__1 = dlamch_("B");
    i__1 = (integer) (log(safmin / eps) / log(dlamch_("B")) / 2.);
    safmn2 = pow_di(&d__1, &i__1);
    safmx2 = 1. / safmn2;
    /* Computing MAX */
    /* Computing MAX */
    d__7 = (d__1 = f->r, f2c_abs(d__1));
    d__8 = (d__2 = d_imag(f), f2c_abs(d__2)); // , expr subst
    /* Computing MAX */
    d__9 = (d__3 = g->r, f2c_abs(d__3));
    d__10 = (d__4 = d_imag(g), f2c_abs(d__4)); // , expr subst
    d__5 = max(d__7,d__8);
    d__6 = max(d__9,d__10); // , expr subst
    scale = max(d__5,d__6);
    fs.r = f->r;
    fs.i = f->i; // , expr subst
    gs.r = g->r;
    gs.i = g->i; // , expr subst
    count = 0;
    if (scale >= safmx2)
    {
L10:
        ++count;
        z__1.r = safmn2 * fs.r;
        z__1.i = safmn2 * fs.i; // , expr subst
        fs.r = z__1.r;
        fs.i = z__1.i; // , expr subst
        z__1.r = safmn2 * gs.r;
        z__1.i = safmn2 * gs.i; // , expr subst
        gs.r = z__1.r;
        gs.i = z__1.i; // , expr subst
        scale *= safmn2;
        if (scale >= safmx2)
        {
            goto L10;
        }
    }
    else if (scale <= safmn2)
    {
        d__1 = z_abs(g);
        if (g->r == 0. && g->i == 0. || disnan_(&d__1))
        {
            *cs = 1.;
            sn->r = 0., sn->i = 0.;
            r__->r = f->r, r__->i = f->i;
            return 0;
        }
L20:
        --count;
        z__1.r = safmx2 * fs.r;
        z__1.i = safmx2 * fs.i; // , expr subst
        fs.r = z__1.r;
        fs.i = z__1.i; // , expr subst
        z__1.r = safmx2 * gs.r;
        z__1.i = safmx2 * gs.i; // , expr subst
        gs.r = z__1.r;
        gs.i = z__1.i; // , expr subst
        scale *= safmx2;
        if (scale <= safmn2)
        {
            goto L20;
        }
    }
    /* Computing 2nd power */
    d__1 = fs.r;
    /* Computing 2nd power */
    d__2 = d_imag(&fs);
    f2 = d__1 * d__1 + d__2 * d__2;
    /* Computing 2nd power */
    d__1 = gs.r;
    /* Computing 2nd power */
    d__2 = d_imag(&gs);
    g2 = d__1 * d__1 + d__2 * d__2;
    if (f2 <= max(g2,1.) * safmin)
    {
        /* This is a rare case: F is very small. */
        if (f->r == 0. && f->i == 0.)
        {
            *cs = 0.;
            d__2 = g->r;
            d__3 = d_imag(g);
            d__1 = dlapy2_(&d__2, &d__3);
            r__->r = d__1, r__->i = 0.;
            /* Do complex/real division explicitly with two real divisions */
            d__1 = gs.r;
            d__2 = d_imag(&gs);
            d__ = dlapy2_(&d__1, &d__2);
            d__1 = gs.r / d__;
            d__2 = -d_imag(&gs) / d__;
            z__1.r = d__1;
            z__1.i = d__2; // , expr subst
            sn->r = z__1.r, sn->i = z__1.i;
            return 0;
        }
        d__1 = fs.r;
        d__2 = d_imag(&fs);
        f2s = dlapy2_(&d__1, &d__2);
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
        d__3 = (d__1 = f->r, f2c_abs(d__1));
        d__4 = (d__2 = d_imag(f), f2c_abs(d__2)); // , expr subst
        if (max(d__3,d__4) > 1.)
        {
            d__1 = f->r;
            d__2 = d_imag(f);
            d__ = dlapy2_(&d__1, &d__2);
            d__1 = f->r / d__;
            d__2 = d_imag(f) / d__;
            z__1.r = d__1;
            z__1.i = d__2; // , expr subst
            ff.r = z__1.r;
            ff.i = z__1.i; // , expr subst
        }
        else
        {
            dr = safmx2 * f->r;
            di = safmx2 * d_imag(f);
            d__ = dlapy2_(&dr, &di);
            d__1 = dr / d__;
            d__2 = di / d__;
            z__1.r = d__1;
            z__1.i = d__2; // , expr subst
            ff.r = z__1.r;
            ff.i = z__1.i; // , expr subst
        }
        d__1 = gs.r / g2s;
        d__2 = -d_imag(&gs) / g2s;
        z__2.r = d__1;
        z__2.i = d__2; // , expr subst
        z__1.r = ff.r * z__2.r - ff.i * z__2.i;
        z__1.i = ff.r * z__2.i + ff.i * z__2.r; // , expr subst
        sn->r = z__1.r, sn->i = z__1.i;
        z__2.r = *cs * f->r;
        z__2.i = *cs * f->i; // , expr subst
        z__3.r = sn->r * g->r - sn->i * g->i;
        z__3.i = sn->r * g->i + sn->i * g->r; // , expr subst
        z__1.r = z__2.r + z__3.r;
        z__1.i = z__2.i + z__3.i; // , expr subst
        r__->r = z__1.r, r__->i = z__1.i;
    }
    else
    {
        /* This is the most common case. */
        /* Neither F2 nor F2/G2 are less than SAFMIN */
        /* F2S cannot overflow, and it is accurate */
        f2s = sqrt(g2 / f2 + 1.);
        /* Do the F2S(real)*FS(complex) multiply with two real multiplies */
        d__1 = f2s * fs.r;
        d__2 = f2s * d_imag(&fs);
        z__1.r = d__1;
        z__1.i = d__2; // , expr subst
        r__->r = z__1.r, r__->i = z__1.i;
        *cs = 1. / f2s;
        d__ = f2 + g2;
        /* Do complex/real division explicitly with two real divisions */
        d__1 = r__->r / d__;
        d__2 = d_imag(r__) / d__;
        z__1.r = d__1;
        z__1.i = d__2; // , expr subst
        sn->r = z__1.r, sn->i = z__1.i;
        d_cnjg(&z__2, &gs);
        z__1.r = sn->r * z__2.r - sn->i * z__2.i;
        z__1.i = sn->r * z__2.i + sn->i * z__2.r; // , expr subst
        sn->r = z__1.r, sn->i = z__1.i;
        if (count != 0)
        {
            if (count > 0)
            {
                i__1 = count;
                for (i__ = 1;
                        i__ <= i__1;
                        ++i__)
                {
                    z__1.r = safmx2 * r__->r;
                    z__1.i = safmx2 * r__->i; // , expr subst
                    r__->r = z__1.r, r__->i = z__1.i;
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
                    z__1.r = safmn2 * r__->r;
                    z__1.i = safmn2 * r__->i; // , expr subst
                    r__->r = z__1.r, r__->i = z__1.i;
                    /* L40: */
                }
            }
        }
    }
    return 0;
    /* End of ZLARTG */
}
/* zlartg_ */
