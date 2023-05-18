/* zlartg.f -- translated by f2c (version 20190311). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Subroutine */
/* > \brief \b ZLARTG generates a plane rotation with real cosine and complex sine. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZLARTG( F, G, C, S, R ) */

/*       .. Scalar Arguments .. */
/*       REAL(wp)              C */
/*       COMPLEX(wp)           F, G, R, S */
/*       .. */

/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZLARTG generates a plane rotation so that */
/* > */
/* >    [  C         S  ] . [ F ]  =  [ R ] */
/* >    [ -conjg(S)  C  ]   [ G ]     [ 0 ] */
/* > */
/* > where C is real and C**2 + |S|**2 = 1. */
/* > */
/* > The mathematical formulas used for C and S are */
/* > */
/* >    sgn(x) = {  x / |x|,   x != 0 */
/* >             {  1,         x  = 0 */
/* > */
/* >    R = sgn(F) * sqrt(|F|**2 + |G|**2) */
/* > */
/* >    C = |F| / sqrt(|F|**2 + |G|**2) */
/* > */
/* >    S = sgn(F) * conjg(G) / sqrt(|F|**2 + |G|**2) */
/* > */
/* > Special conditions: */
/* >    If G=0, then C=1 and S=0. */
/* >    If F=0, then C=0 and S is chosen so that R is real. */
/* > */
/* > When F and G are real, the formulas simplify to C = F/R and */
/* > S = G/R, and the returned values of C, S, and R should be */
/* > identical to those returned by DLARTG. */
/* > */
/* > The algorithm used to compute these quantities incorporates scaling */
/* > to avoid overflow or underflow in computing the square root of the */
/* > sum of squares. */
/* > */
/* > This is the same routine ZROTG fom BLAS1, except that */
/* > F and G are unchanged on return. */
/* > */
/* > Below, wp=>dp stands for double precision from LA_CONSTANTS module. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] F */
/* > \verbatim */
/* >          F is COMPLEX(wp) */
/* >          The first component of vector to be rotated. */
/* > \endverbatim */
/* > */
/* > \param[in] G */
/* > \verbatim */
/* >          G is COMPLEX(wp) */
/* >          The second component of vector to be rotated. */
/* > \endverbatim */
/* > */
/* > \param[out] C */
/* > \verbatim */
/* >          C is REAL(wp) */
/* >          The cosine of the rotation. */
/* > \endverbatim */
/* > */
/* > \param[out] S */
/* > \verbatim */
/* >          S is COMPLEX(wp) */
/* >          The sine of the rotation. */
/* > \endverbatim */
/* > */
/* > \param[out] R */
/* > \verbatim */
/* >          R is COMPLEX(wp) */
/* >          The nonzero component of the rotated vector. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Weslley Pereira, University of Colorado Denver, USA */

/* > \date December 2021 */

/* > \ingroup OTHERauxiliary */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* > Based on the algorithm from */
/* > */
/* >  Anderson E. (2017) */
/* >  Algorithm 978: Safe Scaling in the Level 1 BLAS */
/* >  ACM Trans Math Softw 44:1--28 */
/* >  https://doi.org/10.1145/3061665 */
/* > */
/* > \endverbatim */
int zlartg_(doublecomplex *f, doublecomplex *g, doublereal * c__, doublecomplex *s, doublecomplex *r__)
{
    AOCL_DTL_TRACE_ENTRY_INDENT
    doublecomplex z__1, z__2, z__3;
    /* Builtin functions */
    double sqrt(doublereal), d_imag(doublecomplex *);
    void d_cnjg(doublecomplex *, doublecomplex *), z_div(doublecomplex *, doublecomplex *, doublecomplex *);
    /* Local variables */
    doublereal d__, u, v, w, f1, f2, g1, g2, h2;
    doublereal d__1, d__2, d__3, d__4;
    doublecomplex fs, gs;
    doublereal rtmin, rtmax, safmin, safmax;
    /* ...Translated by Pacific-Sierra Research vf90 Personal 3.4N3 04:17:29 1/20/23 */
    /*  -- LAPACK auxiliary routine -- */
    /*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
    /*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /*     February 2021 */
    /* .. Constants .. */
    safmin = 2.2250738585072014e-308;
    safmax = 4.4942328371557898e307;
    rtmin = sqrt(safmin);
    /* .. */
    /* .. Executable Statements .. */
    if (g->r == 0. && g->i == 0.)
    {
        *c__ = 1.;
        s->r = 0., s->i = 0.;
        r__->r = f->r, r__->i = f->i;
    }
    else if (f->r == 0. && f->i == 0.)
    {
        *c__ = 0.;
        if (g->r == 0.)
        {
            d__2 = (d__1 = d_imag(g), f2c_dabs(d__1));
            r__->r = d__2, r__->i = 0.;
            d_cnjg(&z__2, g);
            z_div(&z__1, &z__2, r__);
            s->r = z__1.r, s->i = z__1.i;
        }
        else if (d_imag(g) == 0.)
        {
            d__2 = (d__1 = g->r, f2c_dabs(d__1));
            r__->r = d__2, r__->i = 0.;
            d_cnjg(&z__2, g);
            z_div(&z__1, &z__2, r__);
            s->r = z__1.r, s->i = z__1.i;
        }
        else
        {
            /* Computing MAX */
            d__3 = (d__1 = g->r, f2c_dabs(d__1));
            d__4 = (d__2 = d_imag(g), f2c_dabs( d__2)); // , expr subst
            g1 = fla_max(d__3,d__4);
            rtmax = sqrt(safmax / 2);
            if (g1 > rtmin && g1 < rtmax)
            {
                /* Use unscaled algorithm */
                /* The following two lines can be replaced by `d = f2c_dabs( g )`. */
                /* This algorithm do not use the intrinsic complex abs. */
                /* Computing 2nd power */
                d__1 = g->r;
                /* Computing 2nd power */
                d__2 = d_imag(g);
                g2 = d__1 * d__1 + d__2 * d__2;
                d__ = sqrt(g2);
                d_cnjg(&z__2, g);
                z__1.r = z__2.r / d__;
                z__1.i = z__2.i / d__; // , expr subst
                s->r = z__1.r, s->i = z__1.i;
                r__->r = d__, r__->i = 0.;
            }
            else
            {
                /* Use scaled algorithm */
                /* Computing MIN */
                d__1 = safmax;
                d__2 = fla_max(safmin,g1); // , expr subst
                u = fla_min(d__1,d__2);
                z__1.r = g->r / u;
                z__1.i = g->i / u; // , expr subst
                gs.r = z__1.r;
                gs.i = z__1.i; // , expr subst
                /* The following two lines can be replaced by `d = f2c_dabs( gs )`. */
                /* This algorithm do not use the intrinsic complex abs. */
                /* Computing 2nd power */
                d__1 = gs.r;
                /* Computing 2nd power */
                d__2 = d_imag(&gs);
                g2 = d__1 * d__1 + d__2 * d__2;
                d__ = sqrt(g2);
                d_cnjg(&z__2, &gs);
                z__1.r = z__2.r / d__;
                z__1.i = z__2.i / d__; // , expr subst
                s->r = z__1.r, s->i = z__1.i;
                d__1 = d__ * u;
                r__->r = d__1, r__->i = 0.;
            }
        }
    }
    else
    {
        /* Computing MAX */
        d__3 = (d__1 = f->r, f2c_dabs(d__1));
        d__4 = (d__2 = d_imag(f), f2c_dabs(d__2)); // , expr subst
        f1 = fla_max(d__3,d__4);
        /* Computing MAX */
        d__3 = (d__1 = g->r, f2c_dabs(d__1));
        d__4 = (d__2 = d_imag(g), f2c_dabs(d__2)); // , expr subst
        g1 = fla_max(d__3,d__4);
        rtmax = sqrt(safmax / 4);
        if (f1 > rtmin && f1 < rtmax && g1 > rtmin && g1 < rtmax)
        {
            /* Use unscaled algorithm */
            /* Computing 2nd power */
            d__1 = f->r;
            /* Computing 2nd power */
            d__2 = d_imag(f);
            f2 = d__1 * d__1 + d__2 * d__2;
            /* Computing 2nd power */
            d__1 = g->r;
            /* Computing 2nd power */
            d__2 = d_imag(g);
            g2 = d__1 * d__1 + d__2 * d__2;
            h2 = f2 + g2;
            /* safmin <= f2 <= h2 <= safmax */
            if (f2 >= h2 * safmin)
            {
                /* safmin <= f2/h2 <= 1, and h2/f2 is finite */
                *c__ = sqrt(f2 / h2);
                z__1.r = f->r / *c__;
                z__1.i = f->i / *c__; // , expr subst
                r__->r = z__1.r, r__->i = z__1.i;
                rtmax *= 2;
                if (f2 > rtmin && h2 < rtmax)
                {
                    /* safmin <= sqrt( f2*h2 ) <= safmax */
                    d_cnjg(&z__2, g);
                    d__1 = sqrt(f2 * h2);
                    z__3.r = f->r / d__1;
                    z__3.i = f->i / d__1; // , expr subst
                    z__1.r = z__2.r * z__3.r - z__2.i * z__3.i;
                    z__1.i = z__2.r * z__3.i + z__2.i * z__3.r; // , expr subst
                    s->r = z__1.r, s->i = z__1.i;
                }
                else
                {
                    d_cnjg(&z__2, g);
                    z__3.r = r__->r / h2;
                    z__3.i = r__->i / h2; // , expr subst
                    z__1.r = z__2.r * z__3.r - z__2.i * z__3.i;
                    z__1.i = z__2.r * z__3.i + z__2.i * z__3.r; // , expr subst
                    s->r = z__1.r, s->i = z__1.i;
                }
            }
            else
            {
                /* f2/h2 <= safmin may be subnormal, and h2/f2 may overflow. */
                /* Moreover, */
                /* safmin <= f2*f2 * safmax < f2 * h2 < h2*h2 * safmin <= sa */
                /* sqrt(safmin) <= sqrt(f2 * h2) <= sqrt(safmax). */
                /* Also, */
                /* g2 >> f2, which means that h2 = g2. */
                d__ = sqrt(f2 * h2);
                *c__ = f2 / d__;
                if (*c__ >= safmin)
                {
                    z__1.r = f->r / *c__;
                    z__1.i = f->i / *c__; // , expr subst
                    r__->r = z__1.r, r__->i = z__1.i;
                }
                else
                {
                    /* f2 / sqrt(f2 * h2) < safmin, then */
                    /* sqrt(safmin) <= f2 * sqrt(safmax) <= h2 / sqrt(f2 * h2 */
                    d__1 = h2 / d__;
                    z__1.r = d__1 * f->r;
                    z__1.i = d__1 * f->i; // , expr subst
                    r__->r = z__1.r, r__->i = z__1.i;
                }
                d_cnjg(&z__2, g);
                z__3.r = f->r / d__;
                z__3.i = f->i / d__; // , expr subst
                z__1.r = z__2.r * z__3.r - z__2.i * z__3.i;
                z__1.i = z__2.r * z__3.i + z__2.i * z__3.r; // , expr subst
                s->r = z__1.r, s->i = z__1.i;
            }
        }
        else
        {
            /* Use scaled algorithm */
            /* Computing MIN */
            /* Computing MAX */
            d__3 = fla_max(safmin,f1);
            d__1 = safmax;
            d__2 = fla_max(d__3,g1); // , expr subst
            u = fla_min(d__1,d__2);
            z__1.r = g->r / u;
            z__1.i = g->i / u; // , expr subst
            gs.r = z__1.r;
            gs.i = z__1.i; // , expr subst
            /* Computing 2nd power */
            d__1 = gs.r;
            /* Computing 2nd power */
            d__2 = d_imag(&gs);
            g2 = d__1 * d__1 + d__2 * d__2;
            if (f1 / u < rtmin)
            {
                /* f is not well-scaled when scaled by g1. */
                /* Use a different scaling for f. */
                /* Computing MIN */
                d__1 = safmax;
                d__2 = fla_max(safmin,f1); // , expr subst
                v = fla_min(d__1,d__2);
                w = v / u;
                z__1.r = f->r / v;
                z__1.i = f->i / v; // , expr subst
                fs.r = z__1.r;
                fs.i = z__1.i; // , expr subst
                /* Computing 2nd power */
                d__1 = fs.r;
                /* Computing 2nd power */
                d__2 = d_imag(&fs);
                f2 = d__1 * d__1 + d__2 * d__2;
                /* Computing 2nd power */
                d__1 = w;
                h2 = f2 * (d__1 * d__1) + g2;
            }
            else
            {
                /* Otherwise use the same scaling for f and g. */
                w = 1.;
                z__1.r = f->r / u;
                z__1.i = f->i / u; // , expr subst
                fs.r = z__1.r;
                fs.i = z__1.i; // , expr subst
                /* Computing 2nd power */
                d__1 = fs.r;
                /* Computing 2nd power */
                d__2 = d_imag(&fs);
                f2 = d__1 * d__1 + d__2 * d__2;
                h2 = f2 + g2;
            }
            /* safmin <= f2 <= h2 <= safmax */
            if (f2 >= h2 * safmin)
            {
                /* safmin <= f2/h2 <= 1, and h2/f2 is finite */
                *c__ = sqrt(f2 / h2);
                z__1.r = fs.r / *c__;
                z__1.i = fs.i / *c__; // , expr subst
                r__->r = z__1.r, r__->i = z__1.i;
                rtmax *= 2;
                if (f2 > rtmin && h2 < rtmax)
                {
                    /* safmin <= sqrt( f2*h2 ) <= safmax */
                    d_cnjg(&z__2, &gs);
                    d__1 = sqrt(f2 * h2);
                    z__3.r = fs.r / d__1;
                    z__3.i = fs.i / d__1; // , expr subst
                    z__1.r = z__2.r * z__3.r - z__2.i * z__3.i;
                    z__1.i = z__2.r * z__3.i + z__2.i * z__3.r; // , expr subst
                    s->r = z__1.r, s->i = z__1.i;
                }
                else
                {
                    d_cnjg(&z__2, &gs);
                    z__3.r = r__->r / h2;
                    z__3.i = r__->i / h2; // , expr subst
                    z__1.r = z__2.r * z__3.r - z__2.i * z__3.i;
                    z__1.i = z__2.r * z__3.i + z__2.i * z__3.r; // , expr subst
                    s->r = z__1.r, s->i = z__1.i;
                }
            }
            else
            {
                /* f2/h2 <= safmin may be subnormal, and h2/f2 may overflow. */
                /* Moreover, */
                /* safmin <= f2*f2 * safmax < f2 * h2 < h2*h2 * safmin <= sa */
                /* sqrt(safmin) <= sqrt(f2 * h2) <= sqrt(safmax). */
                /* Also, */
                /* g2 >> f2, which means that h2 = g2. */
                d__ = sqrt(f2 * h2);
                *c__ = f2 / d__;
                if (*c__ >= safmin)
                {
                    z__1.r = fs.r / *c__;
                    z__1.i = fs.i / *c__; // , expr subst
                    r__->r = z__1.r, r__->i = z__1.i;
                }
                else
                {
                    /* f2 / sqrt(f2 * h2) < safmin, then */
                    /* sqrt(safmin) <= f2 * sqrt(safmax) <= h2 / sqrt(f2 * h2 */
                    d__1 = h2 / d__;
                    z__1.r = d__1 * fs.r;
                    z__1.i = d__1 * fs.i; // , expr subst
                    r__->r = z__1.r, r__->i = z__1.i;
                }
                d_cnjg(&z__2, &gs);
                z__3.r = fs.r / d__;
                z__3.i = fs.i / d__; // , expr subst
                z__1.r = z__2.r * z__3.r - z__2.i * z__3.i;
                z__1.i = z__2.r * z__3.i + z__2.i * z__3.r; // , expr subst
                s->r = z__1.r, s->i = z__1.i;
            }
            /* Rescale c and r */
            *c__ *= w;
            z__1.r = u * r__->r;
            z__1.i = u * r__->i; // , expr subst
            r__->r = z__1.r, r__->i = z__1.i;
        }
    }
    AOCL_DTL_TRACE_EXIT_INDENT
    return 0;
}
/* zlartg_ */
