/* clartg.f -- translated by f2c (version 20190311). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Subroutine */
/* > \brief \b CLARTG generates a plane rotation with real cosine and complex sine. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CLARTG( F, G, C, S, R ) */

/*       .. Scalar Arguments .. */
/*       REAL(wp)              C */
/*       COMPLEX(wp)           F, G, R, S */
/*       .. */

/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CLARTG generates a plane rotation so that */
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
/* > identical to those returned by SLARTG. */
/* > */
/* > The algorithm used to compute these quantities incorporates scaling */
/* > to avoid overflow or underflow in computing the square root of the */
/* > sum of squares. */
/* > */
/* > This is the same routine CROTG fom BLAS1, except that */
/* > F and G are unchanged on return. */
/* > */
/* > Below, wp=>sp stands for single precision from LA_CONSTANTS module. */
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
int clartg_(complex *f, complex *g, real *c__, complex *s, complex *r__)
{
    /* System generated locals */
    real r__1, r__2, r__3, r__4;
    complex q__1, q__2, q__3;
    /* Builtin functions */
    double sqrt(doublereal);
    void c_div(complex *, complex *, complex *);
    /* Local variables */
    real d__, u, v, w, f1, f2, g1, g2, h2;
    complex fs, gs;
    real rtmin, rtmax, safmin, safmax;
    /* ...Translated by Pacific-Sierra Research vf90 Personal 3.4N3 00:37:48 1/24/23 */
    /* ...Switches: */
    /*  -- LAPACK auxiliary routine -- */
    /*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
    /*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /*     February 2021 */    
    /* .. */
    /* .. Local Scalars .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Statement Functions .. */
    /* .. */
    /* .. Statement Function definitions .. */
    /* .. */
    /* .. Constants .. */
    safmin = 1.1754943508222875e-38f;
    safmax = 8.5070591730234616e37f;
    rtmin = sqrt(safmin);
    /* .. */
    /* .. Executable Statements .. */
    if (g->r == 0.f && g->i == 0.f)
    {
        *c__ = 1.f;
        s->r = 0.f, s->i = 0.f;
        r__->r = f->r, r__->i = f->i;
    }
    else if (f->r == 0.f && f->i == 0.f)
    {
        *c__ = 0.f;
        if (g->r == 0.f)
        {
            r__2 = (r__1 = g->i, f2c_abs(r__1));
            r__->r = r__2, r__->i = 0.f;
            q__2.r = g->r;
            q__2.i = -g->i;
            c_div(&q__1, &q__2, r__);
            s->r = q__1.r, s->i = q__1.i;
        }
        else if (g->i == 0.f)
        {
            r__2 = (r__1 = g->r, f2c_abs(r__1));
            r__->r = r__2, r__->i = 0.f;
            q__2.r = g->r;
            q__2.i = -g->i;
            c_div(&q__1, &q__2, r__);
            s->r = q__1.r, s->i = q__1.i;
        }
        else
        {
            /* Computing MAX */
            r__3 = (r__1 = g->r, f2c_abs(r__1));
            r__4 = (r__2 = g->i, f2c_abs( r__2)); // , expr subst
            g1 = fla_max(r__3,r__4);
            rtmax = sqrt(safmax / 2);
            if (g1 > rtmin && g1 < rtmax)
            {
                /* Use unscaled algorithm */
                /* The following two lines can be replaced by `d = f2c_abs( g )`. */
                /* This algorithm do not use the intrinsic complex abs. */
                /* Computing 2nd power */
                r__1 = g->r;
                /* Computing 2nd power */
                r__2 = g->i;
                g2 = r__1 * r__1 + r__2 * r__2;
                d__ = sqrt(g2);
                q__2.r = g->r;
                q__2.i = -g->i;
                q__1.r = q__2.r / d__;
                q__1.i = q__2.i / d__; // , expr subst
                s->r = q__1.r, s->i = q__1.i;
                r__->r = d__, r__->i = 0.f;
            }
            else
            {
                /* Use scaled algorithm */
                /* Computing MIN */
                r__1 = safmax;
                r__2 = fla_max(safmin,g1); // , expr subst
                u = fla_min(r__1,r__2);
                q__1.r = g->r / u;
                q__1.i = g->i / u; // , expr subst
                gs.r = q__1.r;
                gs.i = q__1.i; // , expr subst
                /* The following two lines can be replaced by `d = f2c_abs( gs )`. */
                /* This algorithm do not use the intrinsic complex abs. */
                /* Computing 2nd power */
                r__1 = gs.r;
                /* Computing 2nd power */
                r__2 = gs.i;
                g2 = r__1 * r__1 + r__2 * r__2;
                d__ = sqrt(g2);
                q__2.r = gs.r;
                q__2.i = -gs.i;
                q__1.r = q__2.r / d__;
                q__1.i = q__2.i / d__; // , expr subst
                s->r = q__1.r, s->i = q__1.i;
                r__1 = d__ * u;
                r__->r = r__1, r__->i = 0.f;
            }
        }
    }
    else
    {
        /* Computing MAX */
        r__3 = (r__1 = f->r, f2c_abs(r__1));
        r__4 = (r__2 = f->i, f2c_abs(r__2)); // , expr subst
        f1 = fla_max(r__3,r__4);
        /* Computing MAX */
        r__3 = (r__1 = g->r, f2c_abs(r__1));
        r__4 = (r__2 = g->i, f2c_abs(r__2)); // , expr subst
        g1 = fla_max(r__3,r__4);
        rtmax = sqrt(safmax / 4);
        if (f1 > rtmin && f1 < rtmax && g1 > rtmin && g1 < rtmax)
        {
            /* Use unscaled algorithm */
            /* Computing 2nd power */
            r__1 = f->r;
            /* Computing 2nd power */
            r__2 = f->i;
            f2 = r__1 * r__1 + r__2 * r__2;
            /* Computing 2nd power */
            r__1 = g->r;
            /* Computing 2nd power */
            r__2 = g->i;
            g2 = r__1 * r__1 + r__2 * r__2;
            h2 = f2 + g2;
            /* safmin <= f2 <= h2 <= safmax */
            if (f2 >= h2 * safmin)
            {
                /* safmin <= f2/h2 <= 1, and h2/f2 is finite */
                *c__ = sqrt(f2 / h2);
                q__1.r = f->r / *c__;
                q__1.i = f->i / *c__; // , expr subst
                r__->r = q__1.r, r__->i = q__1.i;
                rtmax *= 2;
                if (f2 > rtmin && h2 < rtmax)
                {
                    /* safmin <= sqrt( f2*h2 ) <= safmax */
                    q__2.r = g->r;
                    q__2.i = -g->i;
                    r__1 = sqrt(f2 * h2);
                    q__3.r = f->r / r__1;
                    q__3.i = f->i / r__1; // , expr subst
                    q__1.r = q__2.r * q__3.r - q__2.i * q__3.i;
                    q__1.i = q__2.r * q__3.i + q__2.i * q__3.r; // , expr subst
                    s->r = q__1.r, s->i = q__1.i;
                }
                else
                {
                    q__2.r = g->r;
                    q__2.i = -g->i;
                    q__3.r = r__->r / h2;
                    q__3.i = r__->i / h2; // , expr subst
                    q__1.r = q__2.r * q__3.r - q__2.i * q__3.i;
                    q__1.i = q__2.r * q__3.i + q__2.i * q__3.r; // , expr subst
                    s->r = q__1.r, s->i = q__1.i;
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
                    q__1.r = f->r / *c__;
                    q__1.i = f->i / *c__; // , expr subst
                    r__->r = q__1.r, r__->i = q__1.i;
                }
                else
                {
                    /* f2 / sqrt(f2 * h2) < safmin, then */
                    /* sqrt(safmin) <= f2 * sqrt(safmax) <= h2 / sqrt(f2 * h2 */
                    r__1 = h2 / d__;
                    q__1.r = r__1 * f->r;
                    q__1.i = r__1 * f->i; // , expr subst
                    r__->r = q__1.r, r__->i = q__1.i;
                }
                q__2.r = g->r;
                q__2.i = -g->i;
                q__3.r = f->r / d__;
                q__3.i = f->i / d__; // , expr subst
                q__1.r = q__2.r * q__3.r - q__2.i * q__3.i;
                q__1.i = q__2.r * q__3.i + q__2.i * q__3.r; // , expr subst
                s->r = q__1.r, s->i = q__1.i;
            }
        }
        else
        {
            /* Use scaled algorithm */
            /* Computing MIN */
            /* Computing MAX */
            r__3 = fla_max(safmin,f1);
            r__1 = safmax;
            r__2 = fla_max(r__3,g1); // , expr subst
            u = fla_min(r__1,r__2);
            q__1.r = g->r / u;
            q__1.i = g->i / u; // , expr subst
            gs.r = q__1.r;
            gs.i = q__1.i; // , expr subst
            /* Computing 2nd power */
            r__1 = gs.r;
            /* Computing 2nd power */
            r__2 = gs.i;
            g2 = r__1 * r__1 + r__2 * r__2;
            if (f1 / u < rtmin)
            {
                /* f is not well-scaled when scaled by g1. */
                /* Use a different scaling for f. */
                /* Computing MIN */
                r__1 = safmax;
                r__2 = fla_max(safmin,f1); // , expr subst
                v = fla_min(r__1,r__2);
                w = v / u;
                q__1.r = f->r / v;
                q__1.i = f->i / v; // , expr subst
                fs.r = q__1.r;
                fs.i = q__1.i; // , expr subst
                /* Computing 2nd power */
                r__1 = fs.r;
                /* Computing 2nd power */
                r__2 = fs.i;
                f2 = r__1 * r__1 + r__2 * r__2;
                /* Computing 2nd power */
                r__1 = w;
                h2 = f2 * (r__1 * r__1) + g2;
            }
            else
            {
                /* Otherwise use the same scaling for f and g. */
                w = 1.f;
                q__1.r = f->r / u;
                q__1.i = f->i / u; // , expr subst
                fs.r = q__1.r;
                fs.i = q__1.i; // , expr subst
                /* Computing 2nd power */
                r__1 = fs.r;
                /* Computing 2nd power */
                r__2 = fs.i;
                f2 = r__1 * r__1 + r__2 * r__2;
                h2 = f2 + g2;
            }
            /* safmin <= f2 <= h2 <= safmax */
            if (f2 >= h2 * safmin)
            {
                /* safmin <= f2/h2 <= 1, and h2/f2 is finite */
                *c__ = sqrt(f2 / h2);
                q__1.r = fs.r / *c__;
                q__1.i = fs.i / *c__; // , expr subst
                r__->r = q__1.r, r__->i = q__1.i;
                rtmax *= 2;
                if (f2 > rtmin && h2 < rtmax)
                {
                    /* safmin <= sqrt( f2*h2 ) <= safmax */
                    q__2.r = gs.r;
                    q__2.i = -gs.i;
                    r__1 = sqrt(f2 * h2);
                    q__3.r = fs.r / r__1;
                    q__3.i = fs.i / r__1; // , expr subst
                    q__1.r = q__2.r * q__3.r - q__2.i * q__3.i;
                    q__1.i = q__2.r * q__3.i + q__2.i * q__3.r; // , expr subst
                    s->r = q__1.r, s->i = q__1.i;
                }
                else
                {
                    q__2.r = gs.r;
                    q__2.i = -gs.i;
                    q__3.r = r__->r / h2;
                    q__3.i = r__->i / h2; // , expr subst
                    q__1.r = q__2.r * q__3.r - q__2.i * q__3.i;
                    q__1.i = q__2.r * q__3.i + q__2.i * q__3.r; // , expr subst
                    s->r = q__1.r, s->i = q__1.i;
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
                    q__1.r = fs.r / *c__;
                    q__1.i = fs.i / *c__; // , expr subst
                    r__->r = q__1.r, r__->i = q__1.i;
                }
                else
                {
                    /* f2 / sqrt(f2 * h2) < safmin, then */
                    /* sqrt(safmin) <= f2 * sqrt(safmax) <= h2 / sqrt(f2 * h2 */
                    r__1 = h2 / d__;
                    q__1.r = r__1 * fs.r;
                    q__1.i = r__1 * fs.i; // , expr subst
                    r__->r = q__1.r, r__->i = q__1.i;
                }
                q__2.r = gs.r;
                q__2.i = -gs.i;
                q__3.r = fs.r / d__;
                q__3.i = fs.i / d__; // , expr subst
                q__1.r = q__2.r * q__3.r - q__2.i * q__3.i;
                q__1.i = q__2.r * q__3.i + q__2.i * q__3.r; // , expr subst
                s->r = q__1.r, s->i = q__1.i;
            }
            /* Rescale c and r */
            *c__ *= w;
            q__1.r = u * r__->r;
            q__1.i = u * r__->i; // , expr subst
            r__->r = q__1.r, r__->i = q__1.i;
        }
    }
    return 0;
}
/* clartg_ */
