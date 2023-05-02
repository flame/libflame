/* slartg.f -- translated by f2c (version 20190311). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static real c_b2 = 1.f;
/* Subroutine */
int slartg_(real *f, real *g, real *c__, real *s, real *r__)
{
    AOCL_DTL_TRACE_LOG_INIT
    /* Builtin functions */
    double sqrt(doublereal), r_sign(real *, real *);
    /* System generated locals */
    real r__1, r__2, r__3;
    /* Local variables */
    real d__, u, f1, g1, fs, gs, rtmin, rtmax, safmin, safmax;
    /* ...Translated by Pacific-Sierra Research vf90 Personal 3.4N3 05:33:20 1/24/23 */
    /* ...Switches: */
    /* .. */
    /* .. Local Scalars .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Constants .. */
    safmin = 1.1754943508222875e-38f;
    safmax = 8.5070591730234616e37f;
    rtmin = sqrt(safmin);
    rtmax = sqrt(safmax / 2);
    /* .. */
    /* .. Executable Statements .. */
    f1 = f2c_abs(*f);
    g1 = f2c_abs(*g);
    if (*g == 0.f)
    {
        *c__ = 1.f;
        *s = 0.f;
        *r__ = *f;
    }
    else if (*f == 0.f)
    {
        *c__ = 0.f;
        *s = r_sign(&c_b2, g);
        *r__ = g1;
    }
    else if (f1 > rtmin && f1 < rtmax && g1 > rtmin && g1 < rtmax)
    {
        d__ = sqrt(*f * *f + *g * *g);
        *c__ = f1 / d__;
        *r__ = r_sign(&d__, f);
        *s = *g / *r__;
    }
    else
    {
        /* Computing MIN */
        /* Computing MAX */
        r__3 = fla_max(safmin,f1);
        r__1 = safmax;
        r__2 = fla_max(r__3,g1); // , expr subst
        u = fla_min(r__1,r__2);
        fs = *f / u;
        gs = *g / u;
        d__ = sqrt(fs * fs + gs * gs);
        *c__ = f2c_abs(fs) / d__;
        *r__ = r_sign(&d__, f);
        *s = gs / *r__;
        *r__ *= u;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return 0;
}
/* slartg_ */
