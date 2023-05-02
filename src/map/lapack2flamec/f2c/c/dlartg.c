/* dlartg.f -- translated by f2c (version 20190311). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static doublereal c_b2 = 1.;
/* Subroutine */
int dlartg_(doublereal *f, doublereal *g, doublereal *c__, doublereal *s, doublereal *r__)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("dlartg inputs : f %lf, g %lf", *f, *g);
    /* Builtin functions */
    double sqrt(doublereal), d_sign(doublereal *, doublereal *);
    /* Local variables */
    doublereal d__, u, f1, g1, fs, gs, rtmin, rtmax, safmin, safmax;
    doublereal d__1, d__2, d__3;
    /* ...Translated by Pacific-Sierra Research vf90 Personal 3.4N3 05:19:29 1/25/23 */
    /* ...Switches: */
    /* .. */
    /* .. Local Scalars .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Constants .. */
    safmin = 2.2250738585072014e-308;
    safmax = 4.4942328371557898e307;
    rtmin = sqrt(safmin);
    rtmax = sqrt(safmax / 2);
    /* .. */
    /* .. Executable Statements .. */
    f1 = f2c_dabs(*f);
    g1 = f2c_dabs(*g);
    if (*g == 0.)
    {
        *c__ = 1.;
        *s = 0.;
        *r__ = *f;
    }
    else if (*f == 0.)
    {
        *c__ = 0.;
        *s = d_sign(&c_b2, g);
        *r__ = g1;
    }
    else if (f1 > rtmin && f1 < rtmax && g1 > rtmin && g1 < rtmax)
    {
        d__ = sqrt(*f * *f + *g * *g);
        *c__ = f1 / d__;
        *r__ = d_sign(&d__, f);
        *s = *g / *r__;
    }
    else
    {
        /* Computing MIN */
        /* Computing MAX */
        d__3 = fla_max(safmin,f1);
        d__1 = safmax;
        d__2 = fla_max(d__3,g1); // , expr subst
        u = fla_min(d__1,d__2);
        fs = *f / u;
        gs = *g / u;
        d__ = sqrt(fs * fs + gs * gs);
        *c__ = f2c_dabs(fs) / d__;
        *r__ = d_sign(&d__, f);
        *s = gs / *r__;
        *r__ *= u;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return 0;
}
/* dlartg_ */
