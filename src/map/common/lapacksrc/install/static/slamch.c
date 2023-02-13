#include "FLA_f2c.h"
#include <float.h>

/* Table of constant values */

static const real half  = 0.5f;
static const real one   = 1.f;
static const real zero  = 0.f;

real slamch_(char *cmach)
{
    /* Initialized data */
    static logical first = TRUE_;

    /* System generated locals */
    real ret_val;

    /* Local variables */
    static real  eps, sfmin, base, prec, t, rnd, emin, rmin, emax, rmax;
    real rmach, small;

    extern logical lsame_(char *, char *);

    /*  Purpose */
    /*  ======= */

    /*  SLAMCH determines single precision machine parameters. */

    /*  Arguments */
    /*  ========= */

    /*  CMACH   (input) CHARACTER*1 */
    /*          Specifies the value to be returned by SLAMCH: */
    /*          = 'E' or 'e',   SLAMCH := eps */
    /*          = 'S' or 's ,   SLAMCH := sfmin */
    /*          = 'B' or 'b',   SLAMCH := base */
    /*          = 'P' or 'p',   SLAMCH := eps*base */
    /*          = 'N' or 'n',   SLAMCH := t */
    /*          = 'R' or 'r',   SLAMCH := rnd */
    /*          = 'M' or 'm',   SLAMCH := emin */
    /*          = 'U' or 'u',   SLAMCH := rmin */
    /*          = 'L' or 'l',   SLAMCH := emax */
    /*          = 'O' or 'o',   SLAMCH := rmax */

    /*          where */

    /*          eps   = relative machine precision */
    /*          sfmin = safe minimum, such that 1/sfmin does not overflow */
    /*          base  = base of the machine */
    /*          prec  = eps*base */
    /*          t     = number of (base) digits in the mantissa */
    /*          rnd   = 1.0 when rounding occurs in addition, 0.0 otherwise */
    /*          emin  = minimum exponent before (gradual) underflow */
    /*          rmin  = underflow threshold - base**(emin-1) */
    /*          emax  = largest exponent before overflow */
    /*          rmax  = overflow threshold  - (base**emax)*(1-eps) */

    /* ===================================================================== */

    /* Assume rounding, not chopping. Always. -- This is a comment from LAPACK.
     */
    if (first)
    {
        /* flt_ROUNDS specification
        -1 undetermined
         0 toward zero
         1 to nearest
         2 toward positive infinity
         3 toward negative infinity
            */
        if (FLT_ROUNDS == 1)
        {
            rnd = one;
            eps = FLT_EPSILON * half;
        }
        else
        {
            rnd = zero;
            eps = FLT_EPSILON;
        }
        base  = FLT_RADIX;
        prec  = eps * base;
        sfmin = FLT_MIN;
        small = one / FLT_MAX;
        if ( small >= sfmin)
            sfmin = small * (one + eps);

		// For t, we need the number of base-2 digits, not base-10 digits.
		// Here, we hardcode the value obtained from netlib LAPACK.
        //t    = FLT_DIG;
        //t    = 24;
        t    = FLT_MANT_DIG;
        emin = FLT_MIN_EXP;
        emax = FLT_MAX_EXP;
        rmin = FLT_MIN;
        rmax = FLT_MAX;
    }

    if (lsame_(cmach, "E"))
    {
        rmach = eps;
    }
    else if (lsame_(cmach, "S"))
    {
        rmach = sfmin;
    }
    else if (lsame_(cmach, "B"))
    {
        rmach = base;
    }
    else if (lsame_(cmach, "P"))
    {
        rmach = prec;
    }
    else if (lsame_(cmach, "N"))
    {
        rmach = t;
    }
    else if (lsame_(cmach, "R"))
    {
        rmach = rnd;
    }
    else if (lsame_(cmach, "M"))
    {
        rmach = emin;
    }
    else if (lsame_(cmach, "U"))
    {
        rmach = rmin;
    }
    else if (lsame_(cmach, "L"))
    {
        rmach = emax;
    }
    else if (lsame_(cmach, "O"))
    {
        rmach = rmax;
    }

    ret_val = rmach;
    first = FALSE_;
    return ret_val;

} /* slamch_ */

real slamc3_(real *a, real *b)
{
    /* System generated locals */
    real ret_val;


    /*  -- LAPACK auxiliary routine (version 3.4.0) -- */
    /*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd.. */
    /*     November 2010 */

    /*     .. Scalar Arguments .. */
    /*     .. */

    /*  Purpose */
    /*  ======= */

    /*  SLAMC3  is intended to force  A  and  B  to be stored prior to doing */
    /*  the addition of  A  and  B ,  for use in situations where optimizers */
    /*  might hold one of these in a register. */

    /*  Arguments */
    /*  ========= */

    /*  A       (input) REAL */
    /*  B       (input) REAL */
    /*          The values A and B. */

    /* ===================================================================== */

    /*     .. Executable Statements .. */

    ret_val = *a + *b;

    return ret_val;

    /*     End of SLAMC3 */

} /* slamc3_ */

