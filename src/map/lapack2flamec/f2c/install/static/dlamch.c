#include "FLA_f2c.h"
#include <float.h>

/* Table of constant values */

static doublereal half  = 0.5;
static doublereal one   = 1.0;
static doublereal zero  = 0.0;

doublereal dlamch_(char *cmach)
{
    /* Initialized data */
    static logical first = TRUE_;

    /* System generated locals */
    doublereal ret_val;

    /* Local variables */
    static doublereal  eps, sfmin, base, prec, t, rnd, emin, rmin, emax, rmax;
    doublereal rmach, small;

    extern logical lsame_(char *, char *);

    /*  Purpose */
    /*  ======= */

    /*  DLAMCH determines single precision machine parameters. */

    /*  Arguments */
    /*  ========= */

    /*  CMACH   (input) CHARACTER*1 */
    /*          Specifies the value to be returned by DLAMCH: */
    /*          = 'E' or 'e',   DLAMCH := eps */
    /*          = 'S' or 's ,   DLAMCH := sfmin */
    /*          = 'B' or 'b',   DLAMCH := base */
    /*          = 'P' or 'p',   DLAMCH := eps*base */
    /*          = 'N' or 'n',   DLAMCH := t */
    /*          = 'R' or 'r',   DLAMCH := rnd */
    /*          = 'M' or 'm',   DLAMCH := emin */
    /*          = 'U' or 'u',   DLAMCH := rmin */
    /*          = 'L' or 'l',   DLAMCH := emax */
    /*          = 'O' or 'o',   DLAMCH := rmax */

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
        /* FLT_ROUNDS specification
        -1 undetermined
         0 toward zero
         1 to nearest
         2 toward positive infinity
         3 toward negative infinity
            */
        if (FLT_ROUNDS == 1)
        {
            rnd = one;
            eps = DBL_EPSILON * half;
        }
        else
        {
            rnd = zero;
            eps = DBL_EPSILON;
        }
        base  = FLT_RADIX;
        prec  = eps * base;
        sfmin = DBL_MIN;
        small = one / DBL_MAX;
        if ( small >= sfmin)
            sfmin = small * (one + eps);

		// For t, we need the number of base-2 digits, not base-10 digits.
		// Here, we hardcode the value obtained from netlib LAPACK.
        //t    = DBL_DIG;
        //t    = 53;
        t    = DBL_MANT_DIG;
        emin = DBL_MIN_EXP;
        emax = DBL_MAX_EXP;
        rmin = DBL_MIN;
        rmax = DBL_MAX;
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

} /* dlamch_ */


doublereal dlamc3_(doublereal *a, doublereal *b)
{
    /* System generated locals */
    doublereal ret_val;


    /*  -- LAPACK auxiliary routine (version 3.4.0) -- */
    /*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd.. */
    /*     November 2010 */

    /*     .. Scalar Arguments .. */
    /*     .. */

    /*  Purpose */
    /*  ======= */

    /*  DLAMC3  is intended to force  A  and  B  to be stored prior to doing */
    /*  the addition of  A  and  B ,  for use in situations where optimizers */
    /*  might hold one of these in a register. */

    /*  Arguments */
    /*  ========= */

    /*  A       (input) DOUBLE PRECISION */
    /*  B       (input) DOUBLE PRECISION */
    /*          The values A and B. */

    /* ===================================================================== */

    /*     .. Executable Statements .. */

    ret_val = *a + *b;

    return ret_val;

    /*     End of DLAMC3 */

} /* dlamc3_ */

