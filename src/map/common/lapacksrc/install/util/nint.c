#include "FLA_f2c.h"

#ifdef __cplusplus
extern "C" {
#endif

    double d_nint(doublereal *x)
    {
        return( (*x)>=0 ? floor(*x + .5) : -floor(.5 - *x) );
    }
    shortint h_nint(real *x)
    {
        return (shortint)(*x >= 0 ? floor(*x + .5) : -floor(.5 - *x));
    }
    integer i_nint(real *x)
    {
        return (integer)(*x >= 0 ? floor(*x + .5) : -floor(.5 - *x));
    }
    double r_nint(real *x)
    {
        return( (*x)>=0 ? floor(*x + .5) : -floor(.5 - *x) );
    }

#ifdef __cplusplus
}
#endif
