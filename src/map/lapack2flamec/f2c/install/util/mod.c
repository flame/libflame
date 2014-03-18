#include "FLA_f2c.h"


#ifdef __cplusplus
extern "C" {
#endif

    shortint h_mod(short *a, short *b)
    {
        return( *a % *b);
    }
    integer i_mod(integer *a, integer *b)
    {
        return( *a % *b);
    }
    double r_mod(real *x, real *y)
    {
        double quotient;
        if( (quotient = (double)*x / *y) >= 0)
            quotient = floor(quotient);
        else
            quotient = -floor(-quotient);

        return(*x - (*y) * quotient );
    }
    double d_mod(doublereal *x, doublereal *y)
    {
        double quotient;
        if( (quotient = *x / *y) >= 0)
            quotient = floor(quotient);
        else
            quotient = -floor(-quotient);
        return(*x - (*y) * quotient );
    }


#ifdef __cplusplus
}
#endif
