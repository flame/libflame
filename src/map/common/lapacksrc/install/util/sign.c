#include "FLA_f2c.h"

#ifdef __cplusplus
extern "C" {
#endif

    double d_sign(doublereal *a, doublereal *b)
    {
        double x = (*a >= 0 ? *a : - *a);
        return ( *b >= 0 ? x : -x);
    }
    shortint h_sign(shortint *a, shortint *b)
    {
        shortint x = (*a >= 0 ? *a : - *a);
        return ( *b >= 0 ? x : -x);
    }
    integer i_sign(integer *a, integer *b)
    {
        integer x = (*a >= 0 ? *a : - *a);
        return ( *b >= 0 ? x : -x);
    }
    double r_sign(real *a, real *b)
    {
        double x = (*a >= 0 ? *a : - *a);
        return ( *b >= 0 ? x : -x);
    }

#ifdef __cplusplus
}
#endif
