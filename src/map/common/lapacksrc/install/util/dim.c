#include "FLA_f2c.h"

#ifdef __cplusplus
extern "C" {
#endif

    double d_dim(doublereal *a, doublereal *b)
    {
        return( *a > *b ? *a - *b : 0);
    }
    shortint h_dim(shortint *a, shortint *b)
    {
        return( *a > *b ? *a - *b : 0);
    }
    integer i_dim(integer *a, integer *b)
    {
        return( *a > *b ? *a - *b : 0);
    }
    double r_dim(real *a, real *b)
    {
        return( *a > *b ? *a - *b : 0);
    }

#ifdef __cplusplus
}
#endif

