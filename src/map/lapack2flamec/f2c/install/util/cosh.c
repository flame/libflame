#include "FLA_f2c.h"

#ifdef __cplusplus
extern "C" {
#endif

    double r_cosh(real *x)
    {
        return( cosh(*x) );
    }
    double d_cosh(doublereal *x)
    {
        return( cosh(*x) );
    }

#ifdef __cplusplus
}
#endif
