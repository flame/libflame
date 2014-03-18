#include "FLA_f2c.h"

#ifdef __cplusplus
extern "C" {
#endif

    double d_lg10(doublereal *x)
    {
        return( log10(*x) );
    }
    double r_lg10(real *x)
    {
        return( log10(*x) );
    }

#ifdef __cplusplus
}
#endif
