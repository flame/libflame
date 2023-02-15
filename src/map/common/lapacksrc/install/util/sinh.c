#include "FLA_f2c.h"

#ifdef __cplusplus
extern "C" {
#endif

    double r_sinh(real *x)
    {
        return( sinh(*x) );
    }
    double d_sinh(doublereal *x)
    {
        return( sinh(*x) );
    }

#ifdef __cplusplus
}
#endif
