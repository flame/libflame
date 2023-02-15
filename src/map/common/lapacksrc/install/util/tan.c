#include "FLA_f2c.h"

#ifdef __cplusplus
extern "C" {
#endif

    double d_tan(doublereal *x)
    {
        return( tan(*x) );
    }
    double r_tan(real *x)
    {
        return( tan(*x) );
    }

#ifdef __cplusplus
}
#endif
