#include "FLA_f2c.h"

#ifdef __cplusplus
extern "C" {
#endif

    double d_int(doublereal *x)
    {
        return( (*x>0) ? floor(*x) : -floor(- *x) );
    }
    double r_int(real *x)
    {
        return( (*x>0) ? floor(*x) : -floor(- *x) );
    }

#ifdef __cplusplus
}
#endif
