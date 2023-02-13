#include "FLA_f2c.h"

#ifdef __cplusplus
extern "C" {
#endif

    double d_atn2(doublereal *x, doublereal *y)
    {
        return( atan2(*x,*y) );
    }
    double r_atn2(real *x, real *y)
    {
        return( atan2(*x,*y) );
    }

#ifdef __cplusplus
}
#endif
