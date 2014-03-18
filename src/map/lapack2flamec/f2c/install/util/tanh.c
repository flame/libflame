#include "FLA_f2c.h"

#ifdef __cplusplus
extern "C" {
#endif

    double d_tanh(doublereal *x)
    {
        return( tanh(*x) );
    }
    double r_tanh(real *x)
    {
        return( tanh(*x) );
    }

#ifdef __cplusplus
}
#endif
