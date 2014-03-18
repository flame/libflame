#include "FLA_f2c.h"

#ifdef __cplusplus
extern "C" {
#endif

    double d_imag(doublecomplex *z)
    {
        return(z->i);
    }
    double r_imag(complex *z)
    {
        return(z->i);
    }

#ifdef __cplusplus
}
#endif
