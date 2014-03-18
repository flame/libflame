#include "FLA_f2c.h"

#ifdef __cplusplus
extern "C" {
#endif
    /* Integer */
    shortint h_abs(shortint *x)
    {
        return ( *x >= 0 ? (*x) : (- *x) );
    }
    integer i_abs(integer *x)
    {
        return ( *x >= 0 ? (*x) : (- *x) );
    }

    /* Double */
    double r_abs(real *x)
    {
        return ( *x >= 0 ? (*x) : (- *x) );
    }
    double d_abs(doublereal *x)
    {
        return ( *x >= 0 ? (*x) : (- *x) );
    }

    /* Complex */
    double c_abs(complex *z)
    {
      return  (cabs(z->r + I*z->i));
    }
    double z_abs(doublecomplex *z)
    {
      return  (cabs(z->r + I*z->i));
    }

#ifdef __cplusplus
}
#endif
