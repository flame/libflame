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

#ifdef _WIN32
    /* Complex */
    double c_abs(complex *z)
    {
        _Fcomplex z_ = {z->r, z->i};
        return  (cabsf(z_));
    }
    double z_abs(doublecomplex *z)
    {
        _Dcomplex z_ = { z->r, z->i };
        return  (cabs(z_));
    }
#else
    /* Complex */
    double c_abs(complex *z)
    {
      return  (cabs(z->r + I*z->i));
    }
    double z_abs(doublecomplex *z)
    {
      return  (cabs(z->r + I*z->i));
    }
#endif

#ifdef __cplusplus
}
#endif
