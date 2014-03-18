#include "FLA_f2c.h"

#ifdef __cplusplus
extern "C" {
#endif

    double r_atan(real *x)
    {
        return( atan(*x) );
    }
    double d_atan(doublereal *x)
    {
        return( atan(*x) );
    }
    /*
    void c_atan(complex *r, complex *z)
    {
      double _Complex ret_val = catan(*z);
      r->r = creal(ret_val);
      r->i = cimag(ret_val);
    }
    void z_atan(doublecomplex *r, doublecomplex *z)
    {
      double _Complex ret_val = catan(*z);
      r->r = creal(ret_val);
      r->i = cimag(ret_val);
    }
    */
#ifdef __cplusplus
}
#endif
