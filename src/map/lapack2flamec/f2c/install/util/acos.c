#include "FLA_f2c.h"

#ifdef __cplusplus
extern "C" {
#endif

    double r_acos(real *x)
    {
        return( acos(*x) );
    }
    double d_acos(doublereal *x)
    {
        return( acos(*x) );
    }
    /*
    void c_acos(complex *r, complex *z)
    {
      double _Complex ret_val = cacos(*z);
      r->r = creal(ret_val);
      r->i = cimag(ret_val);
    }
    void z_acos(doublecomplex *r, doublecomplex *z)
    {
      double _Complex ret_val = cacos(*z);
      r->r = creal(ret_val);
      r->i = cimag(ret_val);
    }
    */
#ifdef __cplusplus
}
#endif
