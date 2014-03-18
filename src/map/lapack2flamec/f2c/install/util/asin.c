#include "FLA_f2c.h"

#ifdef __cplusplus
extern "C" {
#endif

    double r_asin(real *x)
    {
        return( asin(*x) );
    }
    double d_asin(doublereal *x)
    {
        return( asin(*x) );
    }
    /*
    void c_asin(complex *r, complex *z)
    {
      double _Complex ret_val = casin(*z);
      r->r = creal(ret_val);
      r->i = cimag(ret_val);
    }
    void z_asin(doublecomplex *r, doublecomplex *z)
    {
      double _Complex ret_val = casin(*z);
      r->r = creal(ret_val);
      r->i = cimag(ret_val);
    }
    */
#ifdef __cplusplus
}
#endif
