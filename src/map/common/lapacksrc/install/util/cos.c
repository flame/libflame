#include "FLA_f2c.h"

#ifdef __cplusplus
extern "C" {
#endif

    double r_cos(real *x)
    {
        return( cos(*x) );
    }
    double d_cos(doublereal *x)
    {
        return( cos(*x) );
    }
    void c_cos(complex *r, complex *z)
    {
        double _Complex ret_val = ccos(z->r + I*z->i);
        r->r = creal(ret_val);
        r->i = cimag(ret_val);
    }
    void z_cos(doublecomplex *r, doublecomplex *z)
    {
        double _Complex ret_val = ccos(z->r + I*z->i);
        r->r = creal(ret_val);
        r->i = cimag(ret_val);
    }

#ifdef __cplusplus
}
#endif
