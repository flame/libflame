#include "FLA_f2c.h"

#ifdef __cplusplus
extern "C" {
#endif

    double r_sin(real *x)
    {
        return( sin(*x) );
    }
    double d_sin(doublereal *x)
    {
        return( sin(*x) );
    }
    void c_sin(complex *r, complex *z)
    {
        double _Complex ret_val = csin(z->r + I*z->i);
        r->r = creal(ret_val);
        r->i = cimag(ret_val);
    }
    void z_sin(doublecomplex *r, doublecomplex *z)
    {
        double _Complex ret_val = csin(z->r + I*z->i);
        r->r = creal(ret_val);
        r->i = cimag(ret_val);
    }

#ifdef __cplusplus
}
#endif
