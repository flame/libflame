#include "FLA_f2c.h"

#ifdef __cplusplus
extern "C" {
#endif

    double r_exp(real *x)
    {
        return( exp(*x) );
    }
    double d_exp(doublereal *x)
    {
        return( exp(*x) );
    }

    void c_exp(complex *r, complex *z)
    {
        double _Complex ret_val = cexp(z->r + I*z->i);
        r->r = creal(ret_val);
        r->i = cimag(ret_val);
    }
    void z_exp(doublecomplex *r, doublecomplex *z)
    {
        double _Complex ret_val = cexp(z->r + I*z->i);
        r->r = creal(ret_val);
        r->i = cimag(ret_val);
    }
#ifdef __cplusplus
}
#endif
