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

#ifdef _WIN32
    void c_exp(complex *r, complex *z)
    {
        _Fcomplex z_ = { z->r, z->i };
        _Fcomplex ret_val = cexpf(z_);
        r->r = crealf(ret_val);
        r->i = cimagf(ret_val);
    }
    void z_exp(doublecomplex *r, doublecomplex *z)
    {
        _Dcomplex z_ = { z->r, z->i };
        _Dcomplex ret_val = cexp(z_);
        r->r = creal(ret_val);
        r->i = cimag(ret_val);
    }
#else
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
#endif

#ifdef __cplusplus
}
#endif
