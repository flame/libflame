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

#ifdef _WIN32
    void c_cos(complex *r, complex *z)
    {
        _Dcomplex z_ = { z->r, z->i };
        _Dcomplex ret_val = ccos(z_);
        r->r = creal(ret_val);
        r->i = cimag(ret_val);
    }
    void z_cos(doublecomplex *r, doublecomplex *z)
    {
        _Dcomplex z_ = { z->r, z->i };
        _Dcomplex ret_val = ccos(z_);
        r->r = creal(ret_val);
        r->i = cimag(ret_val);
    }
#else
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
#endif

#ifdef __cplusplus
}
#endif
