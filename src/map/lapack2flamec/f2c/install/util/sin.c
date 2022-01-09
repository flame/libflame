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
#ifdef _WIN32
    void c_sin(complex *r, complex *z)
    {
        _Fcomplex z_ = { z->r, z->i };
        _Fcomplex ret_val = csinf(z_);
        r->r = crealf(ret_val);
        r->i = cimagf(ret_val);
    }
    void z_sin(doublecomplex *r, doublecomplex *z)
    {
        _Dcomplex z_ = { z->r, z->i };
        _Dcomplex ret_val = csin(z_);
        r->r = creal(ret_val);
        r->i = cimag(ret_val);
    }
#else
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
#endif

#ifdef __cplusplus
}
#endif
