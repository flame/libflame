#include "FLA_f2c.h"

#ifdef __cplusplus
extern "C" {
#endif

    double r_sqrt(real *x)
    {
        return ( sqrt(*x) );
    }
    double d_sqrt(doublereal *x)
    {
        return ( sqrt(*x) );
    }
#ifdef _WIN32
    void c_sqrt(complex *r, complex *z)
    {
        _Dcomplex z_ = { z->r, z->i };
        _Dcomplex ret_val = csqrt(z_);
        r->r = creal(ret_val);
        r->i = cimag(ret_val);
    }
    void z_sqrt(doublecomplex *r, doublecomplex *z)
    {
        _Dcomplex z_ = { z->r, z->i };
        _Dcomplex ret_val = csqrt(z_);
        r->r = creal(ret_val);
        r->i = cimag(ret_val);
    }
#else
    void c_sqrt(complex *r, complex *z)
    {
        double _Complex ret_val = csqrt(z->r + I*z->i);
        r->r = creal(ret_val);
        r->i = cimag(ret_val);
    }
    void z_sqrt(doublecomplex *r, doublecomplex *z)
    {
        double _Complex ret_val = csqrt(z->r + I*z->i);
        r->r = creal(ret_val);
        r->i = cimag(ret_val);
    }
#endif

#ifdef __cplusplus
}
#endif
