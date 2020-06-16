#include "FLA_f2c.h"

#ifdef __cplusplus
extern "C" {
#endif

    double r_log(real *x)
    {
        return( log(*x) );
    }
    double d_log(doublereal *x)
    {
        return( log(*x) );
    }

#ifdef _WIN32
    void c_log(complex *r, complex *z)
    {
        _Dcomplex z_ = { z->r, z->i };
        _Dcomplex ret_val = clog(z_);
        r->r = creal(ret_val);
        r->i = cimag(ret_val);
    }
    void z_log(doublecomplex *r, doublecomplex *z)
    {
        _Dcomplex z_ = { z->r, z->i };
        _Dcomplex ret_val = clog(z_);
        r->r = creal(ret_val);
        r->i = cimag(ret_val);
    }
#else
    void c_log(complex *r, complex *z)
    {
        double _Complex ret_val = clog(z->r + I*z->i);
        r->r = creal(ret_val);
        r->i = cimag(ret_val);
    }
    void z_log(doublecomplex *r, doublecomplex *z)
    {
        double _Complex ret_val = clog(z->r + I*z->i);
        r->r = creal(ret_val);
        r->i = cimag(ret_val);
    }
#endif

#ifdef __cplusplus
}
#endif
