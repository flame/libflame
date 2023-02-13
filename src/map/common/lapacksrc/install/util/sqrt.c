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

#ifdef __cplusplus
}
#endif
