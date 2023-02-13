#include "FLA_f2c.h"

#ifdef __cplusplus
extern "C" {
#endif

    void c_div(complex *c, complex *a, complex *b)
    {
        double _Complex ret_val = (a->r + I*a->i) / (b->r + I*b->i);
        c->r = creal(ret_val);
        c->i = cimag(ret_val);
    }
    void z_div(doublecomplex *c, doublecomplex *a, doublecomplex *b)
    {
        double _Complex ret_val = (a->r + I*a->i) / (b->r + I*b->i);
        c->r = creal(ret_val);
        c->i = cimag(ret_val);
    }

#ifdef __cplusplus
}
#endif
