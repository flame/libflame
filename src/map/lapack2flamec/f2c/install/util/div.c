#include "FLA_f2c.h"

#ifdef __cplusplus
extern "C" {
#endif

#ifdef _WIN32
    void c_div(complex *c, complex *a, complex *b)
    {
         c->r = (a->r * b->r + a->i * b->i) / (b->r * b->r + b->i * b->i);
         c->i = (a->i * b->r - a->r * b->i) / (b->r * b->r + b->i * b->i);
    }
    void z_div(doublecomplex *c, doublecomplex *a, doublecomplex *b)
    {
        c->r = (a->r * b->r + a->i * b->i) / (b->r * b->r + b->i * b->i);
        c->i = (a->i * b->r - a->r * b->i) / (b->r * b->r + b->i * b->i);
    }
#else
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
#endif

#ifdef __cplusplus
}
#endif
