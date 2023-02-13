#include "FLA_f2c.h"

#ifdef __cplusplus
extern "C" {
#endif

    VOID d_cnjg(doublecomplex *r, doublecomplex *z)
    {
        doublereal zi = z->i;
        r->r = z->r;
        r->i = -zi;
    }
    VOID r_cnjg(complex *r, complex *z)
    {
        real zi = z->i;
        r->r = z->r;
        r->i = -zi;
    }

#ifdef __cplusplus
}
#endif
