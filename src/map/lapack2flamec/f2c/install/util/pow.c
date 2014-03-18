#include "FLA_f2c.h"

#ifdef __cplusplus
extern "C" {
#endif

    /* Integer */
    shortint pow_hh(shortint *ap, shortint *bp)
    {
        return (shortint)(pow(*ap, *bp));
    }
    integer pow_ii(integer *ap, integer *bp)
    {
        return (integer)(pow(*ap, *bp));
    }
#ifdef INTEGER_STAR_8
    longint pow_qq(longint *ap, longint *bp)
    {
        return (longint)(pow(*ap, *bp));
    }
#endif

    /* Double */
    double pow_ri(real *ap, integer *bp)
    {
        return (pow(*ap, *bp));
    }
    double pow_dd(doublereal *ap, doublereal *bp)
    {
        return (pow(*ap, *bp));
    }
    double pow_di(doublereal *ap, integer *bp)
    {
        return (pow(*ap, *bp));
    }

    /* Complex */
    void pow_ci(complex *p, complex *a, integer *b)
    {
        double _Complex ret_val = cpow(a->r + I*a->i, *b);
        p->r = creal(ret_val);
        p->i = cimag(ret_val);
    }
    void pow_zz(doublecomplex *r, doublecomplex *a, doublecomplex *b)
    {
        double _Complex ret_val = cpow(a->r + I*a->i, b->r + I*b->i);
        r->r = creal(ret_val);
        r->i = cimag(ret_val);
    }
    void pow_zi(doublecomplex *p, doublecomplex *a, integer *b)
    {
        double _Complex ret_val = cpow(a->r + I*a->i, *b);
        p->r = creal(ret_val);
        p->i = cimag(ret_val);
    }


#ifdef __cplusplus
}
#endif

