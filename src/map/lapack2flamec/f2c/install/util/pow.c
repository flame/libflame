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

#ifdef _WIN32
    /* Complex */
    void pow_ci(complex *p, complex *a, integer *b)
    {
        _Fcomplex a_ = {a->r, a->i};
        _Fcomplex b_ = {*b , 0 };
        _Fcomplex ret_val = cpowf(a_, b_);
        p->r = crealf(ret_val);
        p->i = cimagf(ret_val);
    }
    void pow_zz(doublecomplex *r, doublecomplex *a, doublecomplex *b)
    {
        _Dcomplex a_ = {a->r, a->i};
        _Dcomplex b_ = { a->r, a->i };
        _Dcomplex ret_val = cpow(a_, b_);
        r->r = creal(ret_val);
        r->i = cimag(ret_val);
    }
    void pow_zi(doublecomplex *p, doublecomplex *a, integer *b)
    {
        _Dcomplex a_ = {a->r, a->i};
        _Dcomplex b_ = {*b , 0 };
        _Dcomplex ret_val = cpow(a_, b_);
        p->r = creal(ret_val);
        p->i = cimag(ret_val);
    }
#else
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
#endif

#ifdef __cplusplus
}
#endif

