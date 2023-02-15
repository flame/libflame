#include "FLA_f2c.h"

extern
complex cdotc_(integer *n, complex *cx, integer *incx, complex *cy, integer *incy);

VOID cdotc_f2c_(complex *r, integer *n, complex *cx, integer *incx, complex *cy, integer *incy)
{
  *r = cdotc_(n, cx, incx, cy, incy);
}

extern
doublecomplex zdotc_(integer *n, doublecomplex *zx, integer *incx, doublecomplex *zy, integer *incy);
VOID zdotc_f2c_(doublecomplex *r, integer *n, doublecomplex *cx, integer *incx, doublecomplex *cy, integer *incy)
{
  *r = zdotc_(n, cx, incx, cy, incy);
}
