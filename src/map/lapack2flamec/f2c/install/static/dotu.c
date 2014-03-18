#include "FLA_f2c.h"


extern
complex cdotu_(integer *n, complex *cx, integer *incx, complex *cy, integer *incy);
VOID cdotu_f2c_(complex *r, integer *n, complex *cx, integer *incx, complex *cy, integer *incy)
{
  *r = cdotu_(n, cx, incx, cy, incy);
}

extern
doublecomplex zdotu_(integer *n, doublecomplex *zx, integer *incx, doublecomplex *zy, integer *incy);
VOID zdotu_f2c_(doublecomplex *r, integer *n, doublecomplex *cx, integer *incx, doublecomplex *cy, integer *incy)
{
  *r = zdotu_(n, cx, incx, cy, incy);
}

