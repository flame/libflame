#include "FLA_f2c.h"
 void z_div(doublecomplex *cp, doublecomplex *ap, doublecomplex *bp) {
 doublecomplex a = *ap;
 doublecomplex b = *bp;
 double temp;
 temp = b.r * b.r + b.i * b.i;
 cp->r = ( a.r * b.r + a.i * b.i ) / temp;
 cp->i = ( a.i * b.r - a.r * b.i ) / temp;
 }
 
