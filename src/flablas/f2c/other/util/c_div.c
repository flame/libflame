#include "FLA_f2c.h"
 void c_div(complex *cp, complex *ap, complex *bp) {
 complex a = *ap;
 complex b = *bp;
 real temp;
 temp = b.r * b.r + b.i * b.i;
 cp->r = ( a.r * b.r + a.i * b.i ) / temp;
 cp->i = ( a.i * b.r - a.r * b.i ) / temp;
 }
 
