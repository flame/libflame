#include "FLA_f2c.h"
 void d_cnjg(doublecomplex *dest, doublecomplex *src) {
 dest->r = src->r ;
 dest->i = -(src->i);
 }
 
