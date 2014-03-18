#include "FLA_f2c.h"
 void r_cnjg(complex *dest, complex *src) {
 dest->r = src->r ;
 dest->i = -(src->i);
 }
 
