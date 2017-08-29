/* dcabs1.f -- translated by f2c (version 19991025). You must link the resulting object file with the libraries: -lf2c -lm (in that order) */
#include "FLA_f2c.h"
doublereal dcabs1_(doublecomplex *z__)
{
    /* System generated locals */
    doublereal ret_val;
    static doublecomplex equiv_0[1];
    /* Local variables */
#define t ((doublereal *)equiv_0) 
#define zz (equiv_0)
    zz->r = z__->r, zz->i = z__->i;
    ret_val = f2c_abs(t[0]) + f2c_abs(t[1]);
    return ret_val;
}
/* dcabs1_ */
#undef zz
#undef t

