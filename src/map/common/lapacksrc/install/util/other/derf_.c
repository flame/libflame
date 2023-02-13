#include "FLA_f2c.h"
#ifdef __cplusplus
extern "C" {
#endif

#ifdef KR_headers
double erf();
double derf_(x) doublereal *x;
#else
extern double erf(double);
double derf_(doublereal *x)
#endif
{
return( erf(*x) );
}
#ifdef __cplusplus
}
#endif
