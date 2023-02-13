#include "FLA_f2c.h"

#ifdef __cplusplus
extern "C" {
#endif

  integer i_dnnt(doublereal *x)
  {
    return (integer)(*x >= 0. ? floor(*x + .5) : -floor(.5 - *x));
  }

#ifdef __cplusplus
}
#endif
