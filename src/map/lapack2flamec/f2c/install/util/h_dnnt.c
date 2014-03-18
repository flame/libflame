#include "FLA_f2c.h"

#ifdef __cplusplus
extern "C" {
#endif

  shortint h_dnnt(doublereal *x)
  {
    return (shortint)(*x >= 0. ? floor(*x + .5) : -floor(.5 - *x));
  }

#ifdef __cplusplus
}
#endif
