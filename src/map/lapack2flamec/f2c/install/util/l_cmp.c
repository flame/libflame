#include "FLA_f2c.h"

#ifdef __cplusplus
extern "C" {
#endif

    extern integer s_cmp(char *a0, char *b0, ftnlen la, ftnlen lb);

    logical l_ge(char *a, char *b, ftnlen la, ftnlen lb)
    {
        return(s_cmp(a,b,la,lb) >= 0);
    }
    logical l_gt(char *a, char *b, ftnlen la, ftnlen lb)
    {
        return(s_cmp(a,b,la,lb) > 0);
    }
    logical l_le(char *a, char *b, ftnlen la, ftnlen lb)
    {
        return(s_cmp(a,b,la,lb) <= 0);
    }
    logical l_lt(char *a, char *b, ftnlen la, ftnlen lb)
    {
        return(s_cmp(a,b,la,lb) < 0);
    }

#ifdef __cplusplus
}
#endif
