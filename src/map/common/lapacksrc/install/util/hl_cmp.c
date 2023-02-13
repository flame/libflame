#include "FLA_f2c.h"

#ifdef __cplusplus
extern "C" {
#endif

    extern integer s_cmp(char *, char *, ftnlen, ftnlen);

    shortlogical hl_ge(char *a, char *b, ftnlen la, ftnlen lb)
    {
        return(s_cmp(a,b,la,lb) >= 0);
    }
    shortlogical hl_gt(char *a, char *b, ftnlen la, ftnlen lb)
    {
        return(s_cmp(a,b,la,lb) > 0);
    }
    shortlogical hl_le(char *a, char *b, ftnlen la, ftnlen lb)
    {
        return(s_cmp(a,b,la,lb) <= 0);
    }
    shortlogical hl_lt(char *a, char *b, ftnlen la, ftnlen lb)
    {
        return(s_cmp(a,b,la,lb) < 0);
    }

#ifdef __cplusplus
}
#endif
