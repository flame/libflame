#include "FLA_f2c.h"


#ifdef FLA_ENABLE_YET_LAPACK2FLAME

int s_kat(char *lp, char *rpp[], integer rnp[], integer *np)
{
    extern int s_cat(char *, char **, integer *, integer *, ftnlen);

    ftnlen ll;
    ll = strlen(lp);
    s_cat(lp, rpp, rnp, np, ll);
    return 0;
}

#endif
