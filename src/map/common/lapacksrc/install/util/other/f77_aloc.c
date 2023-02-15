#include "FLA_f2c.h"
#undef abs
#undef min
#undef max
#include "stdio.h"

static integer memfailure = 3;
extern void exit_(integer*);

#ifdef __cplusplus
extern "C" {
#endif

  char * F77_aloc(integer Len, const char *whence)
  {
    char *rv;
    unsigned int uLen = (unsigned int) Len;	/* for K&R C */
    
    if (!(rv = (char*)malloc(uLen)))
      {
        fprintf(stderr, "malloc(%u) failure in %s\n",
                uLen, whence);
        exit_(&memfailure);
      }
    return rv;
  }

#ifdef __cplusplus
}
#endif
