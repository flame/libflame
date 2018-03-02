#include "FLA_f2c.h"
#include <sys/times.h>
//#include <sys/types.h>
//#include <time.h>

#ifndef CLK_TCK
#define CLK_TCK 60
#endif

real second_( void )
{
    struct tms rusage;

    times(&rusage);
    return (real)(rusage.tms_utime) / CLK_TCK;

} /* second_ */
