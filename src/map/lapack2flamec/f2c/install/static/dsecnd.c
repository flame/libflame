#include "FLA_f2c.h"
#ifdef _WIN32
#include <time.h>
#else
#include <sys/times.h>
#endif

#ifndef CLK_TCK
#define CLK_TCK 60
#endif

#ifdef _WIN32
doublereal dsecnd_( void )
{

    clock_t rusage = clock();
    return (doublereal)(rusage) / CLK_TCK;
}
#else
doublereal dsecnd_( void )
{
    struct tms rusage;

    times(&rusage);
    return (doublereal)(rusage.tms_utime) / CLK_TCK;

} 
#endif /* dsecnd_ */
