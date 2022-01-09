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
real second_( void )
{
    clock_t rusage = clock();
    return (real)(rusage) / CLK_TCK;
}
#else
real second_( void )
{
    struct tms rusage;

    times(&rusage);
    return (real)(rusage.tms_utime) / CLK_TCK;

} 
#endif/* second_ */
