/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

/*
    This file was originally authored by Kazushige Goto (TACC, UT-Austin).
    We use it with his permission.
*/

#include "FLAME.h"

double FLA_Clock_helper( void );

double FLA_Clock()
{
  return FLA_Clock_helper();
}

// --- Begin non-Windows build definitions -------------------------------------
#ifndef FLA_ENABLE_WINDOWS_BUILD

// This function is only used internally, which is why its prototype is not
// with the other "util" directory function prototypes.
void   detect_clocks( void );

// A global variable used when FLA_Clock_helper() is defined in terms of
// clock_gettime()/gettimeofday().
double gtod_ref_time_sec = 0.0;


// --- Begin portable FLA_Clock() definitions ----------------------------------
#ifdef FLA_ENABLE_PORTABLE_TIMER

double FLA_Clock_helper()
{
#ifdef FLA_ENABLE_TIDSP
  double the_time;
  the_time = 0.0;
#else
#ifdef FLA_PORTABLE_TIMER_IS_CLOCK_GETTIME
  double          the_time, norm_sec;
  struct timespec tsp;

  clock_gettime( CLOCK_REALTIME, &tsp );

  // If this is the first invocation of through FLA_Clock(), then initialize
  // the "reference time" global variable to the seconds field of the tv
  // struct.
  if ( gtod_ref_time_sec == 0.0 ) gtod_ref_time_sec = ( double ) tsp.tv_sec;

  // Normalize the seconds field of the tv struct so that it is relative to the
  // "reference time" that was recorded during the first invocation of
  // FLA_Clock().
  norm_sec = ( double ) tsp.tv_sec - gtod_ref_time_sec;

  // Compute the number of seconds since the reference time.
  the_time = norm_sec + tsp.tv_nsec * 1.0e-9;

#else
#ifdef FLA_PORTABLE_TIMER_IS_GETTIMEOFDAY
  
  double         the_time, norm_sec;
  struct timeval tv;

  gettimeofday( &tv, NULL );

  // If this is the first invocation of through FLA_Clock(), then initialize
  // the "reference time" global variable to the seconds field of the tv
  // struct.
  if ( gtod_ref_time_sec == 0.0 ) gtod_ref_time_sec = ( double ) tv.tv_sec;

  // Normalize the seconds field of the tv struct so that it is relative to the
  // "reference time" that was recorded during the first invocation of
  // FLA_Clock().
  norm_sec = ( double ) tv.tv_sec - gtod_ref_time_sec;

  // Compute the number of seconds since the reference time.
  the_time = norm_sec + tv.tv_usec * 1.0e-6;

#else //#ifdef FLA_PORTABLE_TIMER_IS_UNKNOWN

  the_time = 0.0;

#endif
#endif
#endif

  return the_time;
}

#else
// --- Begin non-portable FLA_Clock() definitions ------------------------------


// Global variables that are used only for non-portable FLA_Clock() definitions.
static double        clocks             = 0.0;
#ifdef __i386__
static unsigned int  initialclockoffset = 0;
#endif



// --- Begin ia64 section ------------------------------------------------------
#ifdef __ia64__


static __inline unsigned long rdtsc( void );
static __inline unsigned long rdtsc()
{
  unsigned long clocks;

#ifdef __INTEL_COMPILER
  clocks = __getReg(_IA64_REG_AR_ITC);  
#else
  __asm__ __volatile__ ("mov %0=ar.itc" : "=r"(clocks));
#endif

  return clocks;
}


double FLA_Clock_helper()
{
  unsigned long totalclocks;

  if ( clocks == 0.0 ) detect_clocks();
  totalclocks = rdtsc();

  return (double)totalclocks / clocks;
}


// --- End ia64 section --------------------------------------------------------
#else
// --- Begin i386 section ------------------------------------------------------
#ifdef __i386__

// rdtsc was once declared static. neccessary?
inline void rdtsc( unsigned int *high, unsigned int *low );
inline void rdtsc( unsigned int *high, unsigned int *low )
{
  asm("rdtsc" : "=a" (*low), "=d"(*high): : "cc");
}


double FLA_Clock_helper()
{
  unsigned int high, low;
  unsigned long long totalclocks;

  if (!clocks) detect_clocks();
  rdtsc(&high, &low);
  high -= initialclockoffset;
  totalclocks = (((unsigned long long)high) << 32)
                | (unsigned long long)low;

  return (double)totalclocks / clocks;
}


// --- End i386 section --------------------------------------------------------
#else
// --- Begin non-ia64, non-i386 section (default to portable code) -------------


double FLA_Clock_helper()
{
  double         the_time, norm_sec;
  struct timeval tv;

  gettimeofday( &tv, NULL );

  // If this is the first invocation of through FLA_Clock(), then initialize
  // the "reference time" global variable to the seconds field of the tv
  // struct.
  if ( gtod_ref_time_sec == 0.0 )
    gtod_ref_time_sec = ( double ) tv.tv_sec;

  // Normalize the seconds field of the tv struct so that it is relative to the
  // "reference time" that was recorded during the first invocation of
  // FLA_Clock().
  norm_sec = ( double ) tv.tv_sec - gtod_ref_time_sec;

  // Compute the number of seconds since the reference time.
  the_time = norm_sec + tv.tv_usec * 1.0e-6;

  return the_time;
}

#endif
#endif
// --- End non-ia64, non-i386 section ------------------------------------------
// --- Return to non-portable FLA_Clock() definitions --------------------------


// This implementation of detect_clocks() works only on operating systems that
// support the /proc filesystem (in particular, GNU/Linux).
// We place this function here at the end so that rdtsc() has been defined
// (and prototyped) by the time we get here. Otherwise the compiler spits out
// warnings.
void detect_clocks()
{
  FILE *infile;
  char buffer[256], *p;
#ifdef __i386__
  unsigned int high, low;
#endif

  if ( clocks == 0.0 )
  {
    p = (char *)NULL;
    infile = fopen("/proc/cpuinfo", "r");
    while (fgets(buffer, sizeof(buffer), infile))
    {
      if (!strncmp("cpu MHz", buffer, 6))
      {
        p = strchr(buffer, ':') + 1;
        break;
      }
    }
    clocks = 1.e6 * atof(p);
#ifdef __i386__
    rdtsc(&high, &low);
    initialclockoffset = high;
#endif
  }
}

// --- End non-portable FLA_Clock() definitions --------------------------------
#endif

// --- End non-Windows build definitions ---------------------------------------
// --- Begin Windows build definitions -----------------------------------------
#else

#define WIN32_LEAN_AND_MEAN
#define VC_EXTRALEAN
#include <windows.h>


double FLA_Clock_helper()
{
  LARGE_INTEGER clock_freq = {0};
  LARGE_INTEGER clock_val;
  BOOL          r_val;

  r_val = QueryPerformanceFrequency( &clock_freq );

  if ( r_val == 0 )
  {
    FLA_Print_message( "QueryPerformanceFrequency() failed", __FILE__, __LINE__ );
    FLA_Abort();
  }

  r_val = QueryPerformanceCounter( &clock_val );

  if ( r_val == 0 )
  {
    FLA_Print_message( "QueryPerformanceCounter() failed", __FILE__, __LINE__ );
    FLA_Abort();
  }

  return ( ( double) clock_val.QuadPart / ( double) clock_freq.QuadPart );
}

#endif
// --- End Windows build definitions -------------------------------------------
