/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#ifdef FLA_ENABLE_WINDOWS_BUILD
  #include <time.h>
#else
  // Handle the results of checking for time.h and sys/time.h
  #if TIME_WITH_SYS_TIME
    #include <sys/time.h>
    #include <time.h>
  #else
    #if HAVE_SYS_TIME_H
      #include <sys/time.h>
    #else
      #include <time.h>
    #endif
  #endif
#endif

// Handle the results of checking for ia64intrin.h. The contents of this header
// are required by the ia64 sections of FLA_Clock.c.
#ifdef HAVE_IA64INTRIN_H
  #include <ia64intrin.h>
#endif

