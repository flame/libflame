/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

// --- Define Fortran name-mangling macro --------------------------------------

// If the F77_FUNC name-mangling macro is undefined, then we we need to define
// it ourselves.
#ifndef F77_FUNC

  // Case 1: F77_FUNC is undefined because we're building for Windows.
  #ifdef BLIS1_ENABLE_WINDOWS_BUILD

    // Check whether we need to use uppercase Fortran routine names; otherwise
    // default to lowercase.
    #ifdef BLIS1_ENABLE_UPPERCASE_F77

      // Use uppercase routine names (no underscore).
      #define F77_FUNC( name_lower, name_upper ) \
              name_upper
    #else

      // Use lowercase routine names (no underscore).
      #define F77_FUNC( name_lower, name_upper ) \
              name_lower
    #endif

  // Case 2: F77_FUNC is undefined because we're in a Linux-like environment
  // that did not define it for us.
  #else

    // Check whether we need to use uppercase Fortran routine names; otherwise
    // default to lowercase.
    #ifdef BLIS1_ENABLE_UPPERCASE_F77

      // Use uppercase routine names (single underscore).
      #define F77_FUNC( name_lower, name_upper ) \
              name_upper ## _
    #else

      // Use lowercase routine names (single underscore).
      #define F77_FUNC( name_lower, name_upper ) \
              name_lower ## _
    #endif

  #endif // #ifdef BLIS1_ENABLE_WINDOWS_BUILD

#endif // #ifndef F77_FUNC

