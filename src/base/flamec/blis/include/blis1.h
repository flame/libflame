/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#ifndef BLIS1_H
#define BLIS1_H

// Allow C++ users to include this header file in their source code. However,
// we make the extern "C" conditional on whether we're using a C++ compiler,
// since regular C compilers don't understand the extern "C" construct.
#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Determine whether or not we are using BLIS from libflame.
//#define BLIS1_FROM_LIBFLAME

#ifdef BLIS1_FROM_LIBFLAME

  // If using libflame, pull in its header files so that
  // vector intrinsics-related macro constants are set properly.
  //#include "FLAME.h"
  #include "FLA_config.h"
  #include "FLA_macro_defs.h"
  #include "FLA_type_defs.h"

  // --- Pass-through macros for BLIS ---
  #ifdef FLA_ENABLE_CBLAS_INTERFACES
    #define BLIS1_ENABLE_CBLAS_INTERFACES
  #endif
  #ifdef FLA_ENABLE_WINDOWS_BUILD
    #define BLIS1_ENABLE_WINDOWS_BUILD
  #endif
  #ifdef FLA_ENABLE_UPPERCASE_F77
    #define BLIS1_ENABLE_UPPERCASE_F77
  #endif
  #ifdef FLA_ENABLE_VECTOR_INTRINSICS
    #define BLIS1_ENABLE_VECTOR_INTRINSICS
  #endif

  #define BLIS1_VECTOR_INTRINSIC_TYPE FLA_VECTOR_INTRINSIC_TYPE

#else

  // --- BLIS configuration options ---

  // #define BLIS1_ENABLE_USE_OF_FLA_MALLOC
  // #define BLIS1_ENABLE_CBLAS_INTERFACES
  // #define BLIS1_ENABLE_WINDOWS_BUILD
  // #define BLIS1_ENABLE_UPPERCASE_F77
  // #define BLIS1_ENABLE_VECTOR_INTRINSICS
  //   #define BLIS1_VECTOR_INTRINSIC_TYPE BLIS1_NO_INTRINSICS
  //   #define BLIS1_VECTOR_INTRINSIC_TYPE BLIS1_SSE_INTRINSICS

#endif

#include "blis_macro_defs.h"
#include "blis_type_defs.h"

#include "blis_prototypes_util.h"
#include "blis_prototypes_query.h"
#include "blis_prototypes_misc.h"

#include "blis_prototypes_level1.h"
#include "blis_prototypes_level2.h"
#include "blis_prototypes_level3.h"

#include "blis_prototypes_fused1.h"

#include "blis_f77_name_mangling.h"

#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
  #include "blis_prototypes_cblas.h"
#else
  #include "blis_prototypes_blas.h"
#endif

// End extern "C" construct block.
#ifdef __cplusplus
}
#endif

#endif
