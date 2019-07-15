/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

extern __thread fla_axpyt_t* fla_axpyt_cntl_blas;

FLA_Error FLA_Axpyt( FLA_Trans trans, FLA_Obj alpha, FLA_Obj A, FLA_Obj B )
{
  FLA_Error r_val;

#ifdef FLA_ENABLE_BLAS1_FRONT_END_CNTL_TREES
  // Check parameters.
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Axpyt_check( trans, alpha, A, B );

  // Invoke FLA_Axpyt_internal() with flat control tree that simply calls
  // external wrapper.
  r_val = FLA_Axpyt_internal( trans, alpha, A, B, fla_axpyt_cntl_blas );

#else
  r_val = FLA_Axpyt_external( trans, alpha, A, B );
#endif

  return r_val;
}

