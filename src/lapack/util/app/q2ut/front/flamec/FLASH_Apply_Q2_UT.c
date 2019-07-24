/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

extern TLS_CLASS_SPEC fla_apq2ut_t* flash_apq2ut_cntl;
extern TLS_CLASS_SPEC fla_apq2ut_t* flash_apq2ut_cntl_leaf;
extern TLS_CLASS_SPEC fla_apq2ut_t* fla_apq2ut_cntl_leaf;

FLA_Error FLASH_Apply_Q2_UT( FLA_Side side, FLA_Trans trans, FLA_Direct direct, FLA_Store storev,
                             FLA_Obj D, FLA_Obj T, FLA_Obj W, FLA_Obj C,
                                                              FLA_Obj E )
{
  FLA_Error r_val;

  // Check parameters.
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Apply_Q2_UT_check( side, trans, direct, storev, D, T, W, C, E );

  // Begin a parallel region.
  FLASH_Queue_begin();
  
  // Invoke FLA_Apply_Q2_UT_internal() with the standard control tree.
  r_val = FLA_Apply_Q2_UT_internal( side, trans, direct, storev, D, T, W, C, E, flash_apq2ut_cntl );

  // End the parallel region.
  FLASH_Queue_end();

  return r_val;
}

