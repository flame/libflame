
#include "FLAME.h"

extern fla_apq2ut_t* flash_apq2ut_cntl;
extern fla_apq2ut_t* flash_apq2ut_cntl_leaf;
extern fla_apq2ut_t* fla_apq2ut_cntl_leaf;

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

