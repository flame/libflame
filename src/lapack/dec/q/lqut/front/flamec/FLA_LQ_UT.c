
#include "FLAME.h"

extern fla_lqut_t*  fla_lqut_cntl_leaf;

FLA_Error FLA_LQ_UT( FLA_Obj A, FLA_Obj T )
{
  FLA_Error r_val;

  // Check parameters.
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_LQ_UT_check( A, T );

  // Invoke FLA_LQ_UT_internal() with the standard control tree.
  //r_val = FLA_LQ_UT_blk_var1( A, T, fla_lqut_cntl_leaf );
  r_val = FLA_LQ_UT_internal( A, T, fla_lqut_cntl_leaf );

  return r_val;
}

