
#include "FLAME.h"

extern fla_hessut_t* fla_hessut_cntl_leaf;

FLA_Error FLA_Hess_UT( FLA_Obj A, FLA_Obj T )
{
  FLA_Error r_val;

  // Check parameters.
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Hess_UT_check( A, T );

  // Invoke FLA_Hess_UT_internal() with the standard control tree.
  r_val = FLA_Hess_UT_internal( A, T, fla_hessut_cntl_leaf );

  return r_val;
}

