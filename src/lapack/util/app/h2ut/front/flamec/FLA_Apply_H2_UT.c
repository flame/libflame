
#include "FLAME.h"

FLA_Error FLA_Apply_H2_UT( FLA_Side side, FLA_Obj tau, FLA_Obj u2, FLA_Obj a1, FLA_Obj A2 )
{
  FLA_Error r_val;

  // Check parameters.
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Apply_H2_UT_check( side, tau, u2, a1, A2 );

  if ( FLA_Obj_has_zero_dim( a1 ) ) return FLA_SUCCESS;

  // Invoke FLA_Apply_H2_UT_internal() to parse parameters.
  r_val = FLA_Apply_H2_UT_internal( side, tau, u2, a1, A2 );

  return r_val;
}

