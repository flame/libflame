
#include "FLAME.h"

FLA_Error FLA_Accum_T_UT( FLA_Direct direct, FLA_Store storev, FLA_Obj A, FLA_Obj tau, FLA_Obj T )
{
  FLA_Error r_val;

  // Check parameters.
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Accum_T_UT_check( direct, storev, A, tau, T );

  // Invoke FLA_Accum_T_UT_internal().
  r_val = FLA_Accum_T_UT_internal( direct, storev, A, tau, T );

  return r_val;
}

