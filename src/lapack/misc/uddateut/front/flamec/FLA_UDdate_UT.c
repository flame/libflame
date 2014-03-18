
#include "FLAME.h"

extern fla_uddateut_t* fla_uddateut_cntl_leaf;
extern fla_apqudut_t*  fla_apqudut_cntl_leaf;

FLA_Error FLA_UDdate_UT( FLA_Obj R, FLA_Obj C, FLA_Obj D, FLA_Obj T )
{
  FLA_Error r_val;

  // Check parameters.
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_UDdate_UT_check( R, C, D, T );

  // Invoke the _internal() back-end with the standard control tree.
  r_val = FLA_UDdate_UT_internal( R, C, D, T, fla_uddateut_cntl_leaf );

  return r_val;
}

