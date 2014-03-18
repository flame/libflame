
#include "FLAME.h"

extern fla_lyap_t* fla_lyap_cntl;

FLA_Error FLA_Lyap( FLA_Trans trans, FLA_Obj isgn, FLA_Obj A, FLA_Obj C, FLA_Obj scale )
{
  FLA_Error r_val;

  // Check parameters.
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Lyap_check( trans, isgn, A, C, scale );

  // Invoke FLA_Lyap_internal() with the appropriate control tree.
  r_val = FLA_Lyap_internal( trans, isgn, A, C, scale, fla_lyap_cntl );

  return r_val;
}

