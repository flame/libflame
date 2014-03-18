
#include "FLAME.h"

extern fla_spdinv_t* fla_spdinv_cntl;

FLA_Error FLA_SPDinv( FLA_Uplo uplo, FLA_Obj A )
{
  FLA_Error r_val = FLA_SUCCESS;

  // Check parameters.
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_SPDinv_check( uplo, A );

  // Invoke FLA_SPDinv_internal() with an appropriate control tree.
  r_val = FLA_SPDinv_internal( uplo, A, fla_spdinv_cntl );

  return r_val;
}

