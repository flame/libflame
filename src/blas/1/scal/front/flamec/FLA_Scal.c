
#include "FLAME.h"

extern fla_scal_t* fla_scal_cntl_blas;

FLA_Error FLA_Scal( FLA_Obj alpha, FLA_Obj A )
{
  FLA_Error r_val;

#ifdef FLA_ENABLE_BLAS1_FRONT_END_CNTL_TREES
  // Check parameters.
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Scal_check( alpha, A );

  // Invoke FLA_Scal_internal() with flat control tree that simply calls
  // external wrapper.
  r_val = FLA_Scal_internal( alpha, A, fla_scal_cntl_blas );

#else
  r_val = FLA_Scal_external( alpha, A );
#endif

  return r_val;
}

