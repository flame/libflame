
#include "FLAME.h"

extern fla_scalr_t* fla_scalr_cntl_blas;

FLA_Error FLA_Scalr( FLA_Uplo uplo, FLA_Obj alpha, FLA_Obj A )
{
  FLA_Error r_val;

#ifdef FLA_ENABLE_BLAS1_FRONT_END_CNTL_TREES
  // Check parameters.
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Scalr_check( alpha, A );

  // Invoke FLA_Scalr_internal() with flat control tree that simply calls
  // external wrapper.
  r_val = FLA_Scalr_internal( uplo, alpha, A, fla_scalr_cntl_blas );

#else
  r_val = FLA_Scalr_external( uplo, alpha, A );
#endif

  return r_val;
}

