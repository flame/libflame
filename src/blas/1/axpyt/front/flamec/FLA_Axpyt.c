
#include "FLAME.h"

extern fla_axpyt_t* fla_axpyt_cntl_blas;

FLA_Error FLA_Axpyt( FLA_Trans trans, FLA_Obj alpha, FLA_Obj A, FLA_Obj B )
{
  FLA_Error r_val;

#ifdef FLA_ENABLE_BLAS1_FRONT_END_CNTL_TREES
  // Check parameters.
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Axpyt_check( trans, alpha, A, B );

  // Invoke FLA_Axpyt_internal() with flat control tree that simply calls
  // external wrapper.
  r_val = FLA_Axpyt_internal( trans, alpha, A, B, fla_axpyt_cntl_blas );

#else
  r_val = FLA_Axpyt_external( trans, alpha, A, B );
#endif

  return r_val;
}

