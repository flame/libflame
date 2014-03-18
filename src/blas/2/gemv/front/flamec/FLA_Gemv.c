
#include "FLAME.h"

extern fla_gemv_t* fla_gemv_cntl_blas;

FLA_Error FLA_Gemv( FLA_Trans transa, FLA_Obj alpha, FLA_Obj A, FLA_Obj x, FLA_Obj beta, FLA_Obj y )
{
  FLA_Error r_val;

  // Check parameters.
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Gemv_check( transa, alpha, A, x, beta, y );

#ifdef FLA_ENABLE_BLAS2_FRONT_END_CNTL_TREES
  // Invoke FLA_Gemv_internal() with flat control tree that simply calls
  // external wrapper.
  r_val = FLA_Gemv_internal( transa, alpha, A, x, beta, y, fla_gemv_cntl_blas );

#else
  r_val = FLA_Gemv_external( transa, alpha, A, x, beta, y );
#endif

  return r_val;
}

