
#include "FLAME.h"

extern fla_trsv_t* fla_trsv_cntl_blas;

FLA_Error FLA_Trsv( FLA_Uplo uplo, FLA_Trans trans, FLA_Diag diag, FLA_Obj A, FLA_Obj x )
{
  FLA_Error r_val;

  // Check parameters.
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Trsv_check( uplo, trans, diag, A, x );

#ifdef FLA_ENABLE_BLAS2_FRONT_END_CNTL_TREES
  // Invoke FLA_Trsv_internal() with flat control tree that simply calls
  // external wrapper.
  r_val = FLA_Trsv_internal( uplo, trans, diag, A, x, fla_trsv_cntl_blas );

#else
  r_val = FLA_Trsv_external( uplo, trans, diag, A, x );
#endif

  return r_val;
}

