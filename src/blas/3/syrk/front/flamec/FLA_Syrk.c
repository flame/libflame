
#include "FLAME.h"

extern fla_syrk_t* fla_syrk_cntl_mm;

FLA_Error FLA_Syrk( FLA_Uplo uplo, FLA_Trans trans, FLA_Obj alpha, FLA_Obj A, FLA_Obj beta, FLA_Obj C )
{
  FLA_Error r_val = FLA_SUCCESS;

  // Check parameters.
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Syrk_check( uplo, trans, alpha, A, beta, C );

#ifdef FLA_ENABLE_BLAS3_FRONT_END_CNTL_TREES
  r_val = FLA_Syrk_internal( uplo, trans, alpha, A, beta, C, fla_syrk_cntl_mm );
#else
  r_val = FLA_Syrk_external( uplo, trans, alpha, A, beta, C );
#endif
  
  return r_val;
}

