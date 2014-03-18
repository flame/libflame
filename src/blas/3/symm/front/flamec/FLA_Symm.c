
#include "FLAME.h"

extern fla_symm_t* fla_symm_cntl_mm;

FLA_Error FLA_Symm( FLA_Side side, FLA_Uplo uplo, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C )
{
  FLA_Error r_val = FLA_SUCCESS;

  // Check parameters.
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Symm_check( side, uplo, alpha, A, B, beta, C );

#ifdef FLA_ENABLE_BLAS3_FRONT_END_CNTL_TREES
  r_val = FLA_Symm_internal( side, uplo, alpha, A, B, beta, C, fla_symm_cntl_mm );
#else
  r_val = FLA_Symm_external( side, uplo, alpha, A, B, beta, C );
#endif

  return r_val;
}

