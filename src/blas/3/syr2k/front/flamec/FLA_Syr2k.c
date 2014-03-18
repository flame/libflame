
#include "FLAME.h"

extern fla_syr2k_t* fla_syr2k_cntl_mm;

FLA_Error FLA_Syr2k( FLA_Uplo uplo, FLA_Trans trans, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C )
{
  FLA_Error r_val = FLA_SUCCESS;
  
  // Check parameters.
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Syr2k_check( uplo, trans, alpha, A, B, beta, C );

#ifdef FLA_ENABLE_BLAS3_FRONT_END_CNTL_TREES
  r_val = FLA_Syr2k_internal( uplo, trans, alpha, A, B, beta, C, fla_syr2k_cntl_mm );
#else
  r_val = FLA_Syr2k_external( uplo, trans, alpha, A, B, beta, C );
#endif

  return r_val;
}

