
#include "FLAME.h"

extern fla_trsm_t* fla_trsm_cntl_mm;

FLA_Error FLA_Trsm( FLA_Side side, FLA_Uplo uplo, FLA_Trans trans, FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B )
{
  FLA_Error    r_val = FLA_SUCCESS;

  // Check parameters.
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Trsm_check( side, uplo, trans, diag, alpha, A, B );

#ifdef FLA_ENABLE_BLAS3_FRONT_END_CNTL_TREES
  r_val = FLA_Trsm_internal( side, uplo, trans, diag, alpha, A, B, fla_trsm_cntl_mm );
#else
  r_val = FLA_Trsm_external( side, uplo, trans, diag, alpha, A, B );
#endif

  return r_val;
}

