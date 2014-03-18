
#include "FLAME.h"

FLA_Error FLA_Tridiag_UT_scale_diagonals( FLA_Uplo uplo, FLA_Obj alpha, FLA_Obj A )
{
  FLA_Error r_val = FLA_SUCCESS;

  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Tridiag_UT_scale_diagonals_check( uplo, alpha, A );
  
  if ( uplo == FLA_LOWER_TRIANGULAR )
    r_val = FLA_Bidiag_UT_l_scale_diagonals( alpha, A );
  else
    r_val = FLA_Bidiag_UT_u_scale_diagonals( alpha, A );

  return r_val;
}
