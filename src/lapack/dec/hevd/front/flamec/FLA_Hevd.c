
#include "FLAME.h"

FLA_Error FLA_Hevd( FLA_Evd_type jobz, FLA_Uplo uplo, FLA_Obj A, FLA_Obj l )
{
  FLA_Error r_val      = FLA_SUCCESS;
  dim_t     n_iter_max = 30;
  dim_t     k_accum    = 32;
  dim_t     b_alg      = 512;

  // Check parameters.
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Hevd_check( jobz, uplo, A, l );

  // Invoke FLA_Hevd_external() for now.
  if ( jobz == FLA_EVD_WITH_VECTORS )
  {
    if ( uplo == FLA_LOWER_TRIANGULAR )
    {
      r_val = FLA_Hevd_lv_unb_var1( n_iter_max, A, l, k_accum, b_alg );
    }
    else // if ( uplo == FLA_UPPER_TRIANGULAR )
    {
      FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );
    }
  }
  else // if ( jobz == FLA_EVD_WITHOUT_VECTORS )
  {
    if ( uplo == FLA_LOWER_TRIANGULAR )
    {
      FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );
    }
    else // if ( uplo == FLA_UPPER_TRIANGULAR )
    {
      FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );
    }
  }

  return r_val;
}

