
#include "FLAME.h"

extern fla_tridiagut_t* fla_tridiagut_cntl_fused;
extern fla_tridiagut_t* fla_tridiagut_cntl_nofus;
extern fla_tridiagut_t* fla_tridiagut_cntl_plain;

FLA_Error FLA_Tridiag_UT( FLA_Uplo uplo, FLA_Obj A, FLA_Obj T )
{
  FLA_Error r_val;

  // Check parameters.
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Tridiag_UT_check( uplo, A, T );

  if ( FLA_Obj_row_stride( A ) == 1 &&
       FLA_Obj_is_double_precision( A ) )
    // Temporary fix not to use the fused version (numerically unstable).
    r_val = FLA_Tridiag_UT_internal( uplo, A, T, fla_tridiagut_cntl_plain );
  else
    r_val = FLA_Tridiag_UT_internal( uplo, A, T, fla_tridiagut_cntl_nofus );

  return r_val;
}

