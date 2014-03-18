
#include "FLAME.h"

FLA_Error FLASH_LU_nopiv_solve( FLA_Obj A, FLA_Obj B, FLA_Obj X )
{
  // Check parameters.
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_LU_nopiv_solve_check( A, B, X );

  FLASH_Copy( B, X );

  FLASH_Trsm( FLA_LEFT, FLA_LOWER_TRIANGULAR, FLA_NO_TRANSPOSE,
              FLA_UNIT_DIAG, FLA_ONE, A, X );
  FLASH_Trsm( FLA_LEFT, FLA_UPPER_TRIANGULAR, FLA_NO_TRANSPOSE,
              FLA_NONUNIT_DIAG, FLA_ONE, A, X );

  return FLA_SUCCESS;
}


