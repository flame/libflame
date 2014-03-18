
#include "FLAME.h"

FLA_Error FLASH_Chol_solve( FLA_Uplo uplo, FLA_Obj A, FLA_Obj B, FLA_Obj X )
{
  // Check parameters.
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Chol_solve_check( uplo, A, B, X );

  FLASH_Copy( B, X );

  if ( uplo == FLA_LOWER_TRIANGULAR )
  {
      FLASH_Trsm( FLA_LEFT, FLA_LOWER_TRIANGULAR, FLA_NO_TRANSPOSE,
                  FLA_NONUNIT_DIAG, FLA_ONE, A, X );
      FLASH_Trsm( FLA_LEFT, FLA_LOWER_TRIANGULAR, FLA_CONJ_TRANSPOSE,
                  FLA_NONUNIT_DIAG, FLA_ONE, A, X );
  }
  else // if ( uplo == FLA_UPPER_TRIANGULAR )
  {
      FLASH_Trsm( FLA_LEFT, FLA_UPPER_TRIANGULAR, FLA_CONJ_TRANSPOSE,
                  FLA_NONUNIT_DIAG, FLA_ONE, A, X );
      FLASH_Trsm( FLA_LEFT, FLA_UPPER_TRIANGULAR, FLA_NO_TRANSPOSE,
                  FLA_NONUNIT_DIAG, FLA_ONE, A, X );
  }

  return FLA_SUCCESS;
}

