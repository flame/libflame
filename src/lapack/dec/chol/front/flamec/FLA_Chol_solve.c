
#include "FLAME.h"

FLA_Error FLA_Chol_solve( FLA_Uplo uplo, FLA_Obj A, FLA_Obj B, FLA_Obj X )
{
  // Check parameters.
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Chol_solve_check( uplo, A, B, X );

  if ( FLA_Obj_is_identical( B, X ) == FALSE )
    FLA_Copy_external( B, X );

  if ( uplo == FLA_LOWER_TRIANGULAR )
  {
      FLA_Trsm_external( FLA_LEFT, FLA_LOWER_TRIANGULAR, FLA_NO_TRANSPOSE,
                         FLA_NONUNIT_DIAG, FLA_ONE, A, X );
      FLA_Trsm_external( FLA_LEFT, FLA_LOWER_TRIANGULAR, FLA_CONJ_TRANSPOSE,
                         FLA_NONUNIT_DIAG, FLA_ONE, A, X );
  }
  else // if ( uplo == FLA_UPPER_TRIANGULAR )
  {
      FLA_Trsm_external( FLA_LEFT, FLA_UPPER_TRIANGULAR, FLA_CONJ_TRANSPOSE,
                         FLA_NONUNIT_DIAG, FLA_ONE, A, X );
      FLA_Trsm_external( FLA_LEFT, FLA_UPPER_TRIANGULAR, FLA_NO_TRANSPOSE,
                         FLA_NONUNIT_DIAG, FLA_ONE, A, X );
  }

  return FLA_SUCCESS;
}

