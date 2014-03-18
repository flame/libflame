
#include "FLAME.h"

FLA_Error FLA_LU_piv_solve( FLA_Obj A, FLA_Obj p, FLA_Obj B, FLA_Obj X )
{
  // Check parameters.
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_LU_piv_solve_check( A, p, B, X );

  if ( FLA_Obj_is_identical( B, X ) == FALSE ) 
    FLA_Copy_external( B, X );

  FLA_Apply_pivots( FLA_LEFT, FLA_NO_TRANSPOSE, p, X );

  FLA_Trsm_external( FLA_LEFT, FLA_LOWER_TRIANGULAR, FLA_NO_TRANSPOSE,
                     FLA_UNIT_DIAG, FLA_ONE, A, X );
  FLA_Trsm_external( FLA_LEFT, FLA_UPPER_TRIANGULAR, FLA_NO_TRANSPOSE,
                     FLA_NONUNIT_DIAG, FLA_ONE, A, X );

  return FLA_SUCCESS;
}

