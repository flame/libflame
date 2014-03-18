
#include "FLAME.h"

FLA_Error FLA_Random_spd_matrix( FLA_Uplo uplo, FLA_Obj A )
{
  FLA_Obj R;

  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Random_spd_matrix_check( uplo, A );

  // Create a temporary object R conformal to A.
  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, A, &R );

  // Randomize R to be an uplo-triangular matrix. Note that the diagonal of R
  // needs to be positive to ensure that R * R' is SPD/HPD.
  FLA_Random_tri_matrix( uplo, FLA_NONUNIT_DIAG, R );
    
  if ( uplo == FLA_LOWER_TRIANGULAR )
  {
    // A = R * R';
    FLA_Herk_external( uplo, FLA_NO_TRANSPOSE, FLA_ONE, R, FLA_ZERO, A );
  }
  else // if ( uplo == FLA_UPPER_TRIANGULAR )
  {
    // A = R' * R;
    FLA_Herk_external( uplo, FLA_CONJ_TRANSPOSE, FLA_ONE, R, FLA_ZERO, A );
  }

  // Free R.
  FLA_Obj_free( &R );

  return FLA_SUCCESS;
}

