
#include "FLAME.h"

FLA_Error FLA_Set_to_identity( FLA_Obj A )
{
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Set_to_identity_check( A );

  FLA_Set( FLA_ZERO, A );

  FLA_Set_diag( FLA_ONE, A );

  return FLA_SUCCESS;
}

