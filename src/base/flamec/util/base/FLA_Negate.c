
#include "FLAME.h"

FLA_Error FLA_Negate( FLA_Obj x )
{
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Negate_check( x );

  return FLA_Scal( FLA_MINUS_ONE, x );
}

