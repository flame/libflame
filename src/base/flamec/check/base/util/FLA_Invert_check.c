
#include "FLAME.h"

FLA_Error FLA_Invert_check( FLA_Conj conj, FLA_Obj x )
{
  FLA_Error e_val;

  e_val = FLA_Check_valid_conj( conj );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_floating_object( x );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_nonconstant_object( x );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_if_vector( x );
  FLA_Check_error_code( e_val );

  return FLA_SUCCESS;
}

